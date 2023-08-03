#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import healpy as hp
from scipy.stats import binned_statistic
import rubin_sim.maf as maf
from rubin_sim.data import get_baseline
import time
import sys
import sqlite3
import json
import glob


# In[2]:


# Built on Peter and Lynne's original notebooks:
# https://github.com/rhiannonlynne/notebooks/blob/master/Template%20Generation.ipynb
# https://github.com/yoachim/22_Scratch/blob/main/template_metrics/template_stuff.ipynb


# # Set up variables and get the first year of observations

# In[3]:


# select the survey simulation
baseline_db = "baseline_v3.0_10yrs.db"
year1_fname = 'first_year_{}.db'.format(baseline_db.split(".db")[0])

# range of nights to investigate (year 1 of survey)
night_min = 0
night_max = 365


# In[4]:


# set nside for analysis metric (use higher resolution when removing no-template images later)
# nside = 32
nside = 256
s_hp = maf.HealpixSlicer(nside=nside)

# make a time slicer
nights = np.arange(night_min, night_max+1, 1)
s_t = maf.OneDSlicer(slice_col_name="night", bins = nights-0.5) # subtract half a night as slicer bins are left and right inclusive


# In[5]:


# set path for metric output files
save_dir = "remove_no_template_results_{}_{}".format(nside,year1_fname.replace(".","_"))
print(save_dir)

# In[6]:


# use the healpix area to find approximate number of healpixels in a single visit
pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) # square degrees
lsst_footprint = 9.6 # square degrees, GET EXACT NUMBER?
n_pix = lsst_footprint/pix_area # number of healpixels within each visit footprint
print("healpixel area = {} square degrees\n number of healpixels in visit = {}".format(pix_area,n_pix))


# In[7]:


# Set some fraction of visit covered by template required to be analysed
frac = 0.9


# In[8]:


# Choose to ignore deep drilling fields
# This string will be applied to metric sql query and also to the redacted database code
# sqlDD = ''
sqlDD = ' and note not like "%DD%"'

if sqlDD!='':
    save_dir+="_noDD" # change the save directory

print(save_dir)


# In[9]:


# query or load the database of year 1 observations
if not os.path.isfile(year1_fname):
    print("get year 1 observations")
    conn = sqlite3.connect(baseline_db)
    df = pd.read_sql('select * from observations;', conn)
    df_year1 = df[df["night"]<=night_max]
    conn.close()

    # open up a connection to a new database
    conn = sqlite3.connect(year1_fname)
    df_year1.to_sql('observations', conn, index=False, if_exists='replace')
    conn.close()

else:
    conn = sqlite3.connect(year1_fname)
    df_year1 = pd.read_sql('select * from observations;', conn)
    conn.close()


# In[10]:


# set up some filenames
opsim_fname = year1_fname
print(opsim_fname)
# opsdb = maf.OpsimDatabase(opsim_fname)
opsdb = opsim_fname
runName = os.path.split(opsdb)[-1].replace('.db', '')
print(runName)


# In[11]:


# template generation timescales to test
# tscales = [28,28*4,28*6]
tscales = [7]

# store the timescale and template generation nights in a dict
template_timescales = {}

for tscale in tscales:

    # divide year 1 into chunks of a given template_timescale
    template_nights = np.arange(0,night_max+tscale,tscale)
    template_nights[-1] = night_max # consider only the first year

    template_timescales[str(tscale)] = template_nights

template_timescales


# # Remove visits without templates using Peter's code
#
# Add Lynne's requirements on image quality when counting number of available images to generate template with.
#
# Peter's code assumes that a template will be generated once, on a given night.

# In[41]:


# new base class. We might consider moving this into BaseMetric, or maybe doing this in the slicer instead.
class BaseTemplateMetric(maf.metrics.BaseMetric):
    """Make a new base class that will filter out observation from before we have templates available
    """
    def __init__(self,
                col=None,
                n_visits_for_template=3., # note that this is float
                night_template_min=182,
                mjd_col='observationStartMJD',
                night_col='night',
                filter_col='filter',
                 seeing_ratio=2.0, m5_range=0.5,
                 seeing_col='seeingFwhmEff', m5_col='fiveSigmaDepth',
                 **kwargs):
        if col is None:
            col = []
        else:
            col = [col]
        col += [mjd_col, night_col, filter_col]
        super().__init__(col=col, **kwargs)
        self.n_visits_for_template = n_visits_for_template
        self.night_template_min = night_template_min
        self.mjd_col = mjd_col
        self.night_col = night_col
        self.filter_col = filter_col

        self.seeing_col = seeing_col
        self.m5_col = m5_col
        self.seeing_ratio = seeing_ratio
        self.m5_range = m5_range

        # this snippet (from Lynne nb) is required to set the columns that need to be called when running the metric
        if 'metric_name' in kwargs:
            self.metric_name = kwargs['metric_name']
            del kwargs['metric_name']
        else:
            self.metric_name = 'BaseTemplateMetric'
        super().__init__(col=[self.mjd_col, self.seeing_col, self.m5_col, self.night_col, self.filter_col],
                         metric_name=self.metric_name, **kwargs)

        print("night_template_min = {}, seeing_ratio = {}, m5_range = {}".format(self.night_template_min,self.seeing_ratio,self.m5_range))

    def _remove_no_template_visits(self, dataSlice):

        dataSlice.sort(order=self.mjd_col)

        n_visits = int(self.n_visits_for_template) # number of good images required for template

        # Apply the masks to remove images not suitable for templates/science
        # everything starts as true and is switched to false if required
        has_template_indx = np.ones(dataSlice.size, dtype=bool) # science images with templates
        template_img_indx = np.ones(dataSlice.size, dtype=bool) # images used for templates

        # template images must have been taken before the generation date
        template_time = np.where(dataSlice[self.night_col] < self.night_template_min,
                        True, False)
        # science images can only be counted after the generation date
        image_time = np.where(dataSlice[self.night_col] > self.night_template_min,
                        True, False)
        ### should use a greater than sign for nights here? Or allow leeway for template generation time?

        # define list of possible science/template images by removing impossible images
        has_template_indx[~(image_time)] = False
        template_img_indx[~(template_time)] = False

#         print(self.night_template_min,len(dataSlice),sum(template_img_indx),sum(has_template_indx))

        # look at each dataslice by filter, because template images must be in the same filter
        for filtername in np.unique(dataSlice[self.filter_col]):

            # mask for each filter
            infilt = np.where(dataSlice[self.filter_col] == filtername)[0] # index zero required to get just the mask (or use True, False arguments like Lynne above)
            # mask for template images in the filter
            infilt_templates = infilt & template_time[infilt]

            # What if there are zero template images in this filter?
            if sum(infilt_templates)==0:
                has_template_indx[infilt] = False # there can be no science images
                continue

            # Find the best seeing and depth in the available template images
            bench_seeing = np.min(dataSlice[infilt_templates][self.seeing_col])
            bench_m5 = np.max(dataSlice[infilt_templates][self.m5_col])

            # define the masks for template visits in this filter meeting the seeing/depth requirements
            seeing_ok = np.where(dataSlice[infilt_templates][self.seeing_col]/bench_seeing < self.seeing_ratio,
                                True, False)
            m5_ok = np.where(bench_m5 - dataSlice[infilt_templates][self.m5_col] < self.m5_range,
                            True, False)

            # define list of images that are suitable for generating templates by removing "bad" visits
            # N.B. we have accounted for the template generation night above
            template_img_indx[infilt][~(seeing_ok & m5_ok)] = False

            # check if there are is a sufficent number of images for making the templates
            # if not, no science images can be counted in this healpixel (until next template generation time)
            if sum(template_img_indx[infilt])<n_visits:
                has_template_indx[infilt] = False

#         print(sum(has_template_indx))

        # return the mask of science images with templates and the mask of template images
        return has_template_indx, template_img_indx


# In[42]:



# OK, now let's write a function that takes a bunch of observations in,
# then cuts out ones that probably don't have template images

# Class that will return all the id's that we think have templates generated
class HasTemplateIndx(BaseTemplateMetric):
    """Return the obIds that probably have templates
    This then allows you to run metrics directly on the visits of the dataslice
    that had templates and could generate alerts
    """
    def __init__(self, col='observationId', metric_dtype="object", **kwargs):
        super().__init__(col=col, metric_dtype=metric_dtype, **kwargs)
        self.idCol = col
    def run(self, dataSlice, slice_point=None):
        has_template_indx, template_img_indx = self._remove_no_template_visits(dataSlice) # get the mask of science images
        dataSlice = dataSlice[has_template_indx] # apply the mask to the dataSlice
        return dataSlice[self.idCol]


# In[52]:



def remove_no_templates(data_in, nside=32, id_col='observationId',
                        night_template_min = 182,
                        template_col = "npix_template", # column name to store fraction of visit covered by template
                        **kwargs):
    """
    This function removes all visits from data_in that had zero template coverage.
    data_in is some chunk of year 1 with visits between times of template generation.
    Visits that had a template for at least one of its constituent healpixels are retained.
    We add a new column that records the number of healpixels with templates within each visit.
    """

    slicer = maf.slicers.HealpixSlicer(nside=nside, verbose=False)
    metric = HasTemplateIndx(night_template_min = night_template_min, **kwargs)
    print("night_template_min = {}".format(metric.night_template_min))
    sql=None
    bundle = maf.MetricBundle(metric, slicer, sql)
    mbg = maf.MetricBundleGroup([bundle], None, save_early=False,
                                verbose=False)
    mbg.run_current(None, sim_data=data_in)

    # we retrieve an array with the id_col (observationId) of every single visit
    # that had at least one constituent healpixel with a template
    # An id_col value can appear multiple times as there are multiple healpixels within each visit
    all_vals = np.concatenate(bundle.metric_values.data[~bundle.metric_values.mask])

    # we find all unique id_col, these are the visits that had some form of template
    # we count the number of each unique id_col, this is the number of healpixels within the visit that had templates
    valid_ids,count_ids = np.unique(all_vals,return_counts=True)
#     print(valid_ids,count_ids)

    # keep only visits with templates
    indx = np.in1d(data_in[id_col], valid_ids)
    result = data_in[indx]
    ### ALSO KEEP VISITS WHERE TEMPLATE_COL = 0?

    # add a column to track number of healpixels with templates within visits
    # ONLY for visits that had at least one template
    result[template_col] = count_ids

    # return the database of visits with templates
    return result


# # Run metrics on dataslice directly

# In[44]:


class doAllTemplateMetrics(BaseTemplateMetric):

    # run all the metrics on dataSlice so that the the removal of visits only occurs once

    def __init__(self, col='observationId', **kwargs):
        # ensure all kwargs are passed to teh template class (e.g. night_template_min)
        super().__init__(col=col, **kwargs)
        self.idCol = col

    def run(self, dataSlice, slice_point=None):

        # run the template checking code to create a dataSlice with only the visits with templates
        has_template_indx, template_img_indx = self._remove_no_template_visits(dataSlice) # get the mask of science images

        # return a list of both dataSlices, one with science images, one with template images
        # Also pass the dataslice of first night of all the visits in that healpixel
        return [dataSlice[has_template_indx],
                dataSlice[template_img_indx],
                dataSlice[0]] # This list is stored as metricVal

    def reduceCount(self,metricVal):
        """
        Count the number of visits with existing templates within each healpix
        When grouping take the sum
        """
        sci_img = metricVal[0] # science images
        count = len(sci_img)
        if count==0:
            count = self.badval # set badval to help with metric plotting later
        return count

### Move this one to a separate metric as we want to constrain the pairs of filters
#     def reducePairs(self,metricVal):
#         # Get the number of pairs available in each healpix
#         # When grouping take the sum

#         x = metricVal[0] # science images
#         if len(x)==0:
#             pairs = self.badval
#         else:
#             pairs = maf.metrics.pairMetric.PairMetric(mjdCol='observationStartMJD', metricName='Pairs').run(x)
#         return pairs

    def reduceNight(self,metricVal):
        """
        Get the first possible night when templates are available
        Note that this should not count beyond the template generation night
        When grouping select the minimum value
        """

        sci_img = metricVal[0] # science images
        temp_img = metricVal[1] # dataSlice of template images

        # N.B. there is a subtlety between counting when the first science images appear and when the template is first available
#         if (len(sci_img)==0):
        if (len(temp_img)<int(self.n_visits_for_template)): # only count if there are enough images

            night = self.badval
        else:
            # if there is at least one healpixel which satisifies templates
            # then we can say that the template was made on at least the template generation night
            night = self.night_template_min
        return night

    def reduceDeltaNight(self, metricVal):
        """
        This metric returns time since the first visit for the first template to be available.
        Different to Lynne's implementation which did not account for filters.

        When grouping select the first non-badval number (or minimum?)
        """

        sci_img = metricVal[0] # dataSlice of science images
        temp_img = metricVal[1] # dataSlice of template images
        first_night = metricVal[2][self.night_col] # first night of visit

#         # This version of the metric ignores the template generation timescale, to account for distinct nights:
#         # metric_plot = np.ceil(metric_plot/template_timescale)*template_timescale
#         if (len(temp_img)>=int(self.n_visits_for_template)):
#             dn = temp_img[self.nightCol][int(self.n_visits_for_template)-1] - first_night
#         else:
#             dn = self.badval

        # This version of the metric rounds (up) to the template generation night
        # i.e. the night of possible template generation (temp_img[self.nightCol][int(self.n_visits_for_template)-1])
        # is probably less then night_template_min
        if (len(temp_img)>=int(self.n_visits_for_template)):
#         if (len(sci_img)>0): # science images exist
            dn = self.night_template_min - first_night
#             dt = self.night_template_min - temp_img[self.nightCol][0] # optional, get time range of just the template images, basically the same as above

        else:
            dn = self.badval

        return dn

    def reduceNTemplate(self,metricVal):
        """
        Returns the total number of images used to make a template.
        When grouping results, take the first non-badval number (or minimum?)
        """

        sci_img = metricVal[0] # dataSlice of science images
        temp_img = metricVal[1] # dataSlice of template images
        if (len(temp_img)>=int(self.n_visits_for_template)): # only count if there are enough template images
#         if (len(sci_img)>0): # science images exist
            NTemplate = len(temp_img) # count total number of templates
        else:
            NTemplate = self.badval
        return NTemplate

    def reduceSeeingTemplate(self,metricVal):
        """
        Returns the mean seeing of images used to make a template.
        When grouping results, take the mean (accounting for badvals?)
        """

        sci_img = metricVal[0] # dataSlice of science images
        temp_img = metricVal[1] # dataSlice of template images
        if (len(temp_img)<int(self.n_visits_for_template)): # only count if there are enough images
#         if (len(sci_img)==0): # no science images exist
            seeing = self.badval
        else:
            seeing = np.mean(temp_img[self.seeing_col]) # calculate the mean seeing
        return seeing

    def reduceDepthTemplate(self,metricVal):
        """
        Returns the mean depth of images used to make a template.
        When grouping results, take the mean (accounting for badvals?)
        """

        sci_img = metricVal[0] # dataSlice of science images
        temp_img = metricVal[1] # dataSlice of template images
        if (len(temp_img)<int(self.n_visits_for_template)): # only count if there are enough images
#         if (len(sci_img)==0): # no science images exist
            depth = self.badval
        else:
            depth = np.mean(temp_img[self.m5_col]) # calculate the mean depth
        return depth

    ### additional metrics
    # date range of first n_visits_for_template template images?
    # As we increase template generation timescale this will become the same as reduceDeltaNight


# In[45]:


class doPairTemplateMetrics(BaseTemplateMetric):

    # run all the metrics on dataSlice so that the the removal of visits only occurs once

    def __init__(self, col='observationId', **kwargs):
        # ensure all kwargs are passed to teh template class (e.g. night_template_min)
        super().__init__(col=col, **kwargs)
        self.idCol = col

    def run(self, dataSlice, slice_point=None):

        # run the template checking code to create a dataSlice with only the visits with templates
        has_template_indx, template_img_indx = self._remove_no_template_visits(dataSlice) # get the mask of science images

        # return a list of both dataSlices, one with science images, one with template images
        # Also pass the first night of all the visits
        return [dataSlice[has_template_indx],
                dataSlice[template_img_indx],
                dataSlice[0]] # This list is stored as metricVal

    def reducePairs(self,metricVal):
        # Get the number of pairs available in each healpix, for any combination of the queried filters
        # When grouping take the sum

        x = metricVal[0] # science images
        if len(x)==0:
            pairs = self.badval
        else:
            pairs = maf.metrics.PairMetric(mjd_col='observationStartMJD', metric_name='Pairs').run(x)
        return pairs


# In[46]:


# Do regular metrics on the baseline, count visits and pairs

start = time.time()

slicer1 = s_hp # spatial slicer
slicer2 = s_t # time (night) slicer

t_data=night_max

_runName = "{}_nside-{}".format(runName,nside)

bl = []
summary_stats = [maf.MedianMetric()]

for filt in ["all","u","g","r","i","z","y",
                        "r_or_g","r_or_i","r_or_z"]:

# for filt in ["all"]:

    if filt=="all":
        sql = 'night <= {}'.format(t_data)+sqlDD
    elif "_or_" in filt:
        sql = '(filter = "r" or filter = "{}") and night <= {}'.format(filt.split("_or_")[-1],t_data)+sqlDD
    else:
        sql = 'filter = "{}" and night <= {}'.format(filt,t_data)+sqlDD

    # Run the regular metrics without templates
    metric = maf.CountMetric(col='night', metric_name = "CountMetric")
    bl.append(maf.MetricBundle(metric, slicer1, sql, summary_metrics=summary_stats, run_name=_runName))

    metric = maf.metrics.PairMetric(mjd_col='observationStartMJD', metric_name = "PairMetric")
    bl.append(maf.MetricBundle(metric, slicer1, sql, summary_metrics=summary_stats, run_name=_runName))

    # Run the time based metrics
    metric = maf.CountMetric(col='night', metric_name = "CountMetric")
    bl.append(maf.MetricBundle(metric, slicer2, sql, summary_metrics=summary_stats, run_name=_runName))

    metric = maf.metrics.PairMetric(mjd_col='observationStartMJD', metric_name = "PairMetric")
    bl.append(maf.MetricBundle(metric, slicer2, sql, summary_metrics=summary_stats, run_name=_runName))

mg = maf.MetricBundleGroup(bl, opsdb, out_dir=save_dir)

mg.run_all()
mg.plot_all(closefigs=False)

end = time.time()
dt = end-start
print("{}s ({:2f}min or {:2f}hrs)".format(dt,dt/60,dt/60/60))


# In[18]:


# do metrics with templates for a spatial slicer

start = time.time()

slicer = s_hp

for template_timescale in tscales:

    _runName = "{}_tscale-{}_nside-{}".format(runName,template_timescale,nside)

    template_nights = template_timescales[str(template_timescale)]

    for i in range(1,len(template_nights)):

        t_data = template_nights[i] # query all visits up to the time when we must consider new template generation
        t_template = template_nights[i-1] # this is the last date at which templates were generated

        print(template_timescale,t_template,t_data)

        bl = []
        summary_stats = [maf.MedianMetric()]

        for filt in ["all","u","g","r","i","z","y"]:

            if filt=="all":
                sql = 'night <= {}'.format(t_data)+sqlDD
            else:
                sql = 'filter = "{}" and night <= {}'.format(filt,t_data)+sqlDD

            # Our new metric that only counts things after templates have been generated
            metric = doAllTemplateMetrics(units='count', metric_name = "doAllTemplateMetrics",
                                             night_template_min = t_template)
            bl.append(maf.MetricBundle(metric, slicer, sql, summary_metrics=summary_stats, run_name=_runName))

        mg = maf.MetricBundleGroup(bl, opsdb, out_dir=save_dir)

        mg.run_all()
#             mg.plotAll(closefigs=False)

end = time.time()
dt = end-start
print("{}s ({:2f}min or {:2f}hrs)".format(dt,dt/60,dt/60/60))


# In[19]:


# run pair metrics with templates for a spatial slicer
# use all filters or filter pairs: r-g, r-i, r-z

start = time.time()

slicer = s_hp

for template_timescale in tscales:

    _runName = "{}_tscale-{}_nside-{}".format(runName,template_timescale,nside)

    template_nights = template_timescales[str(template_timescale)]

    for i in range(1,len(template_nights)):

        t_data = template_nights[i] # query all visits up to the time when we must consider new template generation
        t_template = template_nights[i-1] # this is the last date at which templates were generated

        print(template_timescale,t_template,t_data)

        bl = []
        summary_stats = [maf.MedianMetric()]

        for filt in ["all","g","i","z"]:

            if filt=="all":
                sql = 'night <= {}'.format(t_data)+sqlDD
            else:
                sql = '(filter = "r" or filter = "{}") and night <= {}'.format(filt,t_data)+sqlDD

            # Our new metric that only counts things after templates have been generated
            metric = doPairTemplateMetrics(units='count', metric_name = "doPairTemplateMetrics",
                                             night_template_min = t_template)
            bl.append(maf.MetricBundle(metric, slicer, sql, summary_metrics=summary_stats, run_name=_runName))

        mg = maf.MetricBundleGroup(bl, opsdb, out_dir=save_dir)

        mg.run_all()
#             mg.plotAll(closefigs=False)

end = time.time()
dt = end-start
print("{}s ({:2f}min or {:2f}hrs)".format(dt,dt/60,dt/60/60))



# # Run the full remove template code on all date ranges

# In[53]:


template_col = "npix_template"

start = time.time()

print(opsdb,nside)

_runName = "{}_nside-{}".format(runName,nside).replace(".","_")
templates_fname = "{}/{}_tscales-{}.json".format(save_dir,_runName,"-".join(str(x) for x in tscales))

if os.path.isfile(templates_fname):
    print("load {}".format(templates_fname))
    template_visits = json.load( open(templates_fname) )
else:
    template_visits = {}
    for template_timescale in tscales:

        n_visits = []
        obsIds = []
        template_frac_list = []

        template_nights = template_timescales[str(template_timescale)]

        for i in range(1,len(template_nights)):
            t_data = template_nights[i] # query all visits up to the time when we must consider new template generation
            t_template = template_nights[i-1] # this is the last date at which templates were generated
            # select visits in chunk from original year 1 database
            data = maf.get_sim_data(opsdb, None, None,
                                  full_sql_query='select * from observations where night <= {}{};'.format(t_data,sqlDD))

            # convert to dataframe
            df = pd.DataFrame(data)

            # add new column
            df[template_col] = np.nan

            # write the dataframe to recarray
            data = df.to_records(index=False)

            # remove templates without visits
            data_w_templates = remove_no_templates(data, night_template_min = t_template, nside = nside,
                                                  template_col = template_col)
            n = data_w_templates.size
            n_visits.append(n)
            obsIds += list(data_w_templates["observationId"])
            print(t_data,t_template,n,data.size,len(obsIds))

            data_w_templates = remove_no_templates(data, night_template_min = t_template, nside = nside)
            template_frac_list+=list(data_w_templates[template_col])

        # store the obsids and template_fraction of each timescale
        template_visits[str(template_timescale)] = {"n_visits":np.array(n_visits).tolist(),
                                              "template_nights":np.array(template_nights[1:]).tolist(),
                                              "npix_template":np.array(template_frac_list).tolist(),
                                                    "obsIds":np.array(obsIds).tolist()}

        print(template_frac_list[:50])
        print(sum(n_visits),len(template_frac_list))
#         break

    end = time.time()
    dt = end-start
    print("{}s ({:2f}min or {:2f}hrs)".format(dt,dt/60,dt/60/60))

    json.dump(template_visits, open(templates_fname, 'w' ) )


# In[54]:


# save the redacted database for each timescale
# select visits in chunk from original year 1 database
df_data = df_year1.copy()
print(len(df_data))

for template_timescale in tscales:

    df_data_w_templates = df_data[np.isin(df_data["observationId"],template_visits[str(template_timescale)]["obsIds"])]
    print(template_timescale, len(df_data_w_templates))

    # Insert a column for fraction template coverage
    # keep only visits with a template fractional coverage greater than some number (>0.9?)
    df_data_w_templates.loc[:,"npix_template"] = template_visits[str(template_timescale)]["npix_template"]

    fname = "{}/{}_visit_cut_t-{}d_nside-{}.db".format(save_dir,runName.replace(".","_"),template_timescale,nside)
    print(fname)

    # open up a connection to a new database
    conn = sqlite3.connect(fname)
    # save reduced visit dataframe to sql
    df_data_w_templates.to_sql('observations', conn, index=False, if_exists='replace')
    conn.close()
