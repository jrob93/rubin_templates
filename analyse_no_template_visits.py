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
from natsort import natsorted

from rubin_templates.analyse_no_template_visits_funcs import load_template_metric,skymap_plot_Night,skymap_plot,histogram_plot

from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-d", "--database", dest="baseline_db", help="baseline database file", metavar="DATABASE")
parser.add_option( "-n", "--nside", dest="nside", help="nside for healpix", metavar="NSIDE" )
parser.add_option( "-t", "--tscales", action="append", dest="tscales", help="timescales for template generation (list append)", metavar="TSCALES" )

(options,args)=parser.parse_args()

if options.baseline_db:
    baseline_db=options.baseline_db
else:
    baseline_db="baseline_v3.3_10yrs.db"
if options.nside:
    nside = int(options.nside)
else:
    nside = 256
if options.tscales:
    tscales = options.tscales
else:
    tscales = [3,7,14,28]


# In[2]:


# taken from Peter and Lynne's original notebooks
# https://github.com/rhiannonlynne/notebooks/blob/master/Template%20Generation.ipynb
# https://github.com/yoachim/22_Scratch/blob/main/template_metrics/template_stuff.ipynb


# # get just the first year of observations

# In[3]:


# baseline_db = "one_snap_v4.0_10yrs.db"
year1_fname = 'first_year_{}.db'.format(baseline_db.split(".db")[0])

# nside = 256
# nside = 32

all_filters = ["all","u","g","r","i","z","y",
                        "r_or_g","r_or_i","r_or_z"]

# number of visits required for templates
# n_visits = 3
n_visits = 4


# In[4]:


metric_bundle_savedir = "analysis_metric_bundles"

if n_visits==4:
    metric_bundle_savedir += "_n_visits_4"
    

# In[5]:


# use the healpix area to find approximate number of healpixels in a single visit
pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) # square degrees
lsst_footprint = 9.6 # square degrees, GET EXACT NUMBER?
n_pix = lsst_footprint/pix_area # number of healpixels within each visit footprint
print("healpixel area = {} square degrees\n number of healpixels in visit = {}".format(pix_area,n_pix))


# In[6]:


opsim_fname = year1_fname

# In[7]:


# opsdb = maf.OpsimDatabase(opsim_fname)
opsdb = opsim_fname
runName = os.path.split(opsdb)[-1].replace('.db', '').replace(".","_")

# In[8]:


# range of nights to investigate (year 1 of survey)
night_min = 0
night_max = 365
nights = np.arange(night_min, night_max+1, 1)

# template generation timescales to test
# tscales = [3,7,14,28]
# tscales = [7,14,28]
# tscales = [28]

# store the timescale and template generation nights in a dict
template_timescales = {}

for tscale in tscales:
    
    # divide year 1 into chunks of a given template_timescale
    template_nights = np.arange(0,night_max+tscale,tscale)
    template_nights[-1] = night_max # consider only the first year
        
    template_timescales[str(tscale)] = template_nights
    

# In[9]:


# sql = ""
# sql_DD = ""
# save_dir = "remove_no_template_results_{}_{}".format(nside, year1_fname.replace(".","_"))

# sql = ""
# sqlDD = ' and note not like "%DD%"'
# save_dir = "remove_no_template_results_{}_{}_noDD".format(nside, year1_fname.replace(".","_"))

sql = "_noDD_noTwi"
sqlDD = ' and note not like "%DD%" and note not like "%twilight%"'
save_dir = "remove_no_template_results_{}_{}".format(nside, year1_fname.replace(".","_"))

if n_visits==4:
#     save_dir = save_dir.split("_noDD_noTwi")[0]+"_n_visits_4_noDD_noTwi"
    save_dir = save_dir.split("_noDD_noTwi")[0]+"_n_visits_4"
    
print(save_dir)

# Ignore deep drilling fields
# This string will be applied to metric sql query and also to the redacted database code

if sqlDD!='':
    save_dir+=sql # change the save directory
    print(year1_fname)
    year1_fname = "{}{}.db".format(year1_fname.replace(".","_"),sql) # update the baseline db file
    print(year1_fname)
print(save_dir)


# In[10]:


# sometimes runName needs extra letters
runName = runName + "_db" + sql

# if n_visits==3 and baseline_db=="baseline_v3.3_10yrs.db":
#     runName = runName.split("_noDD_noTwi")[0]


# In[11]:


# # temporary fix for files missing "_noDD_noTwi"

# fix_path = "remove_no_template_results_256_first_year_baseline_v3_3_10yrs_db_noDD_noTwi"
# files = glob.glob("{}/*v3_3*".format(fix_path))

# print(len(files))

# print(files[0].split("/")[-1])

# files_fix = [x for x in files if "_noDD_noTwi" not in x.split("/")[-1]]
# print(len(files_fix))

# str_to_replace = "first_year_baseline_v3_3_10yrs_"
# replace_with_str = "first_year_baseline_v3_3_10yrs_db_noDD_noTwi_"

# with open("mv_cmds.sh","w") as f:
#     f.write("#!/bin/bash\n")
#     for i in range(len(files_fix)):
#         file1 = files_fix[i]
#         file2 = "{}/{}".format(fix_path,file1.split("/")[-1].replace(str_to_replace,replace_with_str))
#         print(file1)
#         print(file2)
#         cmd = "mv {} {}\n".format(file1,file2)
#         f.write(cmd)
# #         break


# In[12]:


# query or load the database of year 1 observations
if not os.path.isfile(year1_fname):
    print("get year 1 observations")
    conn = sqlite3.connect(baseline_db)
    qry = 'select * from observations where night<{}{};'.format(night_max,sqlDD)
    print(qry)
    df_year1 = pd.read_sql(qry, conn)
    conn.close()

    # open up a connection to a new database
    conn = sqlite3.connect(year1_fname)
    df_year1.to_sql('observations', conn, index=False, if_exists='replace')
    conn.close()

else:
    print("load {}".format(year1_fname))
    conn = sqlite3.connect(year1_fname)
    df_year1 = pd.read_sql('select * from observations;', conn)
    conn.close()


# In[13]:



# # Count up area in footprint from the baseline, for reference (varies by filter)

# In[14]:



# In[15]:


### is this total "unique" area or the total area of all observed visits?
m = "CountMetric"
_runName = "{}_nside-{}".format(runName,nside).replace(".","_")

footprint_area = {}
total_pixel_area = {}
total_pixels = {}

for filt in all_filters:
    
    print(_runName,m,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, print_flag=True,
                                      save_dir = "_".join(save_dir.split("_n_visits_4_")),
                                       savedir_out  = metric_bundle_savedir,
                                        sql = sql)
    print(metric_plot.shape)
    
    x = metric_plot.data
    
    footprint_area[filt] = len(x[x>0]) * pix_area
    total_pixel_area[filt] = sum(x[x>0]) * pix_area    
    total_pixels[filt] = sum(x[x>0])


# In[16]:



# In[17]:



# In[18]:



# In[19]:


### count the total number of visits
m = "CountMetric"
_runName = "{}_nside-{}".format(runName,nside).replace(".","_")
slicer="ONED"

total_nvisits = {}
total_nvisits_area = {}

for filt in all_filters:
    
    print(_runName,m,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, 
                                       slicer=slicer,
                                       print_flag=True,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,
                                        save_dir = "_".join(save_dir.split("_n_visits_4_")))
    print(metric_plot.shape)
    
    x = metric_plot.data
    
    total_nvisits[filt] = sum(x[x>0])    
    total_nvisits_area[filt] = sum(x[x>0]) * lsst_footprint


# In[20]:



# In[21]:



# In[22]:


# determine the mean detector coverage per visit
for filt in all_filters:
    print(filt,total_pixel_area[filt]/total_nvisits_area[filt])


# # Plot metrics - counts and area for baseline

# In[23]:


m = "CountMetric"

for filt in all_filters:
    
    _runName = "{}_nside-{}".format(runName,nside).replace(".","_")
    print(_runName,m,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, print_flag=True,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                       
                                        save_dir = "_".join(save_dir.split("_n_visits_4_")))
    print(metric_plot.shape)

    title = "{} {} {}".format(_runName,m,filt)
    skymap_plot(metric_plot, title = title)
    histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)
    
    print(filt,sum(metric_plot.data))


# In[24]:


# do counts with time
m = "CountMetric"
slicer = "ONED"

for filt in all_filters:
    
    _runName = "{}_nside-{}".format(runName,nside).replace(".","_")
    print(_runName,m,slicer,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, 
                                       slicer = slicer,
                                       print_flag=True,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                       
                                      save_dir = "_".join(save_dir.split("_n_visits_4_")))
    print(metric_plot.shape)

    title = "{} {} {}".format(_runName,m,filt)

    fig = plt.figure()
    gs = gridspec.GridSpec(1,2)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    
    # plot the total nvisits per night
    ax1.plot(nights[1:],metric_plot.data)
    
    # plot the cumulative nvisits as a function of night
    ax2.scatter(nights[1:], metric_plot.data.cumsum(),label = m)
    
    # total nvisits at the end of the survey
    ax2.axhline(total_nvisits[filt],c="k",label = total_nvisits[filt])

    # plot year 1 baseline directly as a comparison
    if filt=="all":
        _df = df_year1.copy()
    else:
        _df = df_year1[df_year1["filter"]==filt]
    _df = _df[~_df["note"].str.contains("DD")]
    
    counts, bins = np.histogram(_df["night"], bins=nights)
    ax2.plot(bins[1:], counts.cumsum(), label = "df_year1", c = "r")
        
    ax1.set_xlabel("night(d)")
    ax2.set_xlabel("night(d)")
    ax1.set_ylabel("nvisits")
    ax2.set_ylabel("cumulative nvisits")
    ax2.legend()
    
    plt.suptitle("{} {} {} {}".format(_runName,m,slicer,filt))
    
    plt.tight_layout()
    
    plt.close()
#     break


# In[25]:


# do counts with time
m = "CountMetric"
slicer = "ONED"

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
    
for i,filt in enumerate(all_filters):
    
    _runName = "{}_nside-{}".format(runName,nside).replace(".","_")
    print(_runName,m,slicer,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, 
                                       slicer = slicer,
                                       print_flag=True,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                       
                                      save_dir = "_".join(save_dir.split("_n_visits_4_")))
    print(metric_plot.shape)

    title = "{} {} {}".format(_runName,m,filt)


        
    # plot the cumulative survey area as a function of night
    ax1.plot(nights[1:], metric_plot.data.cumsum() * lsst_footprint, c = "C{}".format(i))
    
    # total area at the end of the survey
    area = total_pixel_area[filt]
    area /= (total_pixel_area[filt]/total_nvisits_area[filt]) # correct for incomplete coverage of each visit
    ax1.axhline(area,label = "{}: total {:.1f} sqdeg".format(filt,area), c = "C{}".format(i))

        
ax1.set_xlabel("night(d)")
ax1.set_ylabel("cumulative area(sq deg)")
ax1.legend()
plt.suptitle("{} {} {}".format(_runName,m,slicer))

plt.tight_layout()

plt.close()


# In[26]:


# do counts with time
m = "CountMetric"
slicer = "ONED"

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
    
for filt in all_filters:
    
    _runName = "{}_nside-{}".format(runName,nside).replace(".","_")
    print(_runName,m,slicer,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, 
                                       slicer = slicer,
                                       print_flag=True,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                       
                                      save_dir = "_".join(save_dir.split("_n_visits_4_")))
    print(metric_plot.shape)

    title = "{} {} {}".format(_runName,m,filt)
    
    # total area at the end of the survey
    area = total_pixel_area[filt]
    area /= (total_pixel_area[filt]/total_nvisits_area[filt]) # correct for incomplete coverage of each visit
    
    # plot the cumulative fractional survey area as a function of night
    ax1.plot(nights[1:], (metric_plot.data.cumsum() * lsst_footprint) / area,label = "{}: total {} sqdeg".format(filt,area))
    
ax1.axhline(1,c="k")
ax1.axvline(365,c="k")
        
ax1.set_xlabel("night(d)")
ax1.set_ylabel("cumulative fraction of total area")
ax1.legend()
plt.suptitle("{} {} {}".format(_runName,m,slicer))

plt.tight_layout()

plt.close()


# In[27]:


# # do counts with time
# m = "CountMetric"
# slicer = "ONED"

# fig = plt.figure()
# gs = gridspec.GridSpec(1,1)
# ax1 = plt.subplot(gs[0,0])
    
# for filt in ["u","g","r","i","z","y"]:
    
#     _runName = "{}_nside-{}".format(runName,nside).replace(".","_")
#     print(_runName,m,slicer,filt)
#     metric_plot = load_template_metric(_runName,
#                                        metric=m,
#                                        filt=filt, 
#                                        slicer = slicer,
#                                        print_flag=True,
#                                         sql = sql,
#                                        savedir_out  = metric_bundle_savedir,
#                                       save_dir = save_dir)
#     print(metric_plot.shape)

#     title = "{} {} {}".format(_runName,m,filt)
    
#     # total area at the end of the survey
#     area = footprint_area[filt]
#     area /= (total_pixel_area[filt]/total_nvisits_area[filt]) # correct for incomplete coverage of each visit
    
#     # plot the cumulative fractional survey area as a function of night
#     ax1.plot(nights[1:], (metric_plot.data.cumsum() * lsst_footprint) / area,label = "{}: footprint {} sqdeg".format(filt,area))
# #     ax1.plot(nights[1:], (metric_plot.filled * lsst_footprint) / area,label = "{}: footprint {} sqdeg".format(filt,area))
    
# ax1.axhline(1,c="k")
# ax1.axvline(365,c="k")

# ax1.set_xlim(0,150)
# ax1.set_ylim(0,1.1)
        
# ax1.set_xlabel("night(d)")
# ax1.set_ylabel("cumulative fraction of footprint area")
# ax1.legend()
# plt.suptitle("{} {} {}".format(_runName,m,slicer))

# plt.tight_layout()

# plt.close()

# ### THIS DOES NOT COUNT UNIQUE VISITS!
# # do we need a new metric?


# # look at pair metrics for the baseline

# In[28]:


# m = "PairMetric"

# for filt in ["all","r_or_g","r_or_i","r_or_z"]:

#     _runName = "{}_nside-{}".format(runName,nside).replace(".","_")
#     print(_runName,m,filt)
#     metric_plot = load_template_metric(_runName,
#                                        metric=m,
#                                        filt=filt, print_flag=True,
#                                        sql = sql,
#                                        savedir_out  = metric_bundle_savedir,                                       
#                                       save_dir = save_dir)
#     print(metric_plot.shape)

#     title = "{} {} {}".format(_runName,m,filt)
#     skymap_plot(metric_plot, title = title)
#     histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)


# In[29]:


### TODO
# double check the area in the pair histogram


# # Metrics for template generation

# In[30]:



# In[31]:


# number of nights between first visit and template generation
# Or until first science image?

m = "doAllTemplateMetrics_reduceDeltaNight" # add a vline at the template timescale to see if it peaks

for filt in ["all","u","g","r","i","z","y"]:
    for tscale in tscales:

        # divide year 1 into chunks of a given template_timescale
        template_nights = np.arange(0,365+tscale,tscale)
        template_nights[-1] = 365 # consider only the first year

#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
        print(_runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           filt=filt,
                                           sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                        
                                           operation="min",save_dir = save_dir)
        print(metric_plot.shape)

        title = "{} {}".format(_runName,m)

    #     skymap_plot_Night(metric_plot, title = title, template_nights = template_nights,
    #                       _min = 0, _max = 365)
        skymap_plot_Night(metric_plot, title = title, template_nights = template_nights,
                          _min = 0, _max = 90)
        histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)

    #     break


# In[32]:


# m = "doPairTemplateMetrics_reducePairs"
# # filt = "all"
# # filt = "r_or_g"

# for filt in ["all","r_or_g","r_or_i","r_or_z"]:
#     for tscale in tscales:

#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
#         print(runName,m,filt)
#         metric_plot = load_template_metric(_runName,
#                                            metric=m,
#                                            sql = sql,
#                                        savedir_out  = metric_bundle_savedir,                                        
#                                            filt=filt,save_dir = save_dir)
#         print(metric_plot.shape)

#         title = "{} {}".format(_runName,m)
#         skymap_plot(metric_plot, title = title)
#         histogram_plot(metric_plot,bins=200,title=title,pix_area=pix_area)

#     #     break


# In[33]:



# In[ ]:


m = "doAllTemplateMetrics_reduceCount"
# filt = "all"

for filt in ["all","u","g","r","i","z","y"]:

    for tscale in tscales:

#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           filt=filt,save_dir = save_dir)
        print(metric_plot.shape)

        title = "{} {}".format(_runName,m)
        skymap_plot(metric_plot, title = title)
        histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)

    #     break


# In[ ]:


m = "doAllTemplateMetrics_reduceCount"
filters = ["u","g","r","i","z","y"]
tscale = 28


for filt in filters:
        
#     _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
    _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
    print(runName,m,filt)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                                                               
                                       filt=filt,save_dir = save_dir)
    print(metric_plot.shape)
    
    title = "{} {} {}".format(_runName,m, filt)
    skymap_plot(metric_plot, title = title)
    histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)
        


# In[ ]:


m = "doAllTemplateMetrics_reduceNTemplate"

for filt in ["u","g","r","i","z","y"]:

    for tscale in tscales:
        
#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           sql = sql,
                                           operation = "first",
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           filt=filt,save_dir = save_dir)
        print(metric_plot.shape)

        title = "{} {} {}".format(_runName,m, filt)
        skymap_plot(metric_plot, title = title)
        histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)
    
        


# In[ ]:


m = "doAllTemplateMetrics_reduceSeeingTemplate"

for filt in ["u","g","r","i","z","y"]:

    for tscale in tscales:
        
#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           sql = sql,
                                           operation = "first",
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           filt=filt,save_dir = save_dir)
        print(metric_plot.shape)

        title = "{} {} {}".format(_runName,m, filt)
        skymap_plot(metric_plot, title = title)
        histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)
    
        


# In[ ]:


m = "doAllTemplateMetrics_reduceDepthTemplate"

for filt in ["u","g","r","i","z","y"]:
    for tscale in tscales:
        
#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           sql = sql,
                                           operation = "first",
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           filt=filt,save_dir = save_dir)
        print(metric_plot.shape)

        title = "{} {} {}".format(_runName,m, filt)
        skymap_plot(metric_plot, title = title)
        histogram_plot(metric_plot,bins = 200,title=title,pix_area=pix_area)
    
        


# # Load all template count metrics and determine values at each template night

# In[ ]:


### do full counts for pair metrics as well? filters ["all","r_or_g","r_or_i","r_or_z"]


# In[ ]:


m = "doAllTemplateMetrics_reduceCount"

# define dicts to store results - this cell can be slow
# count the number of pixels and area covered with templates as a function of time

_runName = "{}_nside-{}_{}".format(runName,nside,m).replace(".","_")
# _runName = "{}_db{}_nside-{}".format(runName,sql,nside).replace(".","_")
templates_fname = "{}/{}{}_tscales-{}.json".format(metric_bundle_savedir,_runName,sql,"-".join(str(x) for x in tscales))
print(_runName)

# try load the dict
if os.path.isfile(templates_fname):
    print("load {}".format(templates_fname))
    template_counts = json.load( open(templates_fname) )
    
# create a new nested dict
else:
    
    template_counts = {}

    for tscale in tscales:
        
        template_counts[str(tscale)] = {}
        
        template_counts[str(tscale)]["footprint_area"] = {} # unique footprint pixel area, no repeat visits (with templates)
        template_counts[str(tscale)]["pixel_area"] = {} # total pixel area, i.e. including repeat visits (with templates)
        template_counts[str(tscale)]["pixels"] = {} # total number of pixels with templates

        template_nights = template_timescales[str(tscale)]
        print(template_nights)

        for filt in ["all","u","g","r","i","z","y"]:

            fprint_area = []
            pixel_area = []
            pixels = []
            
            x_footprint = np.zeros(hp.pixelfunc.nside2npix(nside))

            for t in template_nights[1:]:

#                 _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
                _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
                print(runName,m,filt, t)
                print(_runName)
                metric_plot = load_template_metric(_runName,
                                                   metric=m,
                                                   filt=filt,
                                                   time = t,
                                                   save_dir = save_dir,
                                       savedir_out  = metric_bundle_savedir,                                                                                           
                                                   fill_value=0,
                                                   sql = sql,
                                                  print_flag = True)

                # there may be missing metric files (e.g. no observations of u with t<7d)
                if metric_plot is None:
                    fprint_area.append(0)
                    pixel_area.append(0)
                    pixels.append(0)
                    continue
                
                print(metric_plot.shape)

                x = metric_plot.data
                
                # consider the cumulative unique footprint area
                x_footprint +=x
                
                fp_area = len(x_footprint[x_footprint>0]) * pix_area ### is this correct when looking at each timestamp separately? We should consider the number of pixels cumulatively?
                p_area = sum(x[x>0]) * pix_area
                pix = sum(x[x>0])

                fprint_area.append(fp_area)
                pixel_area.append(p_area)
                pixels.append(pix)

            template_counts[str(tscale)]["footprint_area"][filt] = fprint_area
            template_counts[str(tscale)]["pixel_area"][filt] = pixel_area
            template_counts[str(tscale)]["pixels"][filt] = pixels

    print("save {}".format(templates_fname))
    json.dump(template_counts, open(templates_fname, 'w' ) )


# In[ ]:


m = "doAllTemplateMetrics_reduceCount"

# define dicts to store results - this cell can be slow
# count the number of pixels and area covered with templates as a function of time

_runName = "{}_nside-{}_{}".format(runName,nside,m).replace(".","_")
# _runName = "{}_db{}_nside-{}".format(runName,sql,nside).replace(".","_")
templates_fname = "{}/{}{}_tscales-{}.json".format(metric_bundle_savedir,_runName,sql,"-".join(str(x) for x in tscales))
print(_runName)

for tscale in tscales:
# for tscale in [28]:
    template_nights = template_timescales[str(tscale)]
    print(template_nights)
    
#     for filt in ["all","u","g","r","i","z","y"]:
    for filt in ["y"]:

        x_footprint = np.zeros(hp.pixelfunc.nside2npix(nside))

        for t in template_nights[1:]:

#             _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
            _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")
            print(runName,m,filt, t)
            metric_plot = load_template_metric(_runName,
                                               metric=m,
                                               filt=filt,
                                               time = t,
                                               save_dir = save_dir,
                                       savedir_out  = metric_bundle_savedir,                                                                                       
                                               fill_value=0,
                                               sql = sql,
                                              print_flag = False)

            # there may be missing metric files (e.g. no observations of u with t<7d)
            if metric_plot is None:
                continue

            print(metric_plot.shape)

            x = metric_plot.filled()

            # consider the cumulative unique footprint area
            x_footprint +=x

            x_footprint_plot = x_footprint.copy()
            x_footprint_plot[x_footprint_plot==0] = np.nan
            skymap_plot(x_footprint_plot, title=t)



# In[ ]:



# In[ ]:



# In[ ]:


# sum(template_counts["7"]["pixels"]["y"])


# In[ ]:



# In[ ]:


for tscale in tscales:

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    template_nights = template_timescales[str(tscale)]
        
    for filt in 'ugrizy':
        
        visit_area = np.array(template_counts[str(tscale)]["pixel_area"][filt])
        
        cum_visit_area = visit_area.cumsum()
        
        ax1.plot(template_nights[1:],cum_visit_area)
        ax1.scatter(template_nights[1:],cum_visit_area,label = "{} {}".format(tscale,filt))

    ax1.set_xlabel("night")
    ax1.set_ylabel("cumulative area")
    ax1.legend()

    plt.close()
    


# In[ ]:


for tscale in tscales:

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    template_nights = template_timescales[str(tscale)]
        
    for filt in 'ugrizy':
        
        visit_area = np.array(template_counts[str(tscale)]["pixel_area"][filt])
        
        cum_visit_area = visit_area.cumsum()
        
        ax1.plot(template_nights[1:],cum_visit_area/total_pixel_area[filt])
        ax1.scatter(template_nights[1:],cum_visit_area/total_pixel_area[filt],label = "{} {}".format(tscale,filt))

    ax1.set_xlabel("night")
    ax1.set_ylabel("fraction of total baseline survey area")
    ax1.legend()

    plt.close()


# In[ ]:


for tscale in tscales:

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    template_nights = template_timescales[str(tscale)]
        
    for filt in 'ugrizy':
        
        visit_area = np.array(template_counts[str(tscale)]["footprint_area"][filt])
        
#         cum_visit_area = visit_area.cumsum() ### NO CUMULATIVE! already doen in template_counts
        
        ax1.plot(template_nights[1:],visit_area/footprint_area[filt])
        ax1.scatter(template_nights[1:],visit_area/footprint_area[filt],label = "{} {}".format(tscale,filt))

    ax1.set_xlabel("night")
    ax1.set_ylabel("fraction of footprint area")
    ax1.legend()
    
    ax1.set_ylim(0,1.1)
    ax1.axhline(1,c="k")

    plt.close()


# In[ ]:


# do counts with time
m = "CountMetric"
slicer = "ONED"
_runName = "{}_nside-{}".format(runName,nside).replace(".","_")
    
for filt in 'ugrizy':

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])
        
    # plot the baseline survey (no template generation)
    metric_plot = load_template_metric(_runName,
                                       metric=m,
                                       filt=filt, 
                                       slicer = slicer,
                                       print_flag=True,
                                       sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                                                               
                                      save_dir = save_dir)
    print(metric_plot.shape)
    
    # total area at the end of the survey
    area = total_pixel_area[filt]
    area /= (total_pixel_area[filt]/total_nvisits_area[filt]) # correct for incomplete coverage of each visit
    
    # plot the cumulative fractional survey area as a function of night
    ax1.plot(nights[1:], (metric_plot.data.cumsum() * lsst_footprint) / area,
             c="k",label = "baseline {}: total {} sqdeg".format(filt,area))
    
    for tscale in tscales:

        template_nights = template_timescales[str(tscale)]

        visit_area = np.array(template_counts[str(tscale)]["pixel_area"][filt])
        
        cum_visit_area = visit_area.cumsum()
        
        ax1.plot(template_nights[1:],cum_visit_area/total_pixel_area[filt])
        ax1.scatter(template_nights[1:],cum_visit_area/total_pixel_area[filt],label = "{} {}".format(tscale,filt))
        

    ax1.set_xlabel("night")
    ax1.set_ylabel("fraction of total baseline survey area")
    ax1.legend()

    plt.close()


# # Consider the deltaNight metric

# In[ ]:



# In[ ]:


# number of nights between first visit and template generation
# this metric returns the night EACH HEALPIX first got a template
# I.e. the slice is spatial, but the returned units are in days

m = "doAllTemplateMetrics_reduceDeltaNight"

for tscale in tscales:

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    for filt in 'ugrizy':

        # load the deltaNight metric
#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")

        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           filt=filt,
                                           sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           operation="min",save_dir = save_dir)        

        counts, bins = np.histogram(metric_plot.filled(), bins=nights)
        print(filt,sum(counts)*pix_area)
        area_since_first = counts.cumsum() * pix_area 

        ax1.plot(nights[1:],area_since_first, label = "{} {}".format(tscale,filt))
    
    ax1.set_xlabel("deltaNight")
    ax1.set_ylabel("cumulative area")
    ax1.legend()
    plt.close()
#     break


# In[ ]:


### why does this rise much faster than lynne's notebook? nside difference?
# https://github.com/rhiannonlynne/notebooks/blob/main/Template%20Generation.ipynb


# In[ ]:


# compare to the fraction of footprint area with templates (reduceCount)
# why does that go past fraction of 1 so quickly when deltanight does not?
# this is because deltaNight is time difference, not absolute time!


# In[ ]:


m = "doAllTemplateMetrics_reduceDeltaNight" 

for tscale in tscales:

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    for filt in 'ugrizy':

        # load the deltaNight metric
#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")

        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           filt=filt,
                                           sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           operation="min",save_dir = save_dir)        

        counts, bins = np.histogram(metric_plot.filled(), bins=nights)
        area_since_first = counts.cumsum() * pix_area 

#         ax1.scatter(nights[1:],area_since_first/total_pixel_area[filt], label = "{} {}".format(tscale,filt))
#         ax1.set_ylabel("cumulative fraction of total baseline area")

        ax1.plot(nights[1:],area_since_first/footprint_area[filt], label = "{} {}".format(tscale,filt))
        ax1.set_ylabel("cumulative fraction of footprint area")

    ax1.axhline(1,c="k")
    
    ax1.set_xlim(0,90)
    
    ax1.set_xlabel("deltaNight")
    ax1.legend()
    plt.close()
#     break


# In[ ]:


# why is there a difference between the total footprint area and the area since first visit?
# compare the y filter


# In[ ]:



# In[ ]:



# In[ ]:


# likely due to actual science frames vs first possible template time


# In[ ]:


m = "doAllTemplateMetrics_reduceDeltaNight" 

for filt in 'ugrizy':

    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    for tscale in tscales:
    
        # load the deltaNight metric
#         _runName = "{}_db{}_tscale-{}_nside-{}".format(runName,sql,tscale,nside).replace(".","_")
        _runName = "{}_tscale-{}_nside-{}".format(runName,tscale,nside).replace(".","_")

        print(runName,m,filt)
        metric_plot = load_template_metric(_runName,
                                           metric=m,
                                           filt=filt,
                                           sql = sql,
                                       savedir_out  = metric_bundle_savedir,                                                                                   
                                           operation="min",save_dir = save_dir)        

        counts, bins = np.histogram(metric_plot.filled(), bins=nights)
        area_since_first = counts.cumsum() * pix_area 

#         ax1.scatter(nights[1:],area_since_first/total_pixel_area[filt], label = "{} {}".format(tscale,filt))
#         ax1.set_ylabel("cumulative fraction of total baseline area")

        ax1.plot(nights[1:],area_since_first/footprint_area[filt], label = "{} {}".format(tscale,filt))
        ax1.set_ylabel("cumulative fraction of footprint area")
        
    ax1.axhline(1,c="k")

    ax1.set_xlim(0,90)

    ax1.set_xlabel("deltaNight")
    ax1.legend()
    plt.close()
#     break


# # analyse the visit databases directly

# In[ ]:



# In[ ]:


# db_fname = "{}/visit_cut_t-168d_nside-32.db".format(save_dir)
db_fname = "remove_no_template_results_256_noDD/visit_cut_t-28d_nside-256.db"
print(db_fname)

conn = sqlite3.connect(db_fname)
df = pd.read_sql('select * from observations;', conn)
conn.close()
    


# In[ ]:



# In[ ]:



# In[ ]:


fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

ax1.hist(x,bins = "auto")

ax1.set_xlabel("fractional template coverage")
ax1.set_ylabel("number of visits")

plt.close()


# In[ ]:



# In[ ]:


x = df["npix_template"] / n_pix

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

keep_list = []
frac_list = np.linspace(0.1,1.0,100)
# print(frac_list)    

for frac in frac_list:
    keep_list.append(len(x[x>=frac])/len(x))
    
ax1.scatter(frac_list,keep_list)

ax1.axvline(0.9,c="r")

ax1.set_xlabel("fractional template coverage")
ax1.set_ylabel("fraction of visits to keep")

plt.close()


# In[ ]:



# In[ ]:


# rate of change in fractional template coverage
# pick optimum value from the turning point?

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
    
ax1.scatter(frac_list[1:],np.diff(keep_list))
ax1.plot(frac_list[1:],np.diff(keep_list))

ax1.axvline(0.9,c="r")

ax1.set_xlabel("fractional template coverage")
ax1.set_ylabel("delta(fraction of visits to keep)")

plt.close()


# In[ ]:


# redact this visit database based on fractional template coverage
# save the redacted db to run "generate_ss" command on


# In[ ]:




