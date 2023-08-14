#!/usr/bin/env python
# coding: utf-8

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
from optparse import OptionParser

# import the template functions
from rubin_templates.rubin_templates import BaseTemplateMetric, HasTemplateIndx, remove_no_templates, doAllTemplateMetrics, doPairTemplateMetrics

parser = OptionParser()
parser.add_option( "-d", "--database", dest="baseline_db", help="baseline database file", metavar="DATABASE")
parser.add_option( "-n", "--nside", dest="nside", help="nside for healpix", metavar="NSIDE" )
parser.add_option( "-t", "--tscales", action="append", dest="tscales", help="timescales for template generation (list append)", metavar="TSCALES" )
parser.add_option( "-s", "--save_path", dest="save_path", help="save path for outpur dir", metavar="SAVE_PATH" )
parser.add_option( "-b", "--baseline", action="store_true", dest="baseline", help="flag to run baseline metrics")
parser.add_option( "-m", "--metrics", action="store_true", dest="metrics", help="flag to run template metrics")
parser.add_option( "-p", "--pairs", action="store_true", dest="pairs", help="flag to run pair template metrics")
parser.add_option( "-v", "--visits", action="store_true", dest="visits", help="flag to run visit database template redaction code")

(options,args)=parser.parse_args()

if options.baseline_db:
    baseline_db=options.baseline_db
else:
    baseline_db="baseline_v3.0_10yrs.db"
if options.nside:
    nside = int(options.nside)
else:
    nside = 256
if options.tscales:
    tscales = options.tscales
else:
    tscales = [28]
if options.save_path:
    save_path = options.save_path
else:
    save_path = ""
if options.baseline:
    baseline = True
else:
    baseline = False
if options.metrics:
    metrics = True
else:
    metrics = False
if options.pairs:
    pairs = True
else:
    pairs = False
if options.visits:
    visits = True
else:
    visits = False

print("baseline_db: ",baseline_db)
print("nside: ",nside)
print("tscales: ",tscales)
print("save_path: ",save_path)
print("baseline metrics flag: ", baseline)
print("template metrics flag: ", metrics)
print("template pairs flag: ", pairs)
print("template visits flag: ", visits)


# select the survey simulation
year1_fname = 'first_year_{}.db'.format(baseline_db.split(".db")[0])
# range of nights to investigate (year 1 of survey)
night_min = 0
night_max = 365

# set up spatial slicer with nside
s_hp = maf.HealpixSlicer(nside=nside)

# make a time slicer over the night range
nights = np.arange(night_min, night_max+1, 1)
s_t = maf.OneDSlicer(slice_col_name="night", bins = nights-0.5) # subtract half a night as slicer bins are left and right inclusive


# set path for metric output files
save_dir = "remove_no_template_results_{}_{}".format(nside,year1_fname.replace(".","_"))
if save_path != "":
    save_dir = save_path + "/" + save_dir

print(save_dir)


# use the healpix area to find approximate number of healpixels in a single visit
pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) # square degrees
lsst_footprint = 9.6 # square degrees, GET EXACT NUMBER?
n_pix = lsst_footprint/pix_area # number of healpixels within each visit footprint
print("healpixel area = {} square degrees\n number of healpixels in visit = {}".format(pix_area,n_pix))


# Ignore deep drilling fields
# This string will be applied to metric sql query and also to the redacted database code
# sqlDD = ''
sqlDD = ' and note not like "%DD%"'
if sqlDD!='':
    save_dir+="_noDD" # change the save directory
print(save_dir)


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


# set up some filenames
opsim_fname = year1_fname
print(opsim_fname)
opsdb = opsim_fname
runName = os.path.split(opsdb)[-1].replace('.db', '')
print(runName)


# template generation timescales to test
tscales = [int(t) for t in tscales]
# store the timescale and template generation nights in a dict
template_timescales = {}
for tscale in tscales:

    # divide year 1 into chunks of a given template_timescale
    template_nights = np.arange(0,night_max+tscale,tscale)
    template_nights[-1] = night_max # consider only the first year
    template_timescales[str(tscale)] = template_nights

print(template_timescales)

if baseline:
    # Do regular metrics on the baseline, count visits and pairs
    print("run baseline metrics")

    start = time.time()

    slicer1 = s_hp # spatial slicer
    slicer2 = s_t # time (night) slicer

    t_data=night_max

    _runName = "{}_nside-{}".format(runName,nside)

    bl = []
    summary_stats = [maf.MedianMetric()]

    for filt in ["all","u","g","r","i","z","y",
                            "r_or_g","r_or_i","r_or_z"]:

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

if metrics:
    # do metrics with templates for a spatial slicer
    print("run template metrics")

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


if pairs:

    # run pair metrics with templates for a spatial slicer
    # use all filters or filter pairs: r-g, r-i, r-z
    print("run template pair metrics")

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



if visits:
    # Run the full remove template code on all date ranges
    print("run visit template metrics")

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

        json.dump(template_visits, open(templates_fname, 'w' ) )


    # save the redacted database for each timescale
    # select visits in chunk from original year 1 database
    df_data = df_year1.copy()
    print("save df_data: {} rows".format(len(df_data)))

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

    end = time.time()
    dt = end-start
    print("{}s ({:2f}min or {:2f}hrs)".format(dt,dt/60,dt/60/60))

if not baseline or not metrics and not pairs and not visits:
    print("select a flag: --metrics, --pairs, --visits")
