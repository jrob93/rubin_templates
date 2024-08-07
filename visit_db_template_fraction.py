#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import glob
import healpy as hp
import sqlite3


# In[2]:

# script to read a visit db with template coverage column
# and redact it based on some fractional coverage cut

# template_frac_list = [0.9,0.6]
template_frac_list = [0.9]
lsst_footprint = 9.6 # square degrees, GET EXACT NUMBER?
night_max = 365


# In[3]:


# db_files = glob.glob("remove_no_template_results*/*visit_cut*.db")
# db_files = glob.glob("visit_cut_dbs/remove_no_template_results*/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results*baseline_v3_2*/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results*baseline_v3_3*/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results_256_first_year_*_noDD/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results_256_first_year_*_noDD_noTwi/*visit_cut_t-3d*.db")
# db_files = glob.glob("remove_no_template_results_256_first_year_baseline_v3_3_10yrs_db_noDD_noTwi/*visit_cut_t-*.db")
# db_files = glob.glob("remove_no_template_results*baseline_v3_3*override*/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results*baseline_v3_3*n_visits_4*/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results*baseline_v3_3*override-g_n_visits_4*/*visit_cut*.db")
# db_files = glob.glob("remove_no_template_results*baseline_v3_4*n_visits_4*/*visit_cut*.db")
db_files = glob.glob("remove_no_template_results*baseline_v3_4*/*visit_cut*.db")

# db_files = [x for x in db_files if "frac" not in x]
# db_files = [x for x in db_files if "frac" not in x and "noDD_noTwi" in x]
db_files = [x for x in db_files if "frac" not in x and "noDD" in x and "v3_4" in x and "n_visits_4" not in x]
print(db_files)
# exit()

# In[4]:

sqlDD = ' and note not like "%DD%" and note not like "%twilight%"'
qry = 'select * from observations where night<{}{};'.format(night_max,sqlDD)
print(qry)

for template_frac in template_frac_list:
    for i in range(len(db_files)):
        # get the filename and nside info
        dbf = db_files[i]
        _dbf = dbf.split("_nside-")[-1].split(".db")
        nside = int(_dbf[0])
        # print("file: {}".format(dbf))
        # print("nside: {}".format(nside))
        # print("frac: {}".format(template_frac))

        # use the healpix area to find approximate number of healpixels in a single visit
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) # square degrees
        n_pix = lsst_footprint/pix_area # number of healpixels within each visit footprint
        # print("healpixel area = {} square degrees\n number of healpixels in visit = {}".format(pix_area,n_pix))

        # open the db file
        conn = sqlite3.connect(dbf)
        df = pd.read_sql(qry, conn)
        # df = df[df["night"]<=night_max]
        conn.close()

        print(dbf.split("/")[-1])
        print("number of visits with fraction_template>=0: {}".format(len(df)))
        mask_df = (~((df["note"].str.contains("twilight")) |
        (df["note"].str.contains("DD"))) &
        (df["night"]<365))
        print(len(df[mask_df]))

        # continue

        # calculate the fractional template coverage from npix
        df["fraction_template"] = df["npix_template"] / n_pix

        # remove visits below fractional template threshold
        frac_mask = df["fraction_template"]>=template_frac
        print("with template: {}, without template: {}".format(len(df[frac_mask]),len(df[~frac_mask])))

        # save the db with a new name
        df_out = df[frac_mask].reset_index(drop=True) # reset the dataframe index 
        dbf_out = dbf.split(".db")[0]+"_frac{}.db".format(int(template_frac*100))

        # open up a connection to a new database
        conn = sqlite3.connect(dbf_out)
        # save reduced visit dataframe to sql
        df_out.to_sql('observations', conn, index=False, if_exists='replace')
        conn.close()
        print("save {}\n".format(dbf_out))


# In[ ]:




