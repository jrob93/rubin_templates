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

template_frac_list = [0.9,0.6]
lsst_footprint = 9.6 # square degrees, GET EXACT NUMBER?
night_max = 365


# In[3]:


# db_files = glob.glob("remove_no_template_results*/*visit_cut*.db")
db_files = glob.glob("visit_cut_dbs/remove_no_template_results*/*visit_cut*.db")
db_files = [x for x in db_files if "frac" not in x]
print(db_files)

# In[4]:


for template_frac in template_frac_list:
    for i in range(len(db_files)):
        # get the filename and nside info
        dbf = db_files[i]
        _dbf = dbf.split("_nside-")[-1].split(".db")
        nside = int(_dbf[0])
        print("file: {}".format(dbf))
        print("nside: {}".format(nside))
        print("frac: {}".format(template_frac))

        # use the healpix area to find approximate number of healpixels in a single visit
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) # square degrees
        n_pix = lsst_footprint/pix_area # number of healpixels within each visit footprint
        print("healpixel area = {} square degrees\n number of healpixels in visit = {}".format(pix_area,n_pix))

        # open the db file
        conn = sqlite3.connect(dbf)
        df = pd.read_sql('select * from observations;', conn)
        df = df[df["night"]<=night_max]
        conn.close()

        print("number of visits with fraction_template>0: {}".format(len(df)))

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




