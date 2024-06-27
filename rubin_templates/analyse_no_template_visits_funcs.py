import numpy as np
import glob
import os
from natsort import natsorted
import rubin_sim.maf as maf
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def load_template_metric(runName,
                         metric="doAllTemplateMetrics_Count",
                         slicer = "HEAL",
                         filt="all",save_dir="remove_no_template_results_32",
                         time = False,
                         operation="sum",
                         fill_value=np.nan,
                        print_flag=False,
                        overwrite=False,
                        savedir_out="analysis_metric_bundles",
                        sql=False):
    
    """
    The doAllTemplateMetrics works by running metrics on batches of observations made between nights on which templates were generated.
    This function loads the all the files of a previously run metric.
    An operation can be performed on all of these reloaded metrics,
    e.g. the total number of visits from the "Count" metric can be found by using "sum".
    
    time allows you to select a single metric file for some time batch
    
    The save_dir is searched for files matching some pattern: runName, metric, filt, slicer
    """
        
    if time:
        metric_fname="{}/{}_{}_{}_lt_{}_{}.npz".format(savedir_out,runName,metric,filt,time,slicer)
    else:
        metric_fname="{}/{}_{}_{}_{}.npz".format(savedir_out,runName,metric,filt,slicer)
        
    ### WHERE DOES _runName get defined?
    
    if sql:
        metric_fname = metric_fname.split(".npz")[0]+sql+".npz"
        
    print(metric)
    print(metric_fname)
    print(runName)
        
    if os.path.isfile(metric_fname) and not overwrite:
        print("load {}".format(metric_fname))
        metric_vals = np.load(metric_fname,allow_pickle=True)
        return metric_vals
    
    else:
        # glob all files in the save_dir
        files = glob.glob("{}/*npz".format(save_dir))

        # select the relevant metric files using the patterns
        if filt=="all": # select all filters
            _files = [x for x in files if (runName in x) 
                     & (metric in x)
                     & (slicer in x)
                     & ~(("_u_".format(filt) in x)
                     | ("_g_".format(filt) in x)
                     | ("_r_".format(filt) in x)
                     | ("_i_".format(filt) in x)
                     | ("_z_".format(filt) in x)
                     | ("_y_".format(filt) in x))
                     ]
        elif "_or_" in filt:
            _files = [x for x in files if (runName in x) 
                      & (metric in x) 
                      & (slicer in x)
                      & ("_{}_".format(filt) in x)
                      & ("_or_" in x)]
        else:
            _files = [x for x in files if (runName in x) 
                      & (metric in x) 
                      & (slicer in x)
                      & ("_{}_".format(filt) in x)
                      & ("_or_" not in x)]

        if time: # trim down based on a specific timestamp
            _files = [x for x in _files if ("_lt_{}_".format(str(time)) in x)]        

        _files = natsorted(_files) # ensure files are in numerical order (based on the template generation night considered)
        if print_flag:
            print(_files)
        print(len(_files))

        if len(_files)==0:
            # there are no files
            print("no files")
            return

        # load all these metric files
        metric_bundle = []
        for x in _files:
            metric_bundle.append(maf.MetricBundle.load(x))

        # retrieve data/mask for all the masked arrays
        # each value in data corresponds to the metric value for some slice point
        data = [mb.metric_values.data for mb in metric_bundle]
        mask = [mb.metric_values.mask for mb in metric_bundle]
#         metric_data = np.ma.array(data, mask=mask)
        metric_data = np.ma.array(data, mask=mask, fill_value = fill_value)
        print(metric_data.shape)

        if operation == "min":
            # find the minimum of all constituent metrics
            metric_vals = metric_data.min(axis=0)
        elif operation == "sum":
            # find the sum of all constituent metrics
            ### Sum of just one file should produce just the file?
            metric_vals = metric_data.sum(axis=0)
        ### Add a new operation for DeltaNight and NTemplate
        elif operation == "first":
            # take the first not badval (nan?) value
            # procedure to get first non nan value in each column of the multidimensional array
#             ind_list = []
#             for i in range(metric_data.shape[1]):
#                 _x = metric_data.filled()[:,i]
#                 _x = _x[~np.isnan(_x)]
#                 if len(_x)==0:
#                     ind_list.append(np.nan)
#                 else:
#                     ind_list.append(_x[0])
#             ind_list = np.array(ind_list)
#             metric_vals = np.ma.array(ind_list, mask=np.isnan(ind_list))
            x=metric_data.filled()
            ind_list = [x[:,i][~np.isnan(x[:,i])][0] if len(x[:,i][~np.isnan(x[:,i])])>0 else np.nan for i in range(x.shape[1])]
            metric_vals = np.ma.array(ind_list, mask=np.isnan(ind_list))
        else:
            # return the untouched multidimensional masked array if you want to look at it separately
            print("no operation performed")
            metric_vals = metric_data

        print(metric_vals,len(metric_vals))

        # replace the masked values with nan
        metric_vals.fill_value = fill_value
        
        ### N.B. choose the fill value wisely
        # e.g. must be nan for combining deltaNight

        if not time:
            # only save if time is not specified (otherwise we're just making copies of the files)
            print("save {}".format(metric_fname))
            metric_vals.dump(metric_fname)
            
        return metric_vals
    
def skymap_plot_Night(metric_plot, title, 
                      template_nights,
                      _min=None, _max=None):

    #plot the skymap
    x = hp.mollview(metric_plot, title=title, 
                    min = _min, max=_max,
                    cbar = None)
    hp.graticule()

    # customise the colorbar
    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]
    # cbar = fig.colorbar(image, ax=ax, orientation = "horizontal",aspect = 30, location = "bottom")
    cbar = fig.colorbar(image, ax=ax, orientation = "horizontal", shrink = 0.5, location = "bottom",
                       pad = 0.05)

#     # fix the ticks at 0 nights and at the end?
#     ticks = template_nights[:-1]
#     cbar.set_ticks(ticks)
#     cbar.set_ticklabels(ticks)
    
    # set vmin, vmax to get constant scale for comparision?

#     fname = "{}.png".format("".join(title.split(" ")))
#     plt.savefig(fname, facecolor="w", transparent=True, bbox_inches="tight")
    
    plt.show()

def skymap_plot(metric_plot, title):

    #plot the skymap
    x = hp.mollview(metric_plot, title=title, 
                    cbar = None)
    hp.graticule()

    # customise the colorbar
    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]
    # cbar = fig.colorbar(image, ax=ax, orientation = "horizontal",aspect = 30, location = "bottom")
    cbar = fig.colorbar(image, ax=ax, orientation = "horizontal", shrink = 0.5, location = "bottom",
                       pad = 0.05)

#     fname = "{}.png".format("".join(title.split(" ")))
#     plt.savefig(fname, facecolor="w", transparent=True, bbox_inches="tight")

    plt.show()

def histogram_plot(metric_plot,bins="auto",title="hist_plot",pix_area=None):
    
    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax1 = plt.subplot(gs[0,0])

    n,b,p = ax1.hist(metric_plot, bins = bins, histtype = "step")

    # ax1.axvline(np.median(data),c="C{}".format(i))

    
    ### CONVERT UNITS OR SCALE AXIS?
    if pix_area:
        
        # total area
        print(sum(n)*pix_area)
        
        # scale the y axis to get sky area
        y_vals = ax1.get_yticks()
        ax1.set_yticklabels(['{:3.0f}'.format(x * pix_area) for x in y_vals])
        ax1.set_ylabel("area (square degrees)")
    else:
        ax1.set_ylabel("number of healpixels")

    ax1.set_xlabel("metric number")

    plt.title(title)

    plt.show()