#!/usr/bin/env python
# coding: utf-8

# future note: for the problem of the height bins, a map-roll would be sufficient to roll each profile by the correct number of bins to all be masl-height aligned, and then drop/clip the first/last n,m bins
# # V8.2.1: CHANGES FROM 9-3_build_coverandheights_grid_v8.2.py:
#   * fixed a minor bug (affecting 4/~120 months) where profiles with longitude==180.0 get assigned to the wrong gridbox
# # V8.2: CHANGES FROM 9-2_build_coverandheight_grid_v8.1.py
#   * rigorous selection of profiles in time range. Implemented via:
#      * adding the last granule from the day before the start day into the intermediates
#      * set profiles to NaN which are not within the first day - last day range
#   * changed height levels (rmcp_heights var to height_midpoints) to reflect full data record (survey in 9-2-3 and 9-2-5)
#      * added an extra height bin to end the same or lower than 3S-RMCP (to get to -400 m) (so it's 78 levels now instead of 77)
#      * stored height ranges (min and max) for each bin as coordinates (since it's not a constant 240 m thing)
#   * changed filenames (fname) to look more official
#   * changed 680 mb and 440 mb bin levels to reflect ISA plus new height levels
# # V8.1: CHANGES FROM 9-1_build_coverandheight_grid_v8.0.py
#   * added radar_surface_clutter_counts_on_levels field to output
#   * added Data_quality filter to set all profiles with nonzero Data_quality to NaN
#   * bug fix where apply_all_cloud_defs unique variants return >100% cloud cover
# # V8.0: CHANGES FROM 8-7-6_build_coverandheight_grid_v7.2.1.py
#   * added apply_all_cloud_defs and def='all' mode.
#   * moved v7.2.1 defs to cloudsat_util_08
#   * GOALS FOR V8.1:
#      * careful screening of granules to be in the proper month
#      * filter out data quality flags
#      * change read_convert_merge_granule to rotate profiles where height falls outside the bin range
#      * double-check height coordinates and tweak cloud levels to match ISA 440 and 680 mb
# # V7.2.1: CHANGES FROM 8-7_build_coverandheight_grid_v7.2.py:
#   * just some minor bug fixes. The first two items are one-line modifications already made in CURC v7.2 (but not cires-local v7.2)
#	   * modified to work on alpine directories (i.e. changed scratch/intermediates directory)
#	   * increased number of workers (since on alpine for monthspredoop_olidar_g5 runs I was getting memory overflow + increased mem+cores)
#	   * modified read_GEOPROF to handle an uncommon "start_time" swath attr issue where seconds = 60 (it gets the start time from the filename instead)
#	   * modified predecode_vfm to handle exception when CloudFraction is missing a "_FillValue" attribute for dummy lidar granules (-k radar)
#          * modified main() to manually add CloudFraction and CPR_Cloud_mask attributes if keep_granule (-k) == 'radar' or 'lidar' respectively

print('importing...',end='')
import os
import time #UNUSED
import datetime #UNUSED
import matplotlib as mpl #UNUSED
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd #used for to_datetime features
import xarray as xr
import cartopy.crs as ccrs         # import projections (unused?)
import cartopy.feature as cf       # import features (unsued?)
from datetime import datetime,timedelta      # for formatting names during read/write/plot

from dask.distributed import Client, progress, LocalCluster # used for scheduler/compute cluster
from dask import delayed, compute  # currently unused (1/27, 7-3.ipynb)
import dask.bag as db              # used only for writing nc4s to disk (ibid.)
import dask.multiprocessing        # where does this get used? maybe I saw it in a S/O answer?

from cloudsat_util_08 import *

import sys,getopt                  # for processing command-line arguments

#ONLY NECESSARY WHILE GET_DOOP_FLAG IS IN THIS SCRIPT
import json                        # for reading in doop-spline data
from scipy import interpolate      # for using doop-splines
print('done')

#NEW METHOD v6
def print_cluster_info(client):
    info=client.scheduler_info()
    n_workers=len(info['workers'])
    n_cores=sum(client.ncores().values())
    n_threads=sum(client.nthreads().values())
    memory=sum([w['memory_limit'] for w in info['workers'].values()])/1073741824 #bytes per Gibibyte
    print(f'cluster has {n_cores} cores with {n_workers} workers\
, {n_threads//n_workers} threads per worker, and {memory/n_workers:.2f} GiB memory/worker')
    
#DEV NOTE: CHANGED V8.1
#v8.1: removed setting surface clutter to nan
#v7.2: added del statements
def predecode_cpr(in_da,SurfaceHeightBin,xdim='profile',ydim='bin',clip=True):
    # set surface clutter and subsurface bins to missing_value
    
    #print(in_da.attrs)
    #turn to float (to enable nan)
    da = in_da.astype(float,keep_attrs=True)
    
    scale_factor = da.attrs['factor']
    offset = da.attrs['offset']
    nan = da.attrs['_FillValue']
    valid_min, valid_max = da.attrs['valid_range']
    
    #v8.1 removed the setting of surface clutter to NaN
    
    # set all subsurface bins to nan
    bin_arr    = np.broadcast_to(np.arange(0,125),(1,37100,125))
    bin_da     = xr.DataArray(data=bin_arr,dims=('start_time','profile','bin'))
    surf_da, _ = xr.broadcast(SurfaceHeightBin,bin_da)
    #surf bin is 1-indexed not 0-indexed
    da.data = xr.where(bin_da>=surf_da,nan,da.data) 
    del nan, scale_factor, offset, valid_min, valid_max, bin_arr, bin_da, surf_da
    return da
    
#NOTE DEV: CHANGE V8.1
# added radar surface clutter
# added predecode_vfm
#added predecode vfm and attenuated mask
# changed to assign instead of rebuild dataset
def proc_and_mask(in_cv,radar_thresh=20,lidar_thresh=1,attenuate_lidar=True, count_clutter=True):

    #decode and mask to binary for vfm
    vfm_decoded = decode(predecode_vfm(in_cv['CloudFraction'],in_cv['SurfaceHeightBin']))
    vfm_binary  = mask(vfm_decoded,lidar_thresh=lidar_thresh)

    #decode and mask to binary for cpr
    # not the same because we save surface clutter before we get rid of it
    cpr_decoded = decode(predecode_cpr(in_cv['CPR_Cloud_mask'],in_cv['SurfaceHeightBin'])) #mask subsurface bins, decode to float adding NaNs
    cluttered_bins = cpr_decoded == 5 #save the surface-and-above bins with surface clutter
    cpr_decoded = cpr_decoded.where(cpr_decoded!=5) #set surface clutter to NaN
    cpr_binary  = mask(cpr_decoded,radar_thresh=radar_thresh) #then apply the threshold

    #set regions where lidar is attenuated to nan
    if attenuate_lidar == True:
        combined_vis = calculate_combined_vis(cpr_binary,vfm_binary)
        vfm_binary, attenuated_bins = mask_attenuated_lidar_nods(vfm_binary,combined_vis)

    #combined to dataset
    vfm_binary.name = 'vfm_binary'
    cpr_binary.name = 'cpr_binary'

    vfm_binary.attrs['thresh'] = lidar_thresh
    cpr_binary.attrs['thresh'] = radar_thresh

    new_vars = {'vfm_binary':vfm_binary,
                'cpr_binary':cpr_binary}

    if attenuate_lidar: 
        new_vars['attenuated_lidar_counts'] = attenuated_bins
    if count_clutter:
        new_vars['radar_surface_clutter_counts'] = cluttered_bins

    in_cv = in_cv.drop_vars(('CPR_Cloud_mask','CloudFraction')).assign(new_vars)
    
    del vfm_binary, cpr_binary, cpr_decoded, vfm_decoded, new_vars, cluttered_bins
    return in_cv

#NOTE DEV: CHANGED FOR V8.1
#ADDED V8.1
def print_dataquality(tcv):
    '''takes a Dataset with field 'Data_quality' from 2B-GEOPROF and prints
    the frequency of nonzero (i.e. bad) data quality profiles and what
    fraction of the bad profiles each flag makes up (e.g. which flags dominate the bad data)'''
    
    #load data quality values
    dq = tcv.Data_quality.to_numpy().flatten()

    #get counts of ints for which data_quality is non-nan
    #7-bit int, so possible values are 0 to 255
    counts, bins = np.histogram(dq,bins=np.arange(0,256))

    #get list of arrays showing which ints have which bits true
    #array i in list shows whether array index int j has bit i true or false
    bit_groups = [(np.arange(0,255) & 2**n).astype(bool) for n in range(8)]

    #from 1B-CPR PDICD explaining Data_quality bit flags
    bit_meanings = {0: 'RayStatus_validity not normal',
                    1: 'GPS data not valid',
                    2: 'Temperatures not valid',
                    3: 'Radar telemetry data quality not normal',
                    4: 'Peak power not normal',
                    5: 'CPR calibration maneuver',
                    6: 'Missing frame',
                    7: 'Data advisory, check website'}

    print('PRINTING DATA QUALITY INFORMATION')
    counts_nonzero = np.sum(counts[1:])
    print(f'out of {len(dq)} profiles, {counts_nonzero} are bad, which is {counts_nonzero/len(dq)*100:.02f}%')
    print('the breakdown of bad profiles by bit is')
    for bit,mask in enumerate(bit_groups):
        counts_bit = counts[mask].sum()
        print(f'bit {bit}: {counts_bit} is {counts_bit/counts_nonzero*100:.02f}% ({bit_meanings[bit]})')


#NOTE DEV: CHANGED FOR V8.1

#NOTE V7.2: I might have changed this to handle the new qc vars, but 
# it also might be exactly the same at cloudat_util_07

# v8.2: changed height thresholds to be on the ISA + better height ranges from the survey
# removed min2 for level types
# v8.1: revise unique counts methodology:
#    -- it's a 'unique high' Y/N count if there's high cloud in the relevant region and there ISN'T cloud at the other levels.
#       this makes no stipulation about the surface pressure in the atmospheric column. There is a 'unique mid' cloud if all the cloud in the column
#       is at the stated level. so the condition then is simply that the level of interest be defined.
#    -- so the revision is: total count = current level counts (then this also shouldn't need to do `cloud_counts = cloud_counts & total_counts`)
def apply_all_cloud_defs(ds,bounds={'lo':91,'hi':77}):
    #680 mb: 3238.1 masl
    #440 mb: 6503 masl
    #3238.1 m: bin 91 (2997 to 3237 m) (yes I know it's one lower it's right on the boundary)
    #6503 m: bin 77 (6355 to 6595 m)
    
    # takes counts dataset and applies a cloud definition to classify
    # whether a profile is clear or cloudy.
    
    #takes counts dataset with dims {'start_time','profile','bin'},
    # applied a definition of a cloudy profile given by 'alg',
    # and returns a dataset with dims {'start_time','profile'}
    
    #ds: binary cloud_counts,total_counts with dims {'start_time','profile','bin'}
    #alg: definition of a cloudy profile
    #bounds: bin thresholds for high/middle/low clouds
    
    #helper function for at least 2 contiguous bins
    min2 = lambda cm: cm.shift({'bin':1},fill_value=False) & cm

    #store each counts of type in a dict
    cloud_counts_by_type = { }
    total_counts_by_type = { }
    
    
    #apply different definitions of a cloudy profile
    #minimum N contiguous cloudy bins in a profile
    #just process cloud counts
    for alg in ('any','thick'):
        #definition of thick clouds
        if alg=='thick': alg_loc='min20' #20x240m=4.8km
        if alg=='any': alg_loc='min1' #240 m
        cloud_counts = ds.cloud_counts
        N=int(alg_loc[3:])
        if N==1:
            cloud_counts_by_type[alg] = cloud_counts
        elif N==2:
            cloud_counts_by_type[alg] = min2(cloud_counts)
        elif N>2:
            #bfill method shown to be quickest in 8-3-3 def minN_bfill2
            #requires the package 'bottleneck' to be installed
            fnan = xr.where(cloud_counts==True,np.nan,0) #replace True with NaN
            fill = fnan.bfill(dim='bin',limit=N-1) #fill up to N-1 nans with 0
            cloud_counts_by_type[alg] = xr.where(np.isnan(fill),True,False) #replace NaN with True
        
        total_counts_by_type[alg] = ds.total_counts #for any, thick, say if >=1 bin in profile is defined then it's a valid profile
        
    
    #NB: low cloud cover fraction is a matter of some subtlety, since
    # there are fewer low cloud counts than middle or high. So the issue
    # is whether or not total_counts should be the same or different for 
    # different vertical categories. And when calculating the unique
    # cloud cover fractions, which version of the total counts to use
    # at all. I think it would be best for total_counts to vary by type
    # and just include some duplicates e.g. for any/thick. For the unique
    # clouds, I think the correct definitions are for e.g. unique low
    #   * cloud_counts: low and not (high or mid)
    #   * total_counts: low and high and mid
    # but then this means profiles with thick clouds would be excluded
    # from the fraction if the low bins were not observed.
    
    #non-unique vertical definitions
    #process cloud_counts and total_counts
    for alg,bin_slice in {'high':slice(0,bounds['hi']),
                      'middle':slice(bounds['hi'],bounds['lo']),
                      'low':slice(bounds['lo'],None)}.items():
        #process cloud_counts
        da = ds.cloud_counts
        cloud_counts_by_type[alg] = da.isel({'bin':bin_slice})

        #process total_counts
        da = ds.total_counts
        total_counts_by_type[alg] = da.isel({'bin':bin_slice})
        
    
    #check if any bin in profile satisfies
    #this is needed for the unique calculations
    cloud_counts_by_type = {alg:cc.any(dim='bin') for alg,cc in cloud_counts_by_type.items()}
    total_counts_by_type = {alg:cc.any(dim='bin') for alg,cc in total_counts_by_type.items()}
    
    #unique vertical definitions
    # (note we must return within the case since
    #  we must fully calculate all three options)
    for alg in ('uniquehigh','uniquemiddle','uniquelow'):        
        #get names of current and other two levels
        levels = ['high','middle','low']
        current_level = alg[6:] #e.g. high, middle, low
        other1_level, other2_level = list(set(levels)-set([current_level])) #remove current level
        
        #retrieve cloud counts variables (T/F by profile)
        current = cloud_counts_by_type[current_level]
        other1  = cloud_counts_by_type[other1_level]
        other2  = cloud_counts_by_type[other2_level]
        
        #compare current key against other two keys
        cloud_counts = current & ~(other1 | other2) #unique X if no Y or Z
        
        #retrieve total counts variables
        current = total_counts_by_type[current_level]
        other1  = total_counts_by_type[other1_level]
        other2  = total_counts_by_type[other2_level]
        
        #compare current key against other two keys
        #changed for v8.1
        total_counts = current #obs X possible if the level is defined (same as non-unique total counts)
        
        #put in dict
        cloud_counts_by_type[alg] = cloud_counts
        total_counts_by_type[alg] = total_counts
        
        
    #assign everything as a new data variable back to dict
    #vars of format cloud_counts_high, cloud_counts_any, etc...
    cloud_counts_by_type = {f'cloud_counts_{alg}':ds_i for alg,ds_i in cloud_counts_by_type.items()}
    ds = ds.assign(cloud_counts_by_type)
    total_counts_by_type = {f'total_counts_{alg}':ds_i for alg,ds_i in total_counts_by_type.items()}
    ds = ds.assign(total_counts_by_type)
    
    #attenuated_lidar_counts anywhere in column (still useful)
    #v8.1: put this in a try/except so it doesn't break when attenuate_lidar==False
    try:
        ds['attenuated_lidar_counts'] = ds['attenuated_lidar_counts'].any(dim='bin')
    except KeyError:
        print('apply_cloud_defs sees no attenuated_lidar_counts field')
    
    #remove all unused (profile,bin) variables
    ds = ds.drop_dims('bin')
    
    #remove

    return ds

#DEV NOTE: CHANGED THIS FOR V8.2
#changed v8.2 to get the last granule from the day before the start date
def get_granule_paths(dayfi,dayla,year_path):
    #for cloudsat ftp server filestructure
    #get all granules in day/year range when granules are in folders organized by day
    
    #get directories in time range
    all_daydir_names = list(listdir_nohidden(year_path)) #list all day folders under year
    all_daydir_paths = [os.path.join(year_path,fname) for fname in all_daydir_names] #get full paths
    all_days = np.array([int(i) for i in all_daydir_names]) #convert to int array
    daymap  = (dayfi<=all_days) & (all_days<=dayla) #boolean map of days to keep
    keep_daydir_paths = np.asarray(all_daydir_paths)[daymap] #list of folders to use
    keep_daydir_paths = [str(s) for s in keep_daydir_paths] #use list of strs instead of numpy

    #get filepaths in list of directories
    granule_paths = [[os.path.join(path,name) for name in os.listdir(path)] for path in keep_daydir_paths] #nested list
    granule_paths = flatten_list(granule_paths) #flattened list
    granule_paths = sorted(granule_paths,key=lambda s: os.path.split(s)[1].split('_')[1]) #sort by granule number
    #NB 5/21: I think 'intersection' might not actually require filepath lists sorted by granule number,
    # but I added it here anyway since listdir_fullpaths sorts by granule
    
    #NEW V8.2
    #since granules don't begin exactly at the start of the day, include the last granule from the previous day
    path,year = os.path.split(year_path) #get year back from year_path
    dt_dayfi = pd.to_datetime(f'{year}-{dayfi}',format='%Y-%j') #convert year+dayfi to datetime
    dt_prev = dt_dayfi-pd.Timedelta(1,unit='days') #get the previous day
    daypath_prev = os.path.join(path,str(dt_prev.year),f'{dt_prev.dayofyear:03d}') #get filepath to previous day
    if os.path.exists(daypath_prev): #if there's data for the previous day,
        granulepath_prev = sorted(listdir_fullpaths(daypath_prev))[-1] #get the last granule of that previous day
        granule_paths.insert(0,granulepath_prev) #add it to the start of the list of granules to use
    
    return granule_paths


def main():
    '''Script version of 8-3-3_dev6.1-NbinCode.ipynb. Build gridded column-integrated
    cloud cover product from CloudSat 2B-GEOPROF and CALIPSO 2B-GEOPROF-LIDAR granules.

    Options:
    -h          Print help
    -p          Parent dir (contains radar & lidar folders)
    -r          Radar folder name
    -l          Lidar folder name
    -c          Check intermediates (compare originals to intermediates) #CHANGED
    -k          Granules to keep from input folders ('radar'/'lidar'/'both')
    -o          Use only one kind of the included granules in calculation ('radar'/'lidar'/'both')
    -m          Lidar threshold to be considered cloud (0-100)
    -n          Radar threshold to be considered cloud (0-40)
    -g          Grid spacing in degrees 
    -f          Filename suffix
    -t          Number of diagnostic printouts
    -v          Type of output ('height'/'cover'/'both')
    --na        Don't mask bins where lidar is likely attenuated
    --nb        Don't build intermediates nc4s 
    --nd        Don't delete intermediates on script finish 
    --def=      Definition choice for labeling a profile as cloudy.
                Default is 'min1'. Full options:
                        minN: where N is an int, at least N contiguous cloudy bins
                        high: >7 km asl cloud present+min1 
                        middle: >2.6 <7 km asl cloud present+min1
                        low: <2.6 km asl cloud present+min1 
                        uniquehigh: high but no low or middle 
                        uniquemiddle: mid but no high or low 
                        uniquelow: low but no high or mid 
                        thick: min20 (>4.8 km thick)
                        any: min1 (240 m thick)
                        all: compute all above definitions (except 'minN') and save to a single netCDF
    -s          Summit mode (place first), changes -p,-r,-l defaults and enables:
    --outdir=            base path in which to save output .nc and .pngs
    --scratchdir=        where to save intermediates
    --localdir=          local storage for dask worker files
    --threads=           number of threads for dask cluster
    --month=             MM month to process (incompatible with dayfi/dayla)
    --year=              YY or YYYY year to process (required with -s)
    --dayfi=             DDD first day inclusive (incompatible with month)
    --dayla=             DDD last day inclusive (incompatible with month)
                Note when using -s, '--year' must be specified AND EITHER
                '--month' OR '--dayla' and '--dayfi', but not both

    Examples:
    ./8-4_build_cloudcover_grid-v6.1.py -s --month=06 --year=13 -g 2.5
        Statistics for June 2013 on a 2.5 degree grid when using Summit HPC

    '''

    ### OPTIONS ###
    parent_dir = '/Users/will/Documents/cloudsat_data'
    radar_folder = 'feb2010_R05' #folder with radar HDF4 granules
    lidar_folder = 'feb2010_R05_lidar' #folder with lidar HDF4 granules
    build_intermediates = True #save nc4s to disk
    delete_intermediates = True #delete nc4s
    check_intermediate = False #compare nc4s to originals
    attenuate_lidar = True #mask bins where lidar likely attenuated
    include_only = 'both' #don't include 2B-GEOPROF-lidar in final output
    lidar_thresh = 50 #lidar cloudfraction >= thresh CHANGED FOR v7.2
    radar_thresh = 20 #cpr_cloud_mask >= thresh CHANGED FOR V7.2
    grid_spacing = 2.5 #gridcell size
    fname_suffix = '' #string to add to end of fname
    keep_granule = 'both' #rule for selecting granules from input folders
    num_printouts = 0 #number of diagnostic quicklooks to print out
    cloudy_def = 'all' #condition for considering a profile cloudy (CHANGED v6.2)
    recognized_defs = ('min1','min2','min20','high','middle','low','uniquehigh',
                       'uniquemiddle','uniquelow','thick','all')
    standard_defs = ('any', 'thick', 'high',
             'middle','low', 'uniquehigh',
             'uniquemiddle', 'uniquelow') #defs to compute when def='all'

    vertically_resolved = 'both' #compute 2D cover, 3D occurrence, or both ('cover','height','both')

    #added v7.1: call total/cloud/fraction vars different names for lat/lon vs lat/lon/height outputs
    #v7.2: added attenuated_lidar_counts
    #added v7.1: call total/cloud/fraction vars different names for lat/lon vs lat/lon/height outputs
    #v7.2: added attenuated_lidar_counts
    #v8.1: added radar_surface_clutter_counts
    rename = {'height':{'total_counts':'total_counts_on_levels',
                        'cloud_counts':'cloud_counts_on_levels',
                        'cloud_cover':'cloud_fraction_on_levels',
                        'attenuated_lidar_counts':'attenuated_lidar_counts_on_levels',
                        'radar_surface_clutter_counts':'radar_surface_clutter_counts_on_levels'},
              'cover':{'total_counts':'total_counts_in_column',
                       'cloud_counts':'cloud_counts_in_column',
                       'cloud_cover':'cloud_cover_in_column',
                       'attenuated_lidar_counts':'attenuated_lidar_counts_in_column',
                       'radar_surface_clutter_counts':'radar_surface_clutter_counts_in_column'}} #this is unnecessary because apply_all_cloud_defs doesn't output this for calculation
    
    #added v7.1
    dtypes = { }
    for kind in ('height','cover'):
        dtypes[rename[kind]['cloud_counts']] = int
        dtypes[rename[kind]['total_counts']] = int
        dtypes[rename[kind]['cloud_cover']] = float

    #new v8.2: results of survey of all granules from 8 July 2006 onward (see notebook 9-2-5) 
    height_midpoints = np.array([17987. , 17747. , 17507. , 17267. , 17027. , 16788. , 16548. ,
       16308. , 16068. , 15828. , 15588. , 15349. , 15109. , 14869. ,
       14629. , 14389. , 14149. , 13910. , 13670. , 13430. , 13190. ,
       12950. , 12710.5, 12471. , 12231. , 11991. , 11751. , 11511. ,
       11271.5, 11032. , 10792. , 10552. , 10312. , 10072. ,  9833. ,
        9593. ,  9353. ,  9113. ,  8873. ,  8633. ,  8394. ,  8154. ,
        7914. ,  7674. ,  7434. ,  7194. ,  6955. ,  6715. ,  6475. ,
        6235. ,  5995. ,  5755. ,  5516. ,  5276. ,  5036. ,  4796. ,
        4556. ,  4316. ,  4077. ,  3837. ,  3597. ,  3357. ,  3117. ,
        2877.5,  2638. ,  2398. ,  2158. ,  1918. ,  1678. ,  1438.5,
        1199. ,   959. ,   719. ,   479. ,   239. ,     0. ,  -240. ,
        -480. ])
    
    height_mins = np.array([17867., 17627., 17387., 17147., 16907., 16668., 16428., 16188.,
       15948., 15708., 15468., 15229., 14989., 14749., 14509., 14269.,
       14029., 13790., 13550., 13310., 13070., 12830., 12591., 12351.,
       12111., 11871., 11631., 11391., 11152., 10912., 10672., 10432.,
       10192.,  9952.,  9713.,  9473.,  9233.,  8993.,  8753.,  8513.,
        8274.,  8034.,  7794.,  7554.,  7314.,  7074.,  6835.,  6595.,
        6355.,  6115.,  5875.,  5635.,  5396.,  5156.,  4916.,  4676.,
        4436.,  4196.,  3957.,  3717.,  3477.,  3237.,  2997.,  2758.,
        2518.,  2278.,  2038.,  1798.,  1558.,  1319.,  1079.,   839.,
         599.,   359.,   119.,  -120.,  -360.,  -600.])
    
    height_maxes = np.array([18107., 17867., 17627., 17387., 17147., 16908., 16668., 16428.,
       16188., 15948., 15708., 15469., 15229., 14989., 14749., 14509.,
       14269., 14030., 13790., 13550., 13310., 13070., 12830., 12591.,
       12351., 12111., 11871., 11631., 11391., 11152., 10912., 10672.,
       10432., 10192.,  9953.,  9713.,  9473.,  9233.,  8993.,  8753.,
        8514.,  8274.,  8034.,  7794.,  7554.,  7314.,  7075.,  6835.,
        6595.,  6355.,  6115.,  5875.,  5636.,  5396.,  5156.,  4916.,
        4676.,  4436.,  4197.,  3957.,  3717.,  3477.,  3237.,  2997.,
        2758.,  2518.,  2278.,  2038.,  1798.,  1558.,  1319.,  1079.,
         839.,   599.,   359.,   120.,  -120.,  -360.])

    # summit-related (v6-added) options/variables
    SUMMIT_MODE = None #off by default
    month=None #month of granules to use
    year=None #year of granules to use
    dayfi=None #first day of granules to use
    dayla=None #last day of granules to use
    outdir=os.getcwd() #where to save output files (default=current dir)
    scratch_dir=parent_dir #where to save/read intermediates
    local_dir=outdir

    # defaults activates when SUMMIT_MODE==True
    parent_dir_summit='/pl/active/wbclimate/ftp.cloudsat.cira.colostate.edu'
    radar_folder_summit='2B-GEOPROF.P1_R05'
    lidar_folder_summit='2B-GEOPROF-LIDAR.P2_R05'
    scratch_dir_summit='/scratch/alpine/wibe4964/cloudsat_intermediates'#'/scratch/summit/wibe4964' 

    ### SET OPTIONS FROM CLI ARGS ###
    argv = sys.argv[1:] #TODO NOTE: UNCOMMENT THIS
    print("running main()")
    print(" ".join(sys.argv[:]))
    try:
        opts,args=getopt.getopt(argv,"hp:r:l:ck:o:m:n:g:f:t:v:s",
                                ['month=','year=','dayfi=','dayla=','localdir=','help',
                                 'outdir=','print-defaults','na','nb','nd','def='])
    except getopt.GetoptError:
        print('invalid options(s), get help with flag -h')
        sys.exit(1)
    args = [o[0] for o in opts]
    for opt,arg in opts:
        print(f'{opt}:{arg}')
        if opt == '-p':
            parent_dir = arg
        elif opt=='-r':
            radar_folder=arg
        elif opt=='-l':
            lidar_folder=arg
        elif opt=='-b':
            build_intermediates=(arg=='True')            
        elif opt=='-d':
            delete_intermediates=(arg=='True')
        elif opt=='-c':
            check_intermediates=arg
        elif opt=='-o':
            include_only=arg
        elif opt=='-m':
            lidar_thresh=float(arg)
        elif opt=='-n':
            radar_thresh=int(arg)
        elif opt=='-g':
            grid_spacing=float(arg)
        elif opt=='-f':
            fname_suffix=arg
        elif opt=='-k':
            keep_granule=arg
        elif opt=='-t':
            num_printouts=int(arg)
        elif opt=='-v':
            if arg in ('height','cover','both'):
                vertically_resolved=arg
            else:
                print(f'unrecognized output mode {arg}! Run -h for list of recognized options. Quitting...')
                sys.exit(1)
        elif opt=='-s':
            SUMMIT_MODE=True
            parent_dir = parent_dir_summit
            radar_folder = radar_folder_summit
            lidar_folder = lidar_folder_summit
            scratch_dir = scratch_dir_summit
        elif opt=='--month':
            month=arg
        elif opt=='--year':
            year=arg
        elif opt=='--dayfi':
            dayfi=int(arg)
        elif opt=='--dayla':
            dayla=int(arg)
        elif opt=='--outdir':
            outdir=arg
        elif opt=='--localdir':
            local_dir=arg
        elif opt=='--scratchdir':
            scratch_dir=arg
        elif opt=='--na':
            attenuate_lidar=False
        elif opt=='--nb':
            build_intermediates=False
        elif opt=='--nd':
            delete_intermediates=False
        elif opt=='--def':
            #NB: a little ugly since minN is variable
            if arg in recognized_defs:
                cloudy_def=arg
            else:
                try:
                    if arg[:3] == 'min' and int(arg[3:]):
                        cloudy_def=arg
                    else:
                        print(f'unrecognized cloudy def {arg}! Run -h for list of recognized options. Quitting...')
                        sys.exit(1)
                except ValueError:
                    print(f'unrecognized cloudy def {arg}! Run -h for list of recognized options. Quitting...')
                    sys.exit(1)
        elif opt=='-h' or opt=='--help':
            print(main.__doc__)
            sys.exit(0)

    print(f'local directory: {local_dir}')
    print('initializing cluster...',end='')
    ### INITIALIZE CLUSTER ###
    n_workers = 32 #32 
    threads_per_worker = 2 #2 
    #TODO: set local_directory=scratch_dir where --scratch-dir=$SLURM_SCRATCH
    cluster = LocalCluster(n_workers=n_workers, 
            threads_per_worker=threads_per_worker,local_directory=local_dir)
    #    cluster = LocalCluster() #let dask infer defaults while learning summit
    client = Client(cluster)
    print('done')
    print_cluster_info(client)
    print(f'View the cluster at {client.dashboard_link}')


    ### ENUMERATE INPUT FILE PATHS AND ###
    ### GET INTERMEDIATE FOLDER NAME ###
    print('processing options and enumerating filepaths...',end='')
    if SUMMIT_MODE:
        #check for proper input conditions
        if month and (dayfi or dayla):
            print('ERROR: cannot specify both a month and a day! Quitting...')
            sys.exit(1)
        if (dayfi or dayla) and ((dayfi is None) or (dayla is None)):
            print('ERROR: both --dayfi and --dayla must be specified if not using --month! Quitting...')
            sys.exit(1)
        if not year:
            print('ERROR: --year must be specified when using summit mode (-s)! Quitting...')
            sys.exit(1)

        #process time args
        if len(year)==1:
            year='200'+year #add century+decade
        elif len(year)==2:
            year='20'+year #add century
        if dayfi and dayla:
            dayfi,dayla=int(dayfi),int(dayla)
        if month:
            #get first/last day of year for the month
            t=datetime.strptime(f'{year}{month}','%Y%m')
            dayfi=int(t.strftime('%j'))
            dayla=int(last_day_of_month(t).strftime('%j'))

        #get correct folders to process + get items in folders
        radar_filepaths = get_granule_paths(dayfi,dayla,os.path.join(parent_dir,radar_folder,year))
        lidar_filepaths = get_granule_paths(dayfi,dayla,os.path.join(parent_dir,lidar_folder,year))

        #folder name for intermediates
        t0=datetime.strptime(f'{dayfi:03d}{year}','%j%Y')
        t1=datetime.strptime(f'{dayla:03d}{year}','%j%Y')
        #if dates equal whole month use e.g. 2013jan,
        #otherwise includes months+dates
        if (t0.day==1) and (t1.day==last_day_of_month(t0).day):
            timestr=t0.strftime('%Y-%m').lower()
        else:
            timestr=(t0.strftime('%Y-%m-%d')+'_'+t1.strftime('%m-%d')).lower()
        new_combined_folder = f"nc_{timestr}_combined"
        new_combined_path = os.path.join(scratch_dir,new_combined_folder)

    elif not SUMMIT_MODE:
        #items in folder
        radar_filepaths = tuple(listdir_fullpaths(os.path.join(parent_dir,radar_folder)))
        lidar_filepaths = tuple(listdir_fullpaths(os.path.join(parent_dir,lidar_folder)))

        #folder name for intermediates
        new_combined_folder = 'nc_'+radar_folder+'_combined'
        new_combined_path = os.path.join(parent_dir,new_combined_folder)

    print('done')
    ### MAKE FOLDER FOR SAVING INTERMEDIATES ###
    try:
        os.makedirs(new_combined_path, exist_ok=False)
        print(f'made directory {new_combined_path} to store intermediate files')
    except FileExistsError:
        print(f'intermediates directory {new_combined_path} already present!')
        if build_intermediates:
            print(f'But build_intermediates == {build_intermediates}... Quitting')
            sys.exit(1)

    ### BUILD INTERMEDIATE NETCDF4 FILES FROM HDF4-EOS GRANULES ###
    keep_vars=['CPR_Cloud_mask','SurfaceHeightBin','DEM_elevation','Profile_time', #NOTE: changed v7.2!
               'Data_quality','CloudFraction','Latitude','Longitude','Data_status'] 
    #only keep granules for which there is both lidar and radar data
    #always compute since it's needed for the quicklook printout()
    thin_radar_filepaths, thin_lidar_filepaths = intersection(radar_filepaths,lidar_filepaths,keep_granule=keep_granule)
    if build_intermediates:
        print(f'writing {len(thin_radar_filepaths)} intermediate netcdf4 files at {new_combined_path}...')
        #use dask.bag to parallelize the iteration
        rlb = db.from_sequence(zip(thin_radar_filepaths,thin_lidar_filepaths))
        _   = rlb.map(read_convert_merge_granule,new_combined_path,keep_only=keep_vars).compute()
        print('write complete. Restarting dask client...')
    #        client.restart()

    # if I add `del` statements at the end of `read_convert_merge_granule` will the garbage collection issues improve? <br>

    #NOTE DEV: CHANGED FOR V8.1

    print('opening intermediates as dataset...',end='')
    ### LAZILY LOAD NC4S FROM DISK ###
    tcv = xr.open_mfdataset(new_combined_path+'/*.nc',concat_dim='start_time',combine='nested',
                             decode_cf=True, mask_and_scale=False, decode_times=True, 
                             decode_timedelta=True, use_cftime=False, #modified decode_timedelta v7.2
                             concat_characters=False, decode_coords=True, parallel=True)
    print('done')

    ### CHECK INTERMEDIATE NC4S ###
    if check_intermediate and (build_intermediates==False):
        thin_radar_filepaths, thin_lidar_filepaths = intersection(radar_filepaths,lidar_filepaths,quiet=True)
    if check_intermediate:
        check_i = np.random.randint(len(thin_radar_filepaths))
        check_granule = fpath2granule(thin_radar_filepaths[check_i])
        print(f'checking intermediate vs original for granule {check_granule}...')

        #check that all variables in the original 2B-GEOPROF == nc4 intermediate
        cpr  = read_GEOPROF(thin_radar_filepaths[check_i],nprofiles=37100)
        vfm  = read_GEOPROF(thin_lidar_filepaths[check_i],nprofiles=37100)
        cv   = xr.merge([cpr,vfm],compat='no_conflicts')
        cv2  = tcv.isel(start_time=check_i).compute()
        keys = [*cv2.keys(), *cv2.coords]

        print('radar comparison: ')
        for key in keys:
            comp = compare(cv[key].data,cv2[key].data)
            print('\t {} equal: {}'.format(key,comp))

    #NEW V8.2.1
    ### CHECK IF ANY PROFILES HAVE LONGITUDE==+180 ###
    lon_mask = (tcv.Longitude==180.0).compute()
    if lon_mask.any():
	# subtract a tiny value (equal to a max of 0.1 um) to fix the gridbox assignment
        print(f'correcting {np.sum(lon_mask.values)} profiles with Longitude==180')
        tcv = tcv.assign(Longitude=tcv.Longitude-(1e-12)*lon_mask.astype(int))

    #NEW V7.2
    ### ADD/HANDLE QC VARIABLES ###
    #add in sampling fields
    #days in period (input DayOfYear)
    tcv = tcv.assign(dofy=(tcv.start_time+tcv.Profile_time).dt.dayofyear)
    #unique granules (input start_time)
    tcv = tcv.assign(overpass=tcv.start_time)
    #local time bins
    shape, dims = tcv.Longitude.data.shape,tcv.Longitude.dims
    utc_times_flat = (tcv.start_time+tcv.Profile_time).data.flatten()
    lon_flat = tcv.Longitude.data.flatten()
    tz_offsets_flat = pd.to_timedelta((lon_flat+15/2)//15,unit='h')
    local_times_flat = pd.to_datetime(utc_times_flat + tz_offsets_flat)
    local_hour_flat = np.asarray(local_times_flat.hour) #can't reshape pd.DateTimeIndex
    local_hour = local_hour_flat.reshape(shape)
    # next step: bool for 0<=t<6, 6<=t<12, 12<=t<18, 18<=t<24
    lt_funcs  = {'localhour22':lambda x: np.logical_or(x>=22,x<4),
                 'localhour04':lambda x: np.logical_and(x>=4,x<10),
                 'localhour10':lambda x: np.logical_and(x>=10,x<16),
                 'localhour16':lambda x: np.logical_and(x>=16,x<22)}
    new_vars = {name: (dims,func(local_hour)) for name,func in lt_funcs.items()}
    tcv = tcv.assign(new_vars)
    del shape,dims,utc_times_flat,lon_flat,tz_offsets_flat
    del local_times_flat,local_hour_flat,local_hour,new_vars

    #handle radar data quality
    print_dataquality(tcv) #print out diagnostic info
    tcv = tcv.where(tcv.Data_quality==0) #if any DQ flag is set, set all vars in the profile to NaN
    #note -- this might have to go before or after the sampling vars
    
    #restrict valid profiles to time range
    print()
    print('removing profiles not in month...')
    dstime=tcv.start_time+tcv.Profile_time
    time_mask = (dstime>=pd.to_datetime(t0)) & (dstime<=pd.to_datetime(t1)+pd.Timedelta(1,unit='days'))
    tcv = tcv.where(time_mask,drop=False)

    #print out info on what was masked
    tfmt = lambda da: da.start_time.dt.strftime("%Y-%m-%d %H:%M:%S").values
    nmski = int((~time_mask.isel(start_time=0)).sum(dim='profile') - np.isnan(dstime.isel(start_time=0)).sum(dim='profile'))
    nmskf = int((~time_mask.isel(start_time=-1)).sum(dim='profile') - np.isnan(dstime.isel(start_time=-1)).sum(dim='profile'))
    print(f"masked {nmski} profs (~{nmski*0.16/60:.1f} min) in first granule with start_time {tfmt(time_mask.isel(start_time=0))}")
    print(f"masked {nmskf} profs (~{nmskf*0.16/60:.1f} min) in last granule with start_time {tfmt(time_mask.isel(start_time=-1))}")

    ### ADD ATTRIBUTES IF DUMMY GRANULES ARE PRESENT ###
    # new v7.2.1: this fix should really be in read_convert_merge_granule's dummy features,
    # but I'm putting it here since I don't want to rewrite all the intermediates.
    if keep_granule == 'radar': #if dummy lidar granules are present,
        #then add back the correct attributes (since they may have been stripped)
        tcv['CloudFraction'].attrs = {'_FillValue':-9.0,
                                      'factor':1.0,
                                      'offset':0.0,
                                      'long_name':'Cloud Fraction',
                                      'valid_range':[0,100],
                                      'missing':-9,
                                      'missop':'=='}
    elif keep_granule == 'lidar': #if dummy radar granules are present,
        #then add back the correct attributes (since they may have been stripped)
        tcv['CPR_Cloud_mask'].attrs = {'_FillValue':-9.0,
                                      'factor':1.0,
                                      'offset':0.0,
                                      'long_name':'Cloud Fraction',
                                      'valid_range':[0,100],
                                      'missing':-9,
                                      'missop':'=='}

    ### ASSIGN GRIDBOX NUMBER TO EACH START_TIME,PROFILE COORD ###
    tcv_gb = ds_ll2gbnum(tcv,dx=grid_spacing)
    tcv_gb = tcv_gb.drop_vars(['Latitude','Longitude'])

    ### DECODE/PROCESS RADAR AND LIDAR TO BINARY MASKS ###
    # make a template dataset containing the variables proc_and_mask will return
    template = tcv_gb.rename({'CPR_Cloud_mask':'cpr_binary','CloudFraction':'vfm_binary'})
    template = template.assign(radar_surface_clutter_counts=template.cpr_binary.astype(bool)) #added v8.1 -- return surface clutter occurrence
    if attenuate_lidar: 
        template = template.assign(attenuated_lidar_counts=template.vfm_binary.astype(bool))
    tcm = tcv_gb.map_blocks(proc_and_mask,template=template,kwargs={'lidar_thresh':lidar_thresh,'radar_thresh':radar_thresh})
    tcm = tcm.drop_vars(('DEM_elevation','Data_quality','SurfaceHeightBin')) #modified for v7.2
    # variables needed for final output
    # NEW TO V7.2 since del statements were added
    start_time_initial = tcm.start_time.isel(start_time=0)
    start_time_final = tcm.start_time.isel(start_time=-1)
    n_granules = tcm.dims['start_time']

    ### MERGE MASKS TO TOTAL_COUNTS, CLOUD_COUNTS ###
    if include_only == 'both':
        tcm = get_total_counts_nosum(tcm)
        tcm = get_cloud_counts_nosum(tcm)
    elif include_only == 'lidar':
        tcm = get_total_counts_nosum_lidaronly(tcm)
        tcm = get_cloud_counts_nosum_lidaronly(tcm)
    elif include_only == 'radar':
        tcm = get_total_counts_nosum_radaronly(tcm)
        tcm = get_cloud_counts_nosum_radaronly(tcm)
    tcc = tcm.drop_vars(['vfm_binary','cpr_binary'],errors='ignore')
    tcc = tcc.persist() #added v7.2


    #TODO: filter month by profile here
    #tcc = tcc.drop_vars(['Profile_time'])

    ### PROCESS PROFILES TO BE VERTICALLY-RESOLVED AND/OR COLUMN-INTEGRATED ###
    #remove keys that are not 'cloud_counts' or 'total_counts' from the inner dicts using dict comprehension (py>2.7)
    #do it since other vars like 'cloud_cover' are not present yet and would cause xr.rename to fail
    rename_counts = {k:{i:d[i] for i in d if i in ('cloud_counts',
                                                   'total_counts',
                                                   'attenuated_lidar_counts',
                                                   'radar_surface_clutter_counts')} for k,d in rename.items()}

    #exclude non-counts variables from this list (e.g. doesn't end in "_counts" (new v7.2)
    non_counts_names = list(filter(lambda x: x.split('_')[-1]!='counts',list(tcc.data_vars)))
    s_tcc = tcc.drop_vars(non_counts_names) #s for counts which get summed
    u_tcc = tcc[non_counts_names] #u for non-counts which get the number of unique values

    print(f"computing output fields '{vertically_resolved}' with cloud definition {cloudy_def}...",end='')
    #compute+rename 3d occurrence only
    if vertically_resolved == 'height':
        stcm  = s_tcc.sel(bin=slice(29,107))#.compute() removed memtest1
        stcm  = stcm.rename(rename_counts['height'])

    #compute+rename cloud cover definition only
    elif vertically_resolved == 'cover':
        stcm  = s_tcc.map_blocks(apply_cloud_def,kwargs={'alg':cloudy_def})#.compute() removed memtest1
        stcm  = stcm.rename(rename_counts['cover'])
     #note v8.0: this feature is broken now, since cloudy_def='all' won't be recognized


    #compute+rename both and combine
    #feed one into the other for speedup
    elif vertically_resolved == 'both':
        if cloudy_def=='all':
            stcm  = apply_all_cloud_defs(s_tcc) #this is where radar_surface_clutter_counts_in_column would get computed but doesn't
            hstcm = s_tcc.rename(rename_counts['height']).sel(bin=slice(29,107))
            stcm  = stcm.rename({var:f'{var}_in_column' for var in stcm.data_vars}) #add _in_column to every var
            stcm  = xr.merge([hstcm,stcm])
        else:
            #apply_cloud_def modifies the ds in-place now. .any in util method needs to return .any result into a new array
            stcm  = s_tcc.map_blocks(apply_cloud_def,kwargs={'alg':cloudy_def})#.compute() removed memtest1
            hstcm = s_tcc.rename(rename_counts['height']).sel(bin=slice(29,107))
            stcm  = stcm.rename(rename_counts['cover'])
            stcm  = xr.merge([hstcm,stcm])


    ### COMPUTE BEFORE GROUPY TO AVOID NUMBER OF TASKS OVERFLOW ###
    print('done. Computing mask/process/gridbox_num assign...',end='')
    #u_tcc = u_tcc.compute()
    stcm = xr.merge([stcm,u_tcc]) #put back non-counts variables
    stcm = stcm.compute() # ~7 GB
    #del tcv,tcm,tcc,s_tcc,u_tcc #added for v7.2 (comment out to enable quicklooks)
    #stcm = stcm.compute() 
    print('done')

    ### CREATE TOP-LEVEL GROUPS BY DOOP-FLAG ###
    # v7.2: changed from doing groupby doop and summing because make unique fail (only works for sum)
    # NB: this makes swap usage for 'All cases' spike to ~45 GB during groupby
    # NB: rechunking before groupby slows down operation 6x
    stacked_stcm = stcm.stack(stacked_start_time_profile=('start_time','profile'))
    doop_stcm = stacked_stcm.sel(stacked_start_time_profile=stacked_stcm.doop_flag==True)
    d = {b'DO-OP observable':doop_stcm,
         b'All cases':stacked_stcm}
    d_new = { }
    for key, dstcm in d.items():
        print(f'doing gridding for doop_flag=={key}')
        d_new[key] = do_gridding(dstcm,grid_spacing,height_midpoints)

    ### JOIN RUNS INTO SINGLE DS WITH NEW DIM DOOP ###
    gridded = xr.concat([ds.assign_coords(doop=key) for key,ds in d_new.items()],dim='doop')
    
    ### PUT CLOUD COVER BY TYPE ALONG NEW DIMENSION ###
    # apply_all_cloud_defs returns 8 total_counts_in_column and cloud_counts_in_column by cover type.
    # reorganize these to be in two data variables but with a new dimension 'type'
    if cloudy_def == 'all':
        da_cc = xr.concat([gridded['cloud_counts_'+alg+'_in_column'].assign_coords({'type':alg}) for alg in standard_defs],dim='type')
        gridded = gridded.drop_vars(['cloud_counts_'+alg+'_in_column' for alg in standard_defs])
        gridded = gridded.assign(cloud_counts_in_column=da_cc)

        da_tc = xr.concat([gridded['total_counts_'+alg+'_in_column'].assign_coords({'type':alg}) for alg in standard_defs],dim='type')
        gridded = gridded.drop_vars(['total_counts_'+alg+'_in_column' for alg in standard_defs])
        gridded = gridded.assign(total_counts_in_column=da_tc)

    ### CALCULATE CLOUD COVER ###
    # modified for v7.1:
    # assign 'cloud_cover' variables and names according to specified -v mode and the 'rename' dict
    if vertically_resolved == 'both': 
        for kind in ('height','cover'):
            cloud_counts = gridded[rename[kind]['cloud_counts']]
            total_counts = gridded[rename[kind]['total_counts']]
            gridded = gridded.assign({rename[kind]['cloud_cover']:cloud_counts/total_counts})
    else: #covers both 'height' and 'cover'
        kind = vertically_resolved
        cloud_counts = gridded[rename[kind]['cloud_counts']]
        total_counts = gridded[rename[kind]['total_counts']]
        gridded = gridded.assign({rename[kind]['cloud_cover']:cloud_counts/total_counts})

    # output is done now, so delete other variables
    del total_counts, cloud_counts, d_new, doop_stcm, stacked_stcm


    #NOTE DEV V8.1: CHANGED TO INCLUDE EXTRA DATA VARIABLE
    ### CALCULATE FILENAME ###
    #date_now = datetime.today().strftime('%b%dT%H%M') #unused variable
    try:
        timestr #timestr is undefined outside SUMMIT_MODE
    except NameError:
        timestr = str(start_time_initial.dt.strftime('%Y-%m').data)

    #indicate cloud cover only, 3d occurrence, or cover+3d occurrence in filename ('both' is too ambiguous)
    covertype = 'coverandheight' if vertically_resolved=='both' else vertically_resolved
    
    #new v8.2
    #get flavor from args (e.g. -RO, -LO, -RO2)
    flavor = ''
    if keep_granule == 'both':
        if include_only == 'radar':
            flavor = '-RO'
        elif include_only == 'lidar':
            flavor = '-LO'
    elif keep_granule == 'radar':
        if include_only == 'radar':
            flavor = '-RO2'
    fname = f'{timestr}_CSCAL_3S-GEOPROF-COMB{flavor}_{grid_spacing}x{grid_spacing}_{fname_suffix}'

    ### PRINT OUT QUICKLOOKS FOR SANITY CHECKING ###
    print(f'starting {num_printouts} quicklook printouts... ',end='')
    for i in range(num_printouts):
        n = np.random.randint(tcv.dims['start_time'])
        #printout(n,thin_radar_filepaths,tcv,tcm,tcc,fname,outdir) #removed v7.2-memtest
    print('done')

    ### SAVE PLOT FOR TOTAL CLOUDFRACTION ###
    fig = plt.figure(num=None, figsize=(10, 5), dpi=150, edgecolor='k')
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    doopstate=b'All cases'
    data = gridded.sel(doop=doopstate)[rename['cover']['cloud_cover']]
    if cloudy_def == 'all': data = data.sel(type='any')
    im = ax.imshow(data, origin='lower', extent=[-180,180,-90,90], transform=ccrs.PlateCarree(),cmap='viridis',interpolation='none',vmin=0,vmax=1)
    ax.coastlines()

    t0  = str(start_time_initial.dt.strftime('%Y-%m-%d').data)
    t1  = str(start_time_final.start_time.dt.strftime('%Y-%m-%d').data)
    title = f'{t0} to {t1} Cloudsat+Calipso Cloud Cover doop={doopstate}, def={cloudy_def}, {grid_spacing}$^{{\circ}}$ grid ({n_granules} granules)'
    ax.set_title(title)
    cb=fig.colorbar(im,fraction=0.023, pad=0.04)
    cb.set_label('cloud cover fraction')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(outdir,fname+'.png'),transparent=False,facecolor='white',bbox_inches="tight")

    ### SET GLOBAL ATTRIBUTES ###
    global_attrs = { }

    #time created
    tstr = datetime.today().strftime('%Y-%m-%d %H:%M:%S MST')
    global_attrs['date_created'] = f'{tstr}'

    #created by script
    global_attrs['from_script'] = f'file created by script {os.path.basename(__file__)}'
    global_attrs['from_command'] = " ".join(sys.argv[:])

    #other
    global_attrs['author'] = f'William J Bertrand'
    global_attrs['num_granules'] = n_granules
    global_attrs['keep_type'] = 'calculated only from granules for which both radar and lidar data are available'
    global_attrs['resolution_lon (deg)'] = grid_spacing   
    global_attrs['resolution_lat (deg)'] = grid_spacing
    global_attrs['cloud_type'] = cloudy_def
    global_attrs['radar_thresh'] = radar_thresh
    global_attrs['lidar_thresh'] = lidar_thresh

    #date range
    ti = str(start_time_initial.dt.strftime('%Y-%m-%d %H:%M:%S').data)+' UTC'
    tf = str(start_time_final.dt.strftime('%Y-%m-%d %H:%M:%S').data)+' UTC'
    global_attrs['time_range'] = f'{ti} - {tf}'

    #apply global attrs
    gridded.attrs = global_attrs

    ### SET COORDINATE AND DATA VAR ATTRIBUTES ###
    #specify all coordinate and data var attributes here in the form attribute:varname:val
    #v7.2: added attenuated_lidar_counts_in_column, attenuated_lidar_counts_on_levels, n_overpasses
    #, n_days, localhour22, localhour04, localhour10, localhour16
    descriptions = {rename['height']['total_counts']: 'total number of valid radar-lidar samples (clear sky or cloudy) in lat/lon/height bin',
                    rename['height']['cloud_counts']: 'number of cloudy radar-lidar samples in lat/lon/height bin',
                    rename['height']['cloud_cover']: 'relative frequency (0-1) of cloud in lat/lon/height bin',
                    rename['height']['attenuated_lidar_counts']: 'number of estimated attenuated lidar counts in lat/lon/height bin',
                    rename['height']['radar_surface_clutter_counts']: 'number of times surface clutter prevented radar measurement in lat/lon/height bin',
                    rename['cover']['total_counts']: 'total number of valid radar-lidar profiles (clear sky or cloudy) in lat/lon grid cell',
                    rename['cover']['cloud_counts']: 'number of cloudy profiles in lat/lon grid cell',
                    rename['cover']['cloud_cover']: 'relative frequency (0-1) of cloudy columns in lat/lon grid cell',
                    rename['cover']['attenuated_lidar_counts']: 'profiles which attenuate at some level',
                    'type':'''Different cloudy criteria applied to profiles. High, middle, low are satisfied if any cloud is present in a height range.
    Unique variants are satisfied if all cloud in the profile is within the height range. Thick is satisfied if a cloud layer is >=4.8 km thick.
    Any is satisfied if any cloud is in the profile. Total counts for high/middle/low report the number of profiles with valid observations in the height range.
    Total counts for unique variants report the number of profiles with valid observations in all three height ranges.''',
                    'n_overpasses': 'number of unique visits to/overpasses of the grid cell during the aggregating period',
                    'n_days': 'number of unique days in month grid cell was observed',
                    'localhour22': 'number of profiles with local time >=2200 and <4000',
                    'localhour04': 'number of profiles with local time >=4000 and <1000',
                    'localhour10': 'number of profiles with local time >=1000 and <1600',
                    'localhour16': 'number of profiles with local time >=1600 and <2200'}

    units = {rename['height']['total_counts']: 'number', #or units "count" or blank ""
             rename['height']['cloud_counts']: 'number',
             rename['height']['cloud_cover']: 'ratio',
             rename['height']['attenuated_lidar_counts']: 'number',
             rename['height']['radar_surface_clutter_counts']: 'number',
             rename['cover']['total_counts']: 'number',
             rename['cover']['cloud_counts']: 'number',
             rename['cover']['cloud_cover']: 'ratio',
             rename['cover']['attenuated_lidar_counts']: 'number',
             'lat': 'degrees north',
             'lon': 'degrees east',
             'height': 'm',
             'height_min': 'm',
             'height_max': 'm',
             'n_overpasses':'number',
             'n_days':'number',
             'localhour22': 'number',
             'localhour04': 'number',
             'localhour10': 'number',
             'localhour16': 'number'}

    long_name = {rename['height']['total_counts']: 'Number of observations at altitude',
                 rename['height']['cloud_counts']: 'Number of observations of cloud at altitude',
                 rename['height']['cloud_cover']: 'Cloud counts at level divided by total counts at level',
                 rename['height']['attenuated_lidar_counts']: 'Estimated count of attenuated lidar samples at level',
                 rename['height']['radar_surface_clutter_counts']: 'Count of no radar observation due to surface clutter at level',
                 rename['cover']['total_counts']: 'Number of observations of whole column',
                 rename['cover']['cloud_counts']: 'Number of times cloud of (type) present in column',
                 rename['cover']['cloud_cover']: 'Cloud counts in column of type divided by total counts of type in column',
                 rename['cover']['attenuated_lidar_counts']: 'Count of lidar attenuation estimated anywhere in column',
                 'lat': 'latitude midpoint of grid cell',
                 'lon': 'longitude midpoint of grid cell',
                 'height': 'altitude midpoint meters above mean sea level',
                 'height_min': 'minimum height of bin range (above mean sea level)',
                 'height_max': 'maximum height of bin range (above mean sea level)',
                 'doop':'Daylight-only operations status',
                 'type':'cloud cover definition',
                 'n_overpasses':'number of overpasses of grid cell',
                 'n_days':'number of unique days grid cell observed',
                 'localhour22': 'number of profiles in local time bin',
                 'localhour04': 'number of profiles in local time bin',
                 'localhour10': 'number of profiles in local time bin',
                 'localhour16': 'number of profiles in local time bin'}

    dx = {'lat': grid_spacing, 'lon': grid_spacing}

    #put all attrs together with top level as attr name
    attrs_by_attrname = {'dx':dx,'description':descriptions,'units':units,'long_name':long_name} 

    #rename non-counts variables
    gridded = gridded.rename({'overpass':'n_overpasses','dofy':'n_days'})
    
    #new v8.2: add heights corresponding to bin range maxima and minima (since it's not an even 240 m)
    gridded = gridded.assign({'height_min':height_mins, 'height_max':height_maxes})

    #transpose from attribute:varname:val to varname:attribute:val
    full_vars = [*list(gridded.data_vars),*list(gridded.coords)] #flat list of vars and coords
    attrs_by_varname = {k:{} for k in full_vars} #build nested dict with right keys
    for attrname,d in attrs_by_attrname.items():
        for varname,val in d.items():
            attrs_by_varname[varname][attrname] = val

    #apply attrs to dataset
    for varname in full_vars:
        gridded[varname].attrs = attrs_by_varname[varname]

    #added v7.2
    #re-order data variables in dataset so the listing makes more sense
    sensical_ordering = (rename['height']['cloud_counts'],
                         rename['height']['total_counts'],
                         rename['height']['cloud_cover'],
                         rename['cover']['cloud_counts'],
                         rename['cover']['total_counts'],
                         rename['cover']['cloud_cover'],
                         rename['height']['attenuated_lidar_counts'],
                         rename['cover']['attenuated_lidar_counts'],
                         rename['height']['radar_surface_clutter_counts'],
                         'n_overpasses',
                         'n_days',
                         'localhour22',
                         'localhour04',
                         'localhour10',
                         'localhour16')
    reordered_vars = {name:gridded[name] for name in sensical_ordering}
    gridded = xr.Dataset(data_vars = reordered_vars, coords=gridded.coords, attrs=gridded.attrs)
    

    ### EXPORT GRIDDED DATASET TO DISK ###
    print(f'saving file {fname}.nc')
    gridded.to_netcdf(os.path.join(outdir,fname+'.nc'))


    ### DONE WITH DASK CLUSTER ###
    client.shutdown()

    ### DELETE INTERMEDIATE FILES SAVED TO DISK ###
    if delete_intermediates:
        print('deleting intermediate files...')
        combined_granules = list(listdir_fullpaths(new_combined_path))
        if all([path[-3:] == '.nc' for path in combined_granules]):
            print('all non-hidden files in {} ending in .nc... \ndeleting {} files...'
                  .format(new_combined_path,len(combined_granules)),end='')
            for granule in combined_granules:
                os.remove(granule)
                #print(granule)
            print('done.')
            hidden_files = os.listdir(new_combined_path)
            hidden_filepaths = [os.path.join(new_combined_path,hidden_file) for hidden_file in hidden_files]
            try:
                for hidden_filepath in hidden_filepaths:
                    print('deleting hidden file {}'.format(hidden_filepath))
                    os.remove(hidden_filepath)
                os.rmdir(new_combined_path)
                print(f'deleted directory {new_combined_path}')
            except OSError:
                print(f'ERROR: OS Error when deleting hidden files! Directory {os.path.split(new_combined_path)[-1]} not deleted')
                sys.exit(1)
        else:
            print(f'not all files in {new_combined_path} have .nc file extensions. Are you sure you want to delete these files?')
    else:
        print('delete_intermediates set to False. Leaving intermediate files in place...')

    print('success')
    sys.exit(0)

if __name__ == "__main__":
    main()

