#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#Created on Sat Apr 20 22:46:44 2019

#@author: Chang-Eop Kim
"""

import numpy as np
import matplotlib.pyplot as plt
import os
# os.chdir('/Users/seroh/Desktop/ELIFE/CFPC_R4_2/') # type "os.getcwd()" to check the folder
fnames = 'file_name.tif'


import caiman as cm
from caiman.motion_correction import MotionCorrect
from caiman.source_extraction.cnmf import cnmf as cnmf
from caiman.source_extraction.cnmf import params as params
from caiman.utils.visualization import plot_contours, nb_view_patches, nb_plot_contour


    
##Setup some parameters
# dataset dependent parameters
fr = 32                             # imaging rate in frames per second
decay_time = 0.4                    # length of a typical transient in seconds

# motion correction parameters
strides = (48, 48)          # start a new patch for pw-rigid motion correction every x pixels
overlaps = (24, 24)         # overlap between pathes (size of patch strides+overlaps)
max_shifts = (6,6)          # maximum allowed rigid shifts (in pixels)
max_deviation_rigid = 3     # maximum shifts deviation allowed for patch with respect to rigid shifts
pw_rigid = True             # flag for performing non-rigid motion correction

# parameters for source extraction and deconvolution
p = 1                       # order of the autoregressive system
gnb = 2 #2                     # number of global background components
merge_thr = 0.85 #0.85         # merging threshold, max correlation allowed
rf = 80 #15   250                # half-size of the patches in pixels. e.g., if rf=25, patches are 50x50
stride_cnmf = 3             # amount of overlap between the patches in pixels
K = 6 #4                      # number of components per patch
gSig = [20,20] #[4, 4]               # expected half size of neurons in pixels
method_init = 'sparse_nmf'#'greedy_roi'  # initialization method (if analyzing dendritic data using 'sparse_nmf')
ssub = 1                    # spatial subsampling during initialization
tsub = 1                    # temporal subsampling during intialization

# parameters for component evaluation
min_SNR = 2             # signal to noise ratio for accepting a component
rval_thr = 0.7           # space correlation threshold for accepting a component
cnn_thr = 0  #0.99              # threshold for CNN based classifier
cnn_lowest = 0 #0.1 # neurons with cnn probability lower than this value are rejected

# Create a parameters object
opts_dict = {'fnames': fnames,
            'fr': fr,
            'decay_time': decay_time,
            'strides': strides,
            'overlaps': overlaps,
            'max_shifts': max_shifts,
            'max_deviation_rigid': max_deviation_rigid,
            'pw_rigid': pw_rigid,
            'p': 1,
            'nb': gnb,
            'rf': rf,
            'K': K, 
            'stride': stride_cnmf,
            'method_init': method_init,
            'rolling_sum': True,
            'only_init': True,
            'ssub': ssub,
            'tsub': tsub,
            'merge_thr': merge_thr, 
            'min_SNR': min_SNR,
            'rval_thr': rval_thr,
            'use_cnn': True,
            'min_cnn_thr': cnn_thr,
            'cnn_lowest': cnn_lowest}

opts = params.CNMFParams(params_dict=opts_dict)

## Setup a cluster
# start a cluster for parallel processing (if a cluster already exists it will be closed and a new session will be opened)
if 'dview' in locals():
    cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)

# first we create a motion correction object with the parameters specified
mc = MotionCorrect(fnames, dview=dview, **opts.get_group('motion'))
# note that the file is not loaded in memory


## Motion Correction
#Run piecewise-rigid motion correction using NoRMCorre
mc.motion_correct(save_movie=True)
m_els = cm.load(mc.fname_tot_els)
border_to_0 = 0 if mc.border_nan is 'copy' else mc.border_to_0 
    # maximum shift to be used for trimming against NaNs
    
    
    
## MEMORY MAPPING
# memory map the file in order 'C'
fname_new = cm.save_memmap(mc.mmap_file, base_name='memmap_', order='C',
                           border_to_0=border_to_0) # exclude borders

# now load the file
Yr, dims, T = cm.load_memmap(fname_new)
images = np.reshape(Yr.T, [T] + list(dims), order='F') 
    #load frames in python format (T x X x Y)
    
## restart cluster to clean up memory
cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)


## RUN CNMF ON PATCHES

# First extract spatial and temporal components on patches and combine them
# for this step deconvolution is turned off (p=0)
opts.change_params({'p': 0})
cnm = cnmf.CNMF(n_processes, params=opts, dview=dview)
cnm = cnm.fit(images)

## plot contours of found components
Cn = cm.local_correlations(images.transpose(1,2,0))
Cn[np.isnan(Cn)] = 0
cnm.estimates.plot_contours_nb(img=Cn)

## RE-RUN seeded CNMF on accepted patches to refine and perform deconvolution 
cnm.params.change_params({'p': p})
cnm2 = cnm.refit(images, dview=dview)


## COMPONENT EVALUATION
# the components are evaluated in three ways:
#   a) the shape of each component must be correlated with the data
#   b) a minimum peak SNR is required over the length of a transient
#   c) each shape passes a CNN based classifier

cnm2.estimates.evaluate_components(images, cnm2.params, dview=dview)

## PLOT COMPONENTS
cnm2.estimates.plot_contours_nb(img=Cn, idx=cnm2.estimates.idx_components)

# accepted components
cnm2.estimates.nb_view_components(img=Cn, idx=cnm2.estimates.idx_components)

## Extract DF/F values
cnm2.estimates.detrend_df_f(quantileMin=8, frames_window=250)

# Select only high quality components
cnm2.estimates.select_components(use_object=True)

cnm2.estimates.nb_view_components(img=Cn, denoised_color='red')
print('you may need to change the data rate to generate this one: use jupyter notebook --NotebookApp.iopub_data_rate_limit=1.0e10 before opening jupyter notebook')


## Movie
# View movie with the results
#cnm2.estimates.play_movie(images, q_max=99.9, gain_res=2, magnification=2,bpx=border_to_0,nclude_bck=False)
# reconstruct denoised movie
#denoised = cm.movie(cnm2.estimates.A.dot(cnm2.estimates.C) + \cnm2.estimates.b.dot(cnm2.estimates.f)).reshape(dims + (-1,), order='F').transpose([2, 0, 1])



## save (hdf5) & stop server
save_results = False
if save_results:
    cnm2.save('analysis_results.hdf5')
    
cm.stop_server(dview=dview)



