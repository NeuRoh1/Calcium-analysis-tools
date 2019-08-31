#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 10:29:43 2019

@author: Chang-Eop

caiman_helper.py 파일이 같은 디렉토리, 혹은 pathway에 추가되어있어야 함
2p_analysis.py를 실행시킨 이후에만 실행 가능
ROI의 index: 1부터 순서대로. python의 0부터 시작하는 index를 따르지 않음!
ROI index 확인하는 방법: 
    #footprints of individual components 코드를 실행 후 개별 ROI의 footprint에 커서 위치
    #footprints of summed components에서 커서 위치시 겹치는 영역들 index도 합산되므로 헷갈리지 말것
    
supurious ROI 제거 방법:
    #exclude supurious ROIs 코드 활성화 시키고 실행. 제거할 ROI의 index를 [1,5,6]과 같은 형태로 입력
    F,C,A 모두에 대해(!) 제거코드 실행한 이후, line 42 부터 재실행
"""

#results analysis using estimates
import matplotlib.colors as color
from random import random
import pandas as pd
import caiman_helper_Tdr as helper 
import scipy as sp

#ROI margin detection parameter
neighbor_type = 1 


colors = [(1,1,1)] + [(random(),random(),random()) for i in range(255)]
new_map = color.LinearSegmentedColormap.from_list('new_map', colors, N=256)


F = cnm2.estimates.F_dff #F: DF/F
C = cnm2.estimates.C #C: denoised 
S = cnm2.estimates.S #S: deconvolved spikes
A = cnm2.estimates.A
A = np.array(sp.sparse.csr_matrix.todense(A))



#If you exclude supurious ROIs, execute codes from here again!
Rois = np.zeros((C.shape[0],dims[0],dims[1]))
for i in range(Rois.shape[0]):
    mat = np.reshape(A[:,i], dims, order='F')
    mat[mat>0]=1
    Rois[i,]=helper.margin_detect(mat, neighbor_type)*(i+1)
    
#footprints of summed components
Rois_sum = np.sum(Rois, axis = 0)
plt.figure(); plt.imshow(Rois_sum, cmap = new_map)



#footprints of individual components (각 ROI의 index는 이걸로 확인할것)
for i in range(Rois.shape[0]):
    plt.figure(); plt.imshow(Rois[i,], cmap = new_map)


#construct and visualize correlation matrix
F_corr = np.corrcoef(F)
np.fill_diagonal(F_corr, 0)
plt.figure; plt.imshow(F_corr)



#


#To view the i_th spaital components
i = 5
plt.figure(); plt.imshow(np.reshape(A[:,i-1], dims, order='F')) #cmap = new_map)

#To plot the i_th componet's trace
plt.figure(); plt.plot(C[i-1])
    


##exclude supurious ROIs 
## When want excluding roi #1, #3 just use [1,3]
## must execute F, C, A at once 
#---execute below------
Exc=[1,2,3] # insert ROI numbers to exclude
F = helper.exclude_rois(F, Exc)
C = helper.exclude_rois(C, Exc)
S = helper.exclude_rois(S, Exc)
A = helper.exclude_rois_map(A, Exc)




#save csv files
F_df = pd.DataFrame(F) #F: DF/F
F_df.to_csv('P4_Rg2_100um_F.csv')

C_df = pd.DataFrame(C) #C: denoised 
C_df.to_csv('P4_Rg2_100um_C.csv')

S = pd.DataFrame(S) #S: spikes
S.to_csv('P4_Rg2_100um_S.csv')


