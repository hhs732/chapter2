import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import itertools
import csv
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import os
import rasterio
from scipy import stats
#%% functions
# classification based on lat and long in 100x100 pixels
#veg_snow_temp_sc_df_tst = pd.DataFrame(veg_snow_temp_sc[0:1000,:], columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
#lat_long_df_sc_tst = lat_long_df_sc[['x','y']][0:200]
#cluster100_sc = []
#for row in veg_snow_temp_sc_df.itertuples():
#    indexFinl = np.where(((lat_long_df_sc['x']-row.x<=100) & (lat_long_df_sc['y']-row.y<=100) & (lat_long_df_sc['x']-row.x>=0) & (lat_long_df_sc['y']-row.y>=0)))[0]#
#    if (np.size(indexFinl))>0:
#        cluster100_sc.append([lat_long_df_sc['x'][indexFinl[0]],lat_long_df_sc['y'][indexFinl[0]]])

def classificationElevation(points_df_ls,elevationBand,elev_class_inx):#0.25
    
    for indx in range (len(elev_class_inx)):
        points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[indx])&(points_df_ls[:,2]<elevationBand[indx+1])]=elev_class_inx[indx]
    
    points_df_ls[:,8][(points_df_ls[:,8]<=100)]=110    
    
    return points_df_ls#, pointFile_df
               

def fSCAcalculation_elevNorthClassific(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification

    fSCA_unTr_exp = []
    fSCA_unTr_shl = []
    fSCA_0pen_exp = []
    fSCA_0pen_shl = []
    
    for indx in range (len(elev_class_inx)): 
        #under tree with snow, exposed
        ascii_ut_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #under tree no snow, exposed
        ascii_ut_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #under tree with snow, sheltered
        ascii_ut_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
        #under tree no snow, sheltered
        ascii_ut_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])

        #fSCA under tree, exposed
        fSCA_ut_exp = 100.*ascii_ut_snow_exp/(ascii_ut_snow_exp+ascii_ut_noSnow_exp)
        #fSCA under tree, sheltered
        fSCA_ut_shl = 100.*ascii_ut_snow_shl/(ascii_ut_snow_shl+ascii_ut_nSnow_shl)
    
        fSCA_unTr_exp.append(fSCA_ut_exp)
        fSCA_unTr_shl.append(fSCA_ut_shl)
    
        #0pen with snow, exposed
        ascii_0p_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #0pen no snow, exposed
        ascii_0p_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #0pen with snow, sheltered
        ascii_0p_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
        #0pen no snow, sheltered
        ascii_0p_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])

        #fSCA 0pen, exposed
        fSCA_0p_exp = 100.*ascii_0p_snow_exp/(ascii_0p_snow_exp+ascii_0p_noSnow_exp)
        #fSCA 0pen, sheltered
        fSCA_0p_shl = 100.*ascii_0p_snow_shl/(ascii_0p_snow_shl+ascii_0p_nSnow_shl)
    
        fSCA_0pen_exp.append(fSCA_0p_exp)
        fSCA_0pen_shl.append(fSCA_0p_shl)
    
    return fSCA_unTr_exp, fSCA_unTr_shl, fSCA_0pen_exp, fSCA_0pen_shl

def fsca0pn_fscaUnTr(fSCA_ut_exp,fSCA_ut_shl,fSCA_0p_exp,fSCA_0p_shl,elev_class_inx):
    
    fSCA_0u_exp = []
    fSCA_0u_shl = []
    for indx in range (len(elev_class_inx)): 
        #fSCA in exposed area
        fSCA_0u_exp_each = (fSCA_0p_exp[indx] - fSCA_ut_exp[indx])/fSCA_0p_exp[indx]
        fSCA_0u_exp.append(fSCA_0u_exp_each)
        
        #fSCA in shelter area
        fSCA_0u_shl_each = (fSCA_0p_shl[indx] - fSCA_ut_shl[indx])/fSCA_0p_shl[indx]
        fSCA_0u_shl.append(fSCA_0u_shl_each)
        
    return fSCA_0u_exp, fSCA_0u_shl

def calculationMeanTempforEachElevClassification(veg_snow_temp_elev_clsf_file,elev_class_inx):
    
    meanTemp = []
    for indx in range (len(elev_class_inx)): 
        temp_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])][:,7])
        meanTemp.append(temp_mean)    
        
    return meanTemp

def number0fGridsInEachElevationBand(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification
    
    numGrids = []
    for indx in range (len(elev_class_inx)):
        numGrid = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])])
        numGrids.append(numGrid)
    
    return numGrids
#%% load data sagehen creek 26 march 2010
ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(ascii_grid_veg_indexKrew[:,0]),np.max(ascii_grid_veg_indexKrew[:,0]),np.min(ascii_grid_veg_indexKrew[:,1]),np.max(ascii_grid_veg_indexKrew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

# slope less than 30
ascii_grid_veg_index_lsKrew = ascii_grid_veg_indexKrew[(ascii_grid_veg_indexKrew[:,5]>-90)&(ascii_grid_veg_indexKrew[:,4]<30)]

elev_krew = ascii_grid_veg_index_lsKrew[:,2]
tempDJF_gm = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339
 #R² = 0.5074

veg_snow_temp_Krew = np.vstack([ascii_grid_veg_index_lsKrew[:,0],ascii_grid_veg_index_lsKrew[:,1].astype(int),ascii_grid_veg_index_lsKrew[:,2],
                                ascii_grid_veg_index_lsKrew[:,3],ascii_grid_veg_index_lsKrew[:,4],ascii_grid_veg_index_lsKrew[:,5],
                                ascii_grid_veg_index_lsKrew[:,6],tempDJF_gm,tempDJF_gm]).T

veg_snow_temp_dfKrew = pd.DataFrame(veg_snow_temp_Krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand = np.arange(min(veg_snow_temp_dfKrew['z']),max(veg_snow_temp_dfKrew['z']),67)
elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
veg_snow_temp_clsfKrew = classificationElevation(veg_snow_temp_Krew,elevationBand,elev_class_inx)

aaaa_test2 = veg_snow_temp_clsfKrew[700780:710057,:]

numGrids_krew = number0fGridsInEachElevationBand(veg_snow_temp_clsfKrew,elev_class_inx)

meanTemp_elevClass = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfKrew,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfKrew,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_Krew, fSCA_0u_shl_Krew = fsca0pn_fscaUnTr(fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew, elev_class_inx)

#%% ploting DJF temp
fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_elevClass,fSCA_0u_exp_Krew, label = 'Krew 2010 --- exposed', s=20**2, color = 'orange')#r2=%.6f'%(r_value_tempFscaKrew**2)
ax.scatter(meanTemp_elevClass,fSCA_0u_shl_Krew, label = 'Krew 2010 --- sheltered', s=20**2, color = 'red')


plt.ylabel('delta fSCA/fSCA_op', fontsize=35)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (c)', fontsize=35)
plt.xticks(fontsize=30)

plt.legend(fontsize=20, loc = 'upper right')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrt_krew2.png')











