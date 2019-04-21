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
ascii_grid_veg_indexnrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_nrc.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_nrc = [np.min(ascii_grid_veg_indexnrc[:,0]),np.max(ascii_grid_veg_indexnrc[:,0]),np.min(ascii_grid_veg_indexnrc[:,1]),np.max(ascii_grid_veg_indexnrc[:,1])]
lat_long_degre_nrc = [39.978183, 40.025319, -105.58556840332032, -105.52620245393197]

# slope less than 30
ascii_grid_veg_index_lsnrc = ascii_grid_veg_indexnrc[(ascii_grid_veg_indexnrc[:,5]>-90)&(ascii_grid_veg_indexnrc[:,4]<30)]
aaaa_test = ascii_grid_veg_index_lsnrc[700780:710057,:]

#y = -0.007512x + 14.168475  
  # R² = 0.969091
elev_nrc = ascii_grid_veg_index_lsnrc[:,2]
tempDJF_gm = -0.007512*elev_nrc + 14.168475 

veg_snow_temp_nrc = np.vstack([ascii_grid_veg_index_lsnrc[:,0],ascii_grid_veg_index_lsnrc[:,1].astype(int),ascii_grid_veg_index_lsnrc[:,2],
                                ascii_grid_veg_index_lsnrc[:,3],ascii_grid_veg_index_lsnrc[:,4],ascii_grid_veg_index_lsnrc[:,5],
                                ascii_grid_veg_index_lsnrc[:,6],tempDJF_gm,tempDJF_gm]).T

veg_snow_temp_dfnrc = pd.DataFrame(veg_snow_temp_nrc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand = np.arange(min(veg_snow_temp_dfnrc['z']),max(veg_snow_temp_dfnrc['z']),65)
elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
veg_snow_temp_clsfnrc = classificationElevation(veg_snow_temp_nrc,elevationBand,elev_class_inx)

numGrids_nrc = number0fGridsInEachElevationBand(veg_snow_temp_clsfnrc,elev_class_inx)

meanTemp_elevClass = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfnrc,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfnrc,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_nrc, fSCA_0u_shl_nrc = fsca0pn_fscaUnTr(fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc, elev_class_inx)

#%% ploting DJF temp
fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_elevClass,fSCA_0u_exp_nrc, label = 'NRC 2010 --- exposed', s=20**2, color = 'lightblue')#r2=%.6f'%(r_value_tempFscaKrew**2)
ax.scatter(meanTemp_elevClass,fSCA_0u_shl_nrc, label = 'NRC 2010 --- sheltered', s=20**2, color = 'navy')


plt.ylabel('delta fSCA/fSCA_op', fontsize=35)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (c)', fontsize=35)
plt.xticks(fontsize=30)

plt.legend(fontsize=20, loc = 'lower left')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrt_nrc2.png')











