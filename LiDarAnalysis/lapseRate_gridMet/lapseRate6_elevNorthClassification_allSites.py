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

#%% load data nrc 2010
ascii_grid_veg_indexnrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_nrc.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_nrc = [np.min(ascii_grid_veg_indexnrc[:,0]),np.max(ascii_grid_veg_indexnrc[:,0]),np.min(ascii_grid_veg_indexnrc[:,1]),np.max(ascii_grid_veg_indexnrc[:,1])]
lat_long_degre_nrc = [39.978183, 40.025319, -105.58556840332032, -105.52620245393197]

# slope less than 30
ascii_grid_veg_index_lsnrc = ascii_grid_veg_indexnrc[(ascii_grid_veg_indexnrc[:,5]>-90)&(ascii_grid_veg_indexnrc[:,4]<30)]

#y = -0.007512x + 14.168475  
  # R² = 0.969091
elev_nrc = ascii_grid_veg_index_lsnrc[:,2]
tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475 

veg_snow_temp_nrc = np.vstack([ascii_grid_veg_index_lsnrc[:,0],ascii_grid_veg_index_lsnrc[:,1].astype(int),ascii_grid_veg_index_lsnrc[:,2],
                               ascii_grid_veg_index_lsnrc[:,3],ascii_grid_veg_index_lsnrc[:,4],ascii_grid_veg_index_lsnrc[:,5],
                               ascii_grid_veg_index_lsnrc[:,6],tempDJF_gm_nrc,tempDJF_gm_nrc]).T

veg_snow_temp_dfnrc = pd.DataFrame(veg_snow_temp_nrc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
elevationBand_nrc = np.arange(min(veg_snow_temp_dfnrc['z']),max(veg_snow_temp_dfnrc['z']),65)
veg_snow_temp_clsfnrc = classificationElevation(veg_snow_temp_nrc,elevationBand_nrc,elev_class_inx)

meanTemp_elevClass_nrc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfnrc,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfnrc,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_nrc, fSCA_0u_shl_nrc = fsca0pn_fscaUnTr(fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc, elev_class_inx)

#%% load data Jemez 2010
ascii_grid_veg_indexJmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_jmz.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_Jmz = [np.min(ascii_grid_veg_indexJmz[:,0]),np.max(ascii_grid_veg_indexJmz[:,0]),np.min(ascii_grid_veg_indexJmz[:,1]),np.max(ascii_grid_veg_indexJmz[:,1])]
lat_long_degre_Jmz = [35.863504, 35.918451, -106.60600014573875, -106.54061446340079]

# slope less than 30
ascii_grid_veg_index_lsJmz = ascii_grid_veg_indexJmz[(ascii_grid_veg_indexJmz[:,5]>-90)&(ascii_grid_veg_indexJmz[:,4]<30)]

elev_jmz = ascii_grid_veg_index_lsJmz[:,2]
tempDJF_gm_jmz = -0.0046*elev_jmz + 7.7058 

veg_snow_temp_Jmz = np.vstack([ascii_grid_veg_index_lsJmz[:,0],ascii_grid_veg_index_lsJmz[:,1].astype(int),ascii_grid_veg_index_lsJmz[:,2],
                               ascii_grid_veg_index_lsJmz[:,3],ascii_grid_veg_index_lsJmz[:,4],ascii_grid_veg_index_lsJmz[:,5],
                               ascii_grid_veg_index_lsJmz[:,6],tempDJF_gm_jmz,tempDJF_gm_jmz]).T

veg_snow_temp_dfJmz = pd.DataFrame(veg_snow_temp_Jmz, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_jmz = np.arange(min(veg_snow_temp_dfJmz['z']),max(veg_snow_temp_dfJmz['z']),90)
veg_snow_temp_clsfJmz = classificationElevation(veg_snow_temp_Jmz,elevationBand_jmz,elev_class_inx)

meanTemp_elevClass_jmz = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfJmz,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfJmz,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz = fsca0pn_fscaUnTr(fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz, elev_class_inx)

#%% load data sagehen creek 17 April 2016
ascii_grid_veg_index17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc17A.npy')
#calculating topo_dimension from cut northness

# slope less than 30
ascii_grid_veg_index_ls17a = ascii_grid_veg_index17a[(ascii_grid_veg_index17a[:,5]>-90)&(ascii_grid_veg_index17a[:,4]<30)]

elev_sc17a = ascii_grid_veg_index_ls17a[:,2]
tempDJF_sc17a = -0.0016 * elev_sc17a + 1.802

veg_snow_temp_sc17a = np.vstack([ascii_grid_veg_index_ls17a[:,0],ascii_grid_veg_index_ls17a[:,1].astype(int),ascii_grid_veg_index_ls17a[:,2],
                                 ascii_grid_veg_index_ls17a[:,3],ascii_grid_veg_index_ls17a[:,4],ascii_grid_veg_index_ls17a[:,5],
                                 ascii_grid_veg_index_ls17a[:,6],tempDJF_sc17a,tempDJF_sc17a]).T

veg_snow_temp_sc_df17a = pd.DataFrame(veg_snow_temp_sc17a, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc17a[700780:710057,:]

elevationBand_sc = np.arange(min(veg_snow_temp_sc_df17a['z']),max(veg_snow_temp_sc_df17a['z']),67)
veg_snow_temp_clsfsc17a = classificationElevation(veg_snow_temp_sc17a,elevationBand_sc,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc17a,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfsc17a,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a = fsca0pn_fscaUnTr(fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a, elev_class_inx)


#%% load data sagehen creek 26 march 2010
ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(ascii_grid_veg_indexKrew[:,0]),np.max(ascii_grid_veg_indexKrew[:,0]),np.min(ascii_grid_veg_indexKrew[:,1]),np.max(ascii_grid_veg_indexKrew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

# slope less than 30
ascii_grid_veg_index_lsKrew = ascii_grid_veg_indexKrew[(ascii_grid_veg_indexKrew[:,5]>-90)&(ascii_grid_veg_indexKrew[:,4]<30)]

elev_krew = ascii_grid_veg_index_lsKrew[:,2]
tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339
 #R² = 0.5074

veg_snow_temp_Krew = np.vstack([ascii_grid_veg_index_lsKrew[:,0],ascii_grid_veg_index_lsKrew[:,1].astype(int),ascii_grid_veg_index_lsKrew[:,2],
                                ascii_grid_veg_index_lsKrew[:,3],ascii_grid_veg_index_lsKrew[:,4],ascii_grid_veg_index_lsKrew[:,5],
                                ascii_grid_veg_index_lsKrew[:,6],tempDJF_gm_krew,tempDJF_gm_krew]).T

veg_snow_temp_dfKrew = pd.DataFrame(veg_snow_temp_Krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_krew = np.arange(min(veg_snow_temp_dfKrew['z']),max(veg_snow_temp_dfKrew['z']),67)
veg_snow_temp_clsfKrew = classificationElevation(veg_snow_temp_Krew,elevationBand_krew,elev_class_inx)
aaaa_test2 = veg_snow_temp_clsfKrew[700780:710057,:]

meanTemp_elevClass_krew = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfKrew,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfKrew,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_Krew, fSCA_0u_shl_Krew = fsca0pn_fscaUnTr(fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew, elev_class_inx)

#%% ploting DJF temp
fSCA_0u = [fSCA_0u_exp_nrc, fSCA_0u_shl_nrc, fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz, 
           fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a, fSCA_0u_exp_Krew, fSCA_0u_shl_Krew]
meanTemp = [meanTemp_elevClass_nrc,meanTemp_elevClass_nrc,meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,
            meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew]

label = ['NRC 2010 -- exposed','NRC 2010 -- sheltered','JMZ 2010 -- exposed','JMZ 2010 -- sheltered',
         'SCW 17April2016 -- exposed','SCW 17April2016 -- sheltered','Krew 2010 -- exposed','Krew 2010 -- sheltered']
color = ['lightblue','navy','darkgreen','lightgreen','gold','darkorange','orange','red']

location = ['lower left','lower left','uper left','lower left']
plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (4):
    plt.subplot(221+i)
    plt.scatter(meanTemp[2*i],fSCA_0u[2*i], label = label[2*i], s=20**2, color = color[2*i])
    plt.scatter(meanTemp[2*i+1],fSCA_0u[2*i+1], label = label[2*i+1], s=20**2, color = color[2*i+1])

    plt.ylabel('delta fSCA% / fSCA_open', fontsize=25)
    plt.yticks(fontsize=15)
    plt.xlabel('average temp of DJF (C)', fontsize=25)
    plt.xticks(fontsize=15)
    plt.legend(fontsize=20, loc = location[i])

plt.title('(delta fSCA)/fSCA_op in different temp lapse rate and northness in 4 sites', fontsize=35, y=2.25, x=-0.15) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrth_all2.png')











