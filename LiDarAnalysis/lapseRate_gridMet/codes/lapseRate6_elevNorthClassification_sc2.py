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

def number0fGridsInEachElevationBand(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification
    
    numGrids = []
    for indx in range (len(elev_class_inx)):
        numGrid = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])])
        numGrids.append(numGrid)
    
    return numGrids

#%% load data sagehen creek 26 march 2016
ascii_grid_veg_index26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc26m.npy')

lat_long_nad83_sc = [np.min(ascii_grid_veg_index26m[:,0]),np.max(ascii_grid_veg_index26m[:,0]),np.min(ascii_grid_veg_index26m[:,1]),np.max(ascii_grid_veg_index26m[:,1])]
lat_long_degre_sc = [39.403596, 39.464694, -120.31719287846941, -120.23349053231303]

#ascii_grid_veg_index26m = classificationGridCells(ascii_index_pr26m,gridCells_m)
# slope less than 30
ascii_grid_veg_index_ls26m = ascii_grid_veg_index26m[(ascii_grid_veg_index26m[:,5]>-90)&(ascii_grid_veg_index26m[:,4]<30)]

# temp-elev lapse rate
#y = -0.0014x + 2.0155   R² = 0.8725
elev_sc = ascii_grid_veg_index_ls26m[:,2]

tempDJF_sc = -0.0016 * elev_sc + 1.802 

veg_snow_temp_sc26m = np.vstack([ascii_grid_veg_index_ls26m[:,0],ascii_grid_veg_index_ls26m[:,1].astype(int),ascii_grid_veg_index_ls26m[:,2],
                                 ascii_grid_veg_index_ls26m[:,3],ascii_grid_veg_index_ls26m[:,4],ascii_grid_veg_index_ls26m[:,5],
                                 ascii_grid_veg_index_ls26m[:,6],tempDJF_sc,tempDJF_sc]).T

veg_snow_temp_sc_df26m = pd.DataFrame(veg_snow_temp_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_sc = np.arange(min(veg_snow_temp_sc_df26m['z']),max(veg_snow_temp_sc_df26m['z']),67)
veg_snow_temp_sc_clsf26m = classificationElevation(veg_snow_temp_sc26m,elevationBand_sc)

numGrids_sc26m = number0fGridsInEachElevationBand(veg_snow_temp_sc_clsf26m)
aaaa_test = veg_snow_temp_sc_clsf26m[700780:710057,:]

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(veg_snow_temp_sc_clsf26m)   
# under canopy and 0pen fsca
fSCA_ut_sc26m, fSCA_0p_sc26m = fSCAcalculation_elevClassific(veg_snow_temp_sc_clsf26m)

#fsca open -minus- fsca under canopy
fSCA_0u_sc26m = fsca0pn_fscaUnTr(fSCA_ut_sc26m,fSCA_0p_sc26m)

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

elevationBand = np.arange(min(veg_snow_temp_sc_df17a['z']),max(veg_snow_temp_sc_df17a['z']),67)
elev_class_inx = [101,102,103,104,105,106,107,108,109,110]

veg_snow_temp_clsfsc17a = classificationElevation(veg_snow_temp_sc17a,elevationBand,elev_class_inx)

numGrids_sc17a = number0fGridsInEachElevationBand(veg_snow_temp_clsfsc17a,elev_class_inx)

meanTemp_elevClass = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc17a,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfsc17a,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a = fsca0pn_fscaUnTr(fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a, elev_class_inx)


#%% load data sagehen creek 24 May 2016
ascii_grid_veg_index24m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc24May.npy')
#calculating topo_dimension from cut northness

# slope less than 30
ascii_grid_veg_index_ls24m = ascii_grid_veg_index24m[(ascii_grid_veg_index24m[:,5]>-90)&(ascii_grid_veg_index24m[:,4]<30)]

elev_sc24m = ascii_grid_veg_index_ls24m[:,2]
 
tempDJF_sc24m = -0.0016 * elev_sc24m + 1.802

veg_snow_temp_sc24m = np.vstack([ascii_grid_veg_index_ls24m[:,0],ascii_grid_veg_index_ls24m[:,1].astype(int),ascii_grid_veg_index_ls24m[:,2],
                                 ascii_grid_veg_index_ls24m[:,3],ascii_grid_veg_index_ls24m[:,4],ascii_grid_veg_index_ls24m[:,5],
                                 ascii_grid_veg_index_ls24m[:,6],tempDJF_sc24m,tempDJF_sc24m]).T

veg_snow_temp_sc_df24m = pd.DataFrame(veg_snow_temp_sc24m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc24m[700780:710057,:]

veg_snow_temp_sc_clsf24m = classificationElevation(veg_snow_temp_sc24m,elevationBand_sc)

meanTemp_elevClass_sc24m = calculationMeanTempforEachElevClassification(veg_snow_temp_sc_clsf24m)   

# under canopy and 0pen fsca
fSCA_ut_sc24m, fSCA_0p_sc24m = fSCAcalculation_elevClassific(veg_snow_temp_sc_clsf24m)

#fsca open -minus- fsca under canopy
fSCA_0u_sc24m = fsca0pn_fscaUnTr(fSCA_ut_sc24m,fSCA_0p_sc24m)

#%% ploting DJF temp

fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_elevClass,fSCA_0u_exp_sc17a, label = 'Sagehen 17April2016 --- exposed', s=20**2, color = 'gold')#r2=%.6f'%(r_value_tempFscaKrew**2)
ax.scatter(meanTemp_elevClass,fSCA_0u_shl_sc17a, label = 'Sagehen 17April2016 --- sheltered', s=20**2, color = 'darkorange')

#ax.scatter(meanTemp_elevClass_sc,fSCA_0u_sc26m, label = 'Sagehen Creek 26March2016', s=20**2, color = 'blue')#r2=%.6f'%(r_value_tempFsca26M**2)
#ax.scatter(meanTemp_elevClass_sc17a,fSCA_0u_sc17a, label = 'Sagehen Creek 17April2016', s=20**2, color = 'green')#r2=%.6f'%(r_value_tempFsca26M**2)
#ax.scatter(meanTemp_elevClass_sc24m,fSCA_0u_sc24m, label = 'Sagehen Creek 24May2016', s=20**2, color = 'orange')#r2=%.6f'%(r_value_tempFsca26M**2)
#
#plt.plot(meanTemp_elevClass_sc,fSCA_0u_sc26m, color = 'blue')
#plt.plot(meanTemp_elevClass_sc17a,fSCA_0u_sc17a, color = 'green')
#plt.plot(meanTemp_elevClass_sc24m,fSCA_0u_sc24m, color = 'orange')

plt.ylabel('delta fSCA/fSCA_op', fontsize=35)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (c)', fontsize=35)
plt.xticks(fontsize=30)

plt.legend(fontsize=20, loc = 'upper left')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrth_sc17A.png')











