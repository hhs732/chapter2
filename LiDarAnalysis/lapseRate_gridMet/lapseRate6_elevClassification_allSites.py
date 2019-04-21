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
               

def fSCAcalculation_elevClassific(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification

    fSCA_unTr = []
    fSCA_0pen = []
    
    for indx in range (len(elev_class_inx)): 
        #under tree with snow, exposed
        ascii_ut_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])])
        ascii_ut_nSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])])
        #fSCA under tree
        fSCA_ut = 100.*ascii_ut_snow/(ascii_ut_snow+ascii_ut_nSnow)
        fSCA_unTr.append(fSCA_ut)
    
        #0pen with snow
        ascii_0p_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])])
        #0pen no snow
        ascii_0p_nSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])])
        #fSCA 0pen
        fSCA_0p = 100.*ascii_0p_snow/(ascii_0p_nSnow+ascii_0p_snow)
        fSCA_0pen.append(fSCA_0p)
    
    return fSCA_unTr, fSCA_0pen

def fsca0pn_fscaUnTr(fSCA_unTr, fSCA_0pen,elev_class_inx):
    
    fSCA_0u = []
    for indx in range (len(elev_class_inx)): 
        fSCA_0u_each = (fSCA_0pen[indx] - fSCA_unTr[indx])/fSCA_0pen[indx]
        fSCA_0u.append(fSCA_0u_each)

    return fSCA_0u

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

#%% load data nrc 2010
ascii_grid_veg_indexnrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_nrc.npy')
# slope less than 30
ascii_grid_veg_index_lsnrc = ascii_grid_veg_indexnrc[(ascii_grid_veg_indexnrc[:,5]>-90)&(ascii_grid_veg_indexnrc[:,4]<30)]
  
elev_nrc = ascii_grid_veg_index_lsnrc[:,2]
tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475 #y = -0.007512x + 14.168475  # R² = 0.969091

veg_snow_temp_nrc = np.vstack([ascii_grid_veg_index_lsnrc[:,0],ascii_grid_veg_index_lsnrc[:,1].astype(int),ascii_grid_veg_index_lsnrc[:,2],
                               ascii_grid_veg_index_lsnrc[:,3],ascii_grid_veg_index_lsnrc[:,4],ascii_grid_veg_index_lsnrc[:,5],
                               ascii_grid_veg_index_lsnrc[:,6],tempDJF_gm_nrc,tempDJF_gm_nrc]).T

veg_snow_temp_dfnrc = pd.DataFrame(veg_snow_temp_nrc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
elevationBand_nrc = np.arange(min(veg_snow_temp_dfnrc['z']),max(veg_snow_temp_dfnrc['z']),65)
veg_snow_temp_clsfnrc = classificationElevation(veg_snow_temp_nrc,elevationBand_nrc,elev_class_inx)

meanTemp_elevClass_nrc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfnrc,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_nrc, fSCA_0p_nrc = fSCAcalculation_elevClassific(veg_snow_temp_clsfnrc,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_nrc = fsca0pn_fscaUnTr(fSCA_ut_nrc, fSCA_0p_nrc, elev_class_inx)

#%% load data Jemez 2010
ascii_grid_veg_indexJmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_jmz.npy')
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
fSCA_ut_Jmz, fSCA_0p_Jmz = fSCAcalculation_elevClassific(veg_snow_temp_clsfJmz,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_Jmz = fsca0pn_fscaUnTr(fSCA_ut_Jmz, fSCA_0p_Jmz, elev_class_inx)

#%% load data sagehen creek 26 march 2016
ascii_grid_veg_index26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc26m.npy')
# slope less than 30
ascii_grid_veg_index_ls26m = ascii_grid_veg_index26m[(ascii_grid_veg_index26m[:,5]>-90)&(ascii_grid_veg_index26m[:,4]<30)]

# temp-elev lapse rate
elev_sc_26m = ascii_grid_veg_index_ls26m[:,2] 
tempDJF_sc_26m = -0.0016 * elev_sc_26m + 1.802 

veg_snow_temp_sc26m = np.vstack([ascii_grid_veg_index_ls26m[:,0],ascii_grid_veg_index_ls26m[:,1].astype(int),ascii_grid_veg_index_ls26m[:,2],
                                 ascii_grid_veg_index_ls26m[:,3],ascii_grid_veg_index_ls26m[:,4],ascii_grid_veg_index_ls26m[:,5],
                                 ascii_grid_veg_index_ls26m[:,6],tempDJF_sc_26m,tempDJF_sc_26m]).T

veg_snow_temp_sc_df26m = pd.DataFrame(veg_snow_temp_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_sc26m = np.arange(min(veg_snow_temp_sc_df26m['z']),max(veg_snow_temp_sc_df26m['z']),67)
veg_snow_temp_clsfsc26m = classificationElevation(veg_snow_temp_sc26m,elevationBand_sc26m,elev_class_inx)

meanTemp_elevClass_sc26m = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc26m,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_sc26m, fSCA_0p_sc26m = fSCAcalculation_elevClassific(veg_snow_temp_clsfsc26m,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_sc26m = fsca0pn_fscaUnTr(fSCA_ut_sc26m, fSCA_0p_sc26m, elev_class_inx)

#%% load data sagehen creek 17 April 2016
ascii_grid_veg_index17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc17A.npy')
# slope less than 30
ascii_grid_veg_index_ls17a = ascii_grid_veg_index17a[(ascii_grid_veg_index17a[:,5]>-90)&(ascii_grid_veg_index17a[:,4]<30)]

elev_sc17a = ascii_grid_veg_index_ls17a[:,2]
tempDJF_sc17a = -0.0016 * elev_sc17a + 1.802

veg_snow_temp_sc17a = np.vstack([ascii_grid_veg_index_ls17a[:,0],ascii_grid_veg_index_ls17a[:,1].astype(int),ascii_grid_veg_index_ls17a[:,2],
                                 ascii_grid_veg_index_ls17a[:,3],ascii_grid_veg_index_ls17a[:,4],ascii_grid_veg_index_ls17a[:,5],
                                 ascii_grid_veg_index_ls17a[:,6],tempDJF_sc17a,tempDJF_sc17a]).T

veg_snow_temp_sc_df17a = pd.DataFrame(veg_snow_temp_sc17a, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_sc17a = np.arange(min(veg_snow_temp_sc_df17a['z']),max(veg_snow_temp_sc_df17a['z']),67)
veg_snow_temp_clsfsc17a = classificationElevation(veg_snow_temp_sc17a,elevationBand_sc17a,elev_class_inx)

meanTemp_elevClass_sc17a = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc17a,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_sc17a, fSCA_0p_sc17a = fSCAcalculation_elevClassific(veg_snow_temp_clsfsc17a,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_sc17a = fsca0pn_fscaUnTr(fSCA_ut_sc17a, fSCA_0p_sc17a, elev_class_inx)

#%% load data sagehen creek 24 May 2016
ascii_grid_veg_index24m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc24May.npy')
# slope less than 30
ascii_grid_veg_index_ls24m = ascii_grid_veg_index24m[(ascii_grid_veg_index24m[:,5]>-90)&(ascii_grid_veg_index24m[:,4]<30)]

elev_sc24m = ascii_grid_veg_index_ls24m[:,2]
tempDJF_sc24m = -0.0016 * elev_sc24m + 1.802

veg_snow_temp_sc24m = np.vstack([ascii_grid_veg_index_ls24m[:,0],ascii_grid_veg_index_ls24m[:,1].astype(int),ascii_grid_veg_index_ls24m[:,2],
                                 ascii_grid_veg_index_ls24m[:,3],ascii_grid_veg_index_ls24m[:,4],ascii_grid_veg_index_ls24m[:,5],
                                 ascii_grid_veg_index_ls24m[:,6],tempDJF_sc24m,tempDJF_sc24m]).T

veg_snow_temp_sc_df24m = pd.DataFrame(veg_snow_temp_sc24m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_sc_24m = np.arange(min(veg_snow_temp_sc_df17a['z']),max(veg_snow_temp_sc_df17a['z']),67)
veg_snow_temp_clsfsc24m = classificationElevation(veg_snow_temp_sc24m,elevationBand_sc_24m,elev_class_inx)

meanTemp_elevClass_sc24m = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc24m,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_sc24m, fSCA_0p_sc24m = fSCAcalculation_elevClassific(veg_snow_temp_clsfsc24m,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_sc24m = fsca0pn_fscaUnTr(fSCA_ut_sc24m, fSCA_0p_sc24m, elev_class_inx)


#%% load data sagehen creek 26 march 2010
ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
# slope less than 30
ascii_grid_veg_index_lsKrew = ascii_grid_veg_indexKrew[(ascii_grid_veg_indexKrew[:,5]>-90)&(ascii_grid_veg_indexKrew[:,4]<30)]

elev_krew = ascii_grid_veg_index_lsKrew[:,2]
tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339  #R² = 0.5074
 
veg_snow_temp_Krew = np.vstack([ascii_grid_veg_index_lsKrew[:,0],ascii_grid_veg_index_lsKrew[:,1].astype(int),ascii_grid_veg_index_lsKrew[:,2],
                                ascii_grid_veg_index_lsKrew[:,3],ascii_grid_veg_index_lsKrew[:,4],ascii_grid_veg_index_lsKrew[:,5],
                                ascii_grid_veg_index_lsKrew[:,6],tempDJF_gm_krew,tempDJF_gm_krew]).T

veg_snow_temp_dfKrew = pd.DataFrame(veg_snow_temp_Krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_krew = np.arange(min(veg_snow_temp_dfKrew['z']),max(veg_snow_temp_dfKrew['z']),70)
veg_snow_temp_clsfKrew = classificationElevation(veg_snow_temp_Krew,elevationBand_krew,elev_class_inx)
aaaa_test2 = veg_snow_temp_clsfKrew[700780:710057,:]

meanTemp_elevClass_krew = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfKrew,elev_class_inx)   
# under canopy and 0pen fsca
fSCA_ut_Krew, fSCA_0p_Krew = fSCAcalculation_elevClassific(veg_snow_temp_clsfKrew,elev_class_inx)

#fsca open -minus- fsca under canopy
fSCA_0u_Krew = fsca0pn_fscaUnTr(fSCA_ut_Krew, fSCA_0p_Krew, elev_class_inx)

#%% ploting DJF temp

tempMean = [meanTemp_elevClass_sc26m,meanTemp_elevClass_sc17a,meanTemp_elevClass_sc24m,
            meanTemp_elevClass_krew,meanTemp_elevClass_nrc,meanTemp_elevClass_jmz]

fsca_0u = [fSCA_0u_sc26m,fSCA_0u_sc17a,fSCA_0u_sc24m,
           fSCA_0u_Krew,fSCA_0u_nrc,fSCA_0u_Jmz]
label = ['Sagehen Creek 26March2016','Sagehen Creek 17April2016','Sagehen Creek 18May2016',
         'Krew 2010','NRC 2010','JMZ 2010']
color = ['orange','red','brown','darkcyan','green','navy']

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))
for i in range (6):

    plt.scatter(tempMean[i],fsca_0u[i], label = label[i], s=10**2.2, color = color[i])
    plt.plot(tempMean[i],fsca_0u[i], color = color[i])

plt.ylabel('(delta fSCA)/fSCA_op', fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (C)', fontsize=30)
plt.xticks(fontsize=30)
plt.legend(fontsize=25, loc = 'lower left')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_all2.png')











