import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import scipy
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

def fSCAcalculation_elevNorthRandom(veg_snow_temp_elev_clsf_file,elev_class_inx,sample_size): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification

    randomNum = np.random.randint(0,sample_size,size=len(veg_snow_temp_elev_clsf_file[:,0]))
    veg_snow_temp_elev_clsf_file[:,6]=randomNum
    
    fSCA_0u_exp = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2} 
    fSCA_0u_shl = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
    
    for indx in range (len(elev_class_inx)): 
        
        fSCA_0u_exp_rand = np.zeros(sample_size) 
        fSCA_0u_shl_rand = np.zeros(sample_size)
        
        for rand in range (sample_size):
            
            #under tree with snow, exposed
            ascii_ut_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow, exposed
            ascii_ut_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree with snow, sheltered
            ascii_ut_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow, sheltered
            ascii_ut_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            
            #fSCA under tree, exposed
            #fSCA_ut_exp = 100.*ascii_ut_snow_exp/(ascii_ut_snow_exp+ascii_ut_noSnow_exp)
            #fSCA under tree, sheltered
            #fSCA_ut_shl = 100.*ascii_ut_snow_shl/(ascii_ut_snow_shl+ascii_ut_nSnow_shl)
            
            #0pen with snow, exposed
            ascii_0p_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen no snow, exposed
            ascii_0p_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen with snow, sheltered
            ascii_0p_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen no snow, sheltered
            ascii_0p_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])

            #fSCA 0pen, exposed
            #fSCA_0p_exp = 100.*ascii_0p_snow_exp/(ascii_0p_snow_exp+ascii_0p_noSnow_exp)
            #fSCA 0pen, sheltered
            #fSCA_0p_shl = 100.*ascii_0p_snow_shl/(ascii_0p_snow_shl+ascii_0p_nSnow_shl)
            
#            fSCA_0u_shl_each = (fSCA_0p_shl-fSCA_ut_shl)/fSCA_0p_shl
#            fSCA_0u_exp_each = (fSCA_0p_exp-fSCA_ut_exp)/fSCA_0p_exp
            
            #fSCA 0pen, exposed
            if ascii_0p_snow_exp+ascii_0p_noSnow_exp == 0:
                fSCA_0p_exp = 0
            else:            
                fSCA_0p_exp = 100.*ascii_0p_snow_exp/(ascii_0p_snow_exp+ascii_0p_noSnow_exp)
            #fSCA 0pen, sheltered
            if ascii_0p_snow_shl+ascii_0p_nSnow_shl == 0:
                fSCA_0p_shl = 0
            else:
                fSCA_0p_shl = 100.*ascii_0p_snow_shl/(ascii_0p_snow_shl+ascii_0p_nSnow_shl)
            
            #fSCA under canopy, exposed
            if ascii_ut_snow_exp+ascii_ut_noSnow_exp == 0:
                fSCA_ut_exp = 0
            else:            
                fSCA_ut_exp = 100.*ascii_ut_snow_exp/(ascii_ut_snow_exp+ascii_ut_noSnow_exp)
            #fSCA under canopy, sheltered
            if ascii_ut_snow_shl+ascii_ut_nSnow_shl == 0:
                fSCA_ut_shl = 0
            else:
                fSCA_ut_shl = 100.*ascii_ut_snow_shl/(ascii_ut_snow_shl+ascii_ut_nSnow_shl)

            if fSCA_0p_exp == 0:
                fSCA_0u_exp_each = -1000
            else:
                fSCA_0u_exp_each = (fSCA_0p_exp-fSCA_ut_exp)/fSCA_0p_exp
            
            if fSCA_0p_shl == 0:
                fSCA_0u_shl_each = -1000
            else:
                fSCA_0u_shl_each = (fSCA_0p_shl-fSCA_ut_shl)/fSCA_0p_shl
                
            fSCA_0u_exp_rand[rand] = fSCA_0u_exp_each
            fSCA_0u_shl_rand[rand] = fSCA_0u_shl_each
    
        fSCA_0u_exp[indx] = fSCA_0u_exp_rand
        fSCA_0u_shl[indx] = fSCA_0u_shl_rand
        
    return fSCA_0u_exp, fSCA_0u_shl #fSCA_unTr_exp, fSCA_unTr_shl, fSCA_0pen_exp, fSCA_0pen_shl

def fSCAcalculation_elevRandom(veg_snow_temp_elev_clsf_file,elev_class_inx,sample_size): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification

    randomNum = np.random.randint(0,sample_size,size=len(veg_snow_temp_elev_clsf_file[:,0]))
    veg_snow_temp_elev_clsf_file[:,6]=randomNum
    
    fSCA_0p = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2} 
    fSCA_ut = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
    
    for indx in range (len(elev_class_inx)): 
        
        fSCA_0p_rand = np.zeros(sample_size) 
        fSCA_ut_rand = np.zeros(sample_size)
        
        for rand in range (sample_size):
            
            #under tree with snow
            ascii_ut_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow
            ascii_ut_noSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #fSCA under tree
            fSCA_ut_each = 100.*ascii_ut_snow/(ascii_ut_snow+ascii_ut_noSnow)
            
            #0pen with snow
            ascii_0p_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen no snow
            ascii_0p_noSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #fSCA 0pen
            fSCA_0p_each = 100.*ascii_0p_snow/(ascii_0p_snow+ascii_0p_noSnow)
            

            fSCA_ut_rand[rand] = fSCA_ut_each
            fSCA_0p_rand[rand] = fSCA_0p_each
    
        fSCA_0p[indx] = fSCA_ut_rand
        fSCA_ut[indx] = fSCA_0p_rand
        
    return fSCA_ut, fSCA_0p


def wilcoxonTest(fSCA_0u_exp_samp,fSCA_0u_shl_samp):
    pvalue = []
    statistic = []
    average_sample_shl = []
    average_sample_exp = []
    significant_sign = []
    
    for iii in range (10):
        
        statistic_each, pvalue_each = scipy.stats.wilcoxon(fSCA_0u_exp_samp[iii], fSCA_0u_shl_samp[iii], zero_method='wilcox', correction=False)
        statistic.append(statistic_each)
        pvalue.append(pvalue_each)
        
        if pvalue_each < 0.05:
            significant = '^'
        else:
            significant = 'o'
        
        significant_sign.append(significant)
        
        average_sample_shl.append(np.average(fSCA_0u_shl_samp[iii]))
        average_sample_exp.append(np.average(fSCA_0u_exp_samp[iii]))
        
    return statistic,pvalue,average_sample_exp,average_sample_shl,significant_sign
#%% load data nrc 2010
ascii_grid_veg_indexnrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_nrc.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_nrc = [np.min(ascii_grid_veg_indexnrc[:,0]),np.max(ascii_grid_veg_indexnrc[:,0]),np.min(ascii_grid_veg_indexnrc[:,1]),np.max(ascii_grid_veg_indexnrc[:,1])]
lat_long_degre_nrc = [39.978183, 40.025319, -105.58556840332032, -105.52620245393197]

# slope less than 30
ascii_grid_veg_index_lsnrc = ascii_grid_veg_indexnrc[(ascii_grid_veg_indexnrc[:,5]>-90)&(ascii_grid_veg_indexnrc[:,4]<30)]

#y = -0.007512x + 14.168475  # R² = 0.969091
elev_nrc = ascii_grid_veg_index_lsnrc[:,2]
tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475 

veg_snow_temp_nrc = np.vstack([ascii_grid_veg_index_lsnrc[:,0],ascii_grid_veg_index_lsnrc[:,1].astype(int),ascii_grid_veg_index_lsnrc[:,2],
                               ascii_grid_veg_index_lsnrc[:,3],ascii_grid_veg_index_lsnrc[:,4],ascii_grid_veg_index_lsnrc[:,5],
                               ascii_grid_veg_index_lsnrc[:,6],tempDJF_gm_nrc,tempDJF_gm_nrc]).T

veg_snow_temp_dfnrc = pd.DataFrame(veg_snow_temp_nrc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
elevationBand_nrc = np.arange(min(veg_snow_temp_dfnrc['z']),max(veg_snow_temp_dfnrc['z']),65)
veg_snow_temp_clsfnrc = classificationElevation(veg_snow_temp_nrc,elevationBand_nrc,elev_class_inx)
aaaa_test2 = veg_snow_temp_clsfnrc[0:710,:]

meanTemp_elevClass_nrc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfnrc,elev_class_inx)   

#sampling
sample_size = 100
fsca_ut_nrc, fsca_0p_nrc = fSCAcalculation_elevClassific(veg_snow_temp_clsfnrc,elev_class_inx)
fsca_ut_samp_nrc, fsca_0p_samp_nrc = fSCAcalculation_elevRandom(veg_snow_temp_clsfnrc,elev_class_inx,sample_size)
statistic_nrc1,pvalue_nrc1,average_sample_exp_nrc1,average_sample_shl_nrc1,significant_nrc1 = wilcoxonTest(fsca_ut_samp_nrc, fsca_0p_samp_nrc)

# under canopy and 0pen fsca // exposed vs. sheltered
fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfnrc,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_nrc, fSCA_0u_shl_nrc = fsca0pn_fscaUnTr(fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc, elev_class_inx)
#sampling
randomNum = np.random.randint(0,sample_size,size=len(veg_snow_temp_clsfnrc[:,0]))
fSCA_0u_exp_nrc_samp, fSCA_0u_shl_nrc_samp = fSCAcalculation_elevNorthRandom(veg_snow_temp_clsfnrc,elev_class_inx,sample_size)
statistic_nrc2,pvalue_nrc2,average_sample_exp_nrc2,average_sample_shl_nrc2,significant_nrc2 = wilcoxonTest(fSCA_0u_exp_nrc_samp, fSCA_0u_shl_nrc_samp)

# December
#y = -0.00777x + 17.81629
tempDec_nrc = -0.00777*elevationBand_nrc + 17.81629

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
fsca_ut_jmz, fsca_0p_jmz = fSCAcalculation_elevClassific(veg_snow_temp_clsfJmz,elev_class_inx)
fsca_ut_samp_jmz, fsca_0p_samp_jmz = fSCAcalculation_elevRandom(veg_snow_temp_clsfJmz,elev_class_inx,sample_size)
#statistical test
statistic_jmz1,pvalue_jmz1,average_sample_exp_jmz1,average_sample_shl_jmz1,significant_jmz1 = wilcoxonTest(fsca_ut_samp_jmz, fsca_0p_samp_jmz)

# under canopy and 0pen fsca exposed vs. open
fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfJmz,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz = fsca0pn_fscaUnTr(fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz, elev_class_inx)
#statistical test
fSCA_0u_exp_Jmz_samp, fSCA_0u_shl_Jmz_samp = fSCAcalculation_elevNorthRandom(veg_snow_temp_clsfJmz,elev_class_inx,sample_size)
statistic_jmz2,pvalue_jmz2,average_sample_exp_jmz2,average_sample_shl_jmz2,significant_jmz2 = wilcoxonTest(fSCA_0u_exp_Jmz_samp, fSCA_0u_shl_Jmz_samp)

# y = -0.0043x + 5.6707
tempDec_jmz = -0.0048*elevationBand_jmz + 8.6707
#%% load data sagehen creek 26 March 2016
ascii_grid_veg_index26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc26m.npy')
#calculating topo_dimension from cut northness

# slope less than 30
ascii_grid_veg_index_ls26m = ascii_grid_veg_index26m[(ascii_grid_veg_index26m[:,5]>-90)&(ascii_grid_veg_index26m[:,4]<30)]

elev_sc26m = ascii_grid_veg_index_ls26m[:,2]
tempDJF_sc26m = -0.0016 * elev_sc26m + 1.802

veg_snow_temp_sc26m = np.vstack([ascii_grid_veg_index_ls26m[:,0],ascii_grid_veg_index_ls26m[:,1].astype(int),ascii_grid_veg_index_ls26m[:,2],
                                 ascii_grid_veg_index_ls26m[:,3],ascii_grid_veg_index_ls26m[:,4],ascii_grid_veg_index_ls26m[:,5],
                                 ascii_grid_veg_index_ls26m[:,6],tempDJF_sc26m,tempDJF_sc26m]).T

veg_snow_temp_sc_df26m = pd.DataFrame(veg_snow_temp_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc26m[700780:710057,:]

elevationBand_sc = np.arange(min(veg_snow_temp_sc_df26m['z']),max(veg_snow_temp_sc_df26m['z']),67)
veg_snow_temp_clsfsc26m = classificationElevation(veg_snow_temp_sc26m,elevationBand_sc,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc26m,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc26m, fsca_0p_sc26m = fSCAcalculation_elevClassific(veg_snow_temp_clsfsc26m,elev_class_inx)
fsca_ut_samp_sc26m, fsca_0p_samp_sc26m = fSCAcalculation_elevRandom(veg_snow_temp_clsfsc26m,elev_class_inx,sample_size)
#statistical test
statistic_sc26m1,pvalue_sc26m1,average_sample_exp_sc26m1,average_sample_shl_sc26m1,significant_sc26m1 = wilcoxonTest(fsca_ut_samp_sc26m, fsca_0p_samp_sc26m)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc26m, fSCA_ut_shl_sc26m, fSCA_0p_exp_sc26m, fSCA_0p_shl_sc26m = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfsc26m,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc26m, fSCA_0u_shl_sc26m = fsca0pn_fscaUnTr(fSCA_ut_exp_sc26m, fSCA_ut_shl_sc26m, fSCA_0p_exp_sc26m, fSCA_0p_shl_sc26m, elev_class_inx)
#statistical test
fSCA_0u_exp_sc26m_samp, fSCA_0u_shl_sc26m_samp = fSCAcalculation_elevNorthRandom(veg_snow_temp_clsfsc26m,elev_class_inx,sample_size)
statistic_sc26m2,pvalue_sc26m2,average_sample_exp_sc26m2,average_sample_shl_sc26m2,significant_sc26m2 = wilcoxonTest(fSCA_0u_exp_sc26m_samp, fSCA_0u_shl_sc26m_samp)

#y = y = -0.0013047442x + 2.0480940252
tempDec_sc26m = -0.0013047442 * elevationBand_sc + + 2.0480940252
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

elevationBand_sc_17a = np.arange(min(veg_snow_temp_sc_df17a['z']),max(veg_snow_temp_sc_df17a['z']),67)
veg_snow_temp_clsfsc17a = classificationElevation(veg_snow_temp_sc17a,elevationBand_sc_17a,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc17a,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc17a, fsca_0p_sc17a = fSCAcalculation_elevClassific(veg_snow_temp_clsfsc17a,elev_class_inx)
fsca_ut_samp_sc17a, fsca_0p_samp_sc17a = fSCAcalculation_elevRandom(veg_snow_temp_clsfsc17a,elev_class_inx,sample_size)
#statistical test
statistic_sc17a1,pvalue_sc17a1,average_sample_exp_sc17a1,average_sample_shl_sc17a1,significant_sc17a1 = wilcoxonTest(fsca_ut_samp_sc17a, fsca_0p_samp_sc17a)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfsc17a,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a = fsca0pn_fscaUnTr(fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a, elev_class_inx)
#statistical test
fSCA_0u_exp_sc17a_samp, fSCA_0u_shl_sc17a_samp = fSCAcalculation_elevNorthRandom(veg_snow_temp_clsfsc17a,elev_class_inx,sample_size)
statistic_sc17a2,pvalue_sc17a2,average_sample_exp_sc17a2,average_sample_shl_sc17a2,significant_sc17a2 = wilcoxonTest(fSCA_0u_exp_sc17a_samp, fSCA_0u_shl_sc17a_samp)

#y = -0.0013x + 2.7982
tempDec_sc17a = -0.0013 * elevationBand_sc_17a + 2.7982

#%% load data sagehen creek 17 May 2016
ascii_grid_veg_index18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc24May.npy')
#calculating topo_dimension from cut northness

# slope less than 30
ascii_grid_veg_index_ls18m = ascii_grid_veg_index18m[(ascii_grid_veg_index18m[:,5]>-90)&(ascii_grid_veg_index18m[:,4]<30)]

elev_sc18m = ascii_grid_veg_index_ls18m[:,2]
tempDJF_sc18m = -0.0016 * elev_sc18m + 1.802

veg_snow_temp_sc18m = np.vstack([ascii_grid_veg_index_ls18m[:,0],ascii_grid_veg_index_ls18m[:,1].astype(int),ascii_grid_veg_index_ls18m[:,2],
                                 ascii_grid_veg_index_ls18m[:,3],ascii_grid_veg_index_ls18m[:,4],ascii_grid_veg_index_ls18m[:,5],
                                 ascii_grid_veg_index_ls18m[:,6],tempDJF_sc18m,tempDJF_sc18m]).T

veg_snow_temp_sc_df18m = pd.DataFrame(veg_snow_temp_sc18m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc18m[700780:710057,:]

elevationBand_sc_18m = np.arange(min(veg_snow_temp_sc_df18m['z']),max(veg_snow_temp_sc_df18m['z']),67)
veg_snow_temp_clsfsc18m = classificationElevation(veg_snow_temp_sc18m,elevationBand_sc_18m,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfsc18m,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc18m, fsca_0p_sc18m = fSCAcalculation_elevClassific(veg_snow_temp_clsfsc18m,elev_class_inx)
fsca_ut_samp_sc18m, fsca_0p_samp_sc18m = fSCAcalculation_elevRandom(veg_snow_temp_clsfsc18m,elev_class_inx,sample_size)
#statistical test
statistic_sc18m1,pvalue_sc18m1,average_sample_exp_sc18m1,average_sample_shl_sc18m1,significant_sc18m1 = wilcoxonTest(fsca_ut_samp_sc18m, fsca_0p_samp_sc18m)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc18m, fSCA_ut_shl_sc18m, fSCA_0p_exp_sc18m, fSCA_0p_shl_sc18m = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfsc18m,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc18m, fSCA_0u_shl_sc18m = fsca0pn_fscaUnTr(fSCA_ut_exp_sc18m, fSCA_ut_shl_sc18m, fSCA_0p_exp_sc18m, fSCA_0p_shl_sc18m, elev_class_inx)
#statistical test
fSCA_0u_exp_sc18m_samp, fSCA_0u_shl_sc18m_samp = fSCAcalculation_elevNorthRandom(veg_snow_temp_clsfsc18m,elev_class_inx,sample_size)
statistic_sc18m2,pvalue_sc18m2,average_sample_exp_sc18m2,average_sample_shl_sc18m2,significant_sc18m2 = wilcoxonTest(fSCA_0u_exp_sc18m_samp, fSCA_0u_shl_sc18m_samp)

#y = -0.0020x + 5.2481
tempDec_sc18m = -0.002 * elevationBand_sc_18m + 5.2481
#%% load data KREW 2010
ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(ascii_grid_veg_indexKrew[:,0]),np.max(ascii_grid_veg_indexKrew[:,0]),np.min(ascii_grid_veg_indexKrew[:,1]),np.max(ascii_grid_veg_indexKrew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

# slope less than 30
ascii_grid_veg_index_lsKrew = ascii_grid_veg_indexKrew[(ascii_grid_veg_indexKrew[:,5]>-90)&(ascii_grid_veg_indexKrew[:,4]<30)]

elev_krew = ascii_grid_veg_index_lsKrew[:,2]
tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339 #R² = 0.5074

veg_snow_temp_Krew = np.vstack([ascii_grid_veg_index_lsKrew[:,0],ascii_grid_veg_index_lsKrew[:,1].astype(int),ascii_grid_veg_index_lsKrew[:,2],
                                ascii_grid_veg_index_lsKrew[:,3],ascii_grid_veg_index_lsKrew[:,4],ascii_grid_veg_index_lsKrew[:,5],
                                ascii_grid_veg_index_lsKrew[:,6],tempDJF_gm_krew,tempDJF_gm_krew]).T

veg_snow_temp_dfKrew = pd.DataFrame(veg_snow_temp_Krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand_krew = np.arange(min(veg_snow_temp_dfKrew['z']),max(veg_snow_temp_dfKrew['z']),67)
veg_snow_temp_clsfKrew = classificationElevation(veg_snow_temp_Krew,elevationBand_krew,elev_class_inx)

meanTemp_elevClass_krew = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfKrew,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_krew, fsca_0p_krew = fSCAcalculation_elevClassific(veg_snow_temp_clsfKrew,elev_class_inx)
fsca_ut_samp_krew, fsca_0p_samp_krew = fSCAcalculation_elevRandom(veg_snow_temp_clsfKrew,elev_class_inx,sample_size)
#statistical test
statistic_krew1,pvalue_krew1,average_sample_exp_krew1,average_sample_shl_krew1,significant_krew1 = wilcoxonTest(fsca_ut_samp_krew, fsca_0p_samp_krew)

# under canopy and 0pen fsca----exposed vs. sheltered
fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfKrew,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_Krew, fSCA_0u_shl_Krew = fsca0pn_fscaUnTr(fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew, elev_class_inx)
#statistical test
sample_size2 = 100      
fSCA_0u_exp_Krew_samp, fSCA_0u_shl_Krew_samp = fSCAcalculation_elevNorthRandom(veg_snow_temp_clsfKrew,elev_class_inx,sample_size2)
statistic_Krew2,pvalue_Krew2,average_sample_exp_Krew2,average_sample_shl_Krew2,significant_krew2 = wilcoxonTest(fSCA_0u_exp_Krew_samp, fSCA_0u_shl_Krew_samp)

#y = -0.003300x + 7.287576
tempDec_krew = -0.0033*elevationBand_krew + 7.287576 
#%% ploting fSCA vs. Dec temp lapse rate
fSCA_0u = [fsca_ut_sc26m, fsca_0p_sc26m, fsca_ut_sc17a, fsca_0p_sc17a,
           fsca_ut_sc18m, fsca_0p_sc18m, fsca_ut_krew, fsca_0p_krew,
           fsca_ut_jmz, fsca_0p_jmz, fsca_ut_nrc, fsca_0p_nrc]

meanTempDec = [tempDec_sc26m,tempDec_sc26m,tempDec_sc17a,tempDec_sc17a,
               tempDec_sc18m,tempDec_sc18m,tempDec_krew,tempDec_krew,
               tempDec_jmz,tempDec_jmz,tempDec_nrc,tempDec_nrc]

label = ['SCW 26MAR2016-UnderTree','SCW 26MAR2016-0pen','SCW 17APR2016-UnderTree','SCW 17APR2016-0pen',
         'SCW 18MAY2016-UnderTree','SCW 18MAY2016-0pen','Krew 2010-UnderTree','Krew 2010-0pen',
         'JMZ 2010-UnderTree','JMZ 2010 - 0pen','NRC 2010-UnderTree','NRC 2010 - 0pen']
color = ['darkred','gold','darkred','gold','darkred','gold',
         'red','orange','darkgreen','lightgreen','navy','lightblue']
marker = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

location = ['lower left','lower left','lower left','lower left','lower left','lower left']
xlimit_fsca = [(-1.5,1.5),(-1.5,1.5),(-1.5,1.5),(0.5,3),(-8,-1.5),(-8,-1.5)]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
    for pnt in range (10):
        plt.scatter(meanTempDec[2*i][pnt],fSCA_0u[2*i][pnt], s=25**2, color = color[2*i], marker = marker[i][pnt]) #label = label[2*i],
        plt.scatter(meanTempDec[2*i+1][pnt],fSCA_0u[2*i+1][pnt],  s=25**2, color = color[2*i+1], marker = marker[i][pnt]) #label = label[2*i+1],

    plt.ylabel('fSCA%', fontsize=25)
    plt.yticks(fontsize=15)
    plt.xlabel('average temp from Dec till date of flights(C)', fontsize=20)
    plt.xticks(fontsize=15)
    
    plt.legend([label[2*i], label[2*i+1]], fontsize=20, loc = location[i])
    plt.ylim((0,110))  # adjust the top leaving bottom unchanged
    plt.xlim(xlimit_fsca[i])

plt.title('fSCA under canopy and in open sites in different temp lapse rate in 4 sites', fontsize=35, y=2.25, x=-1) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_dec_fsca_all2.png')
#%% ploting fsca vs. DJF temp lapse rate
fSCA_0u = [fsca_ut_sc26m, fsca_0p_sc26m, fsca_ut_sc17a, fsca_0p_sc17a,
           fsca_ut_sc18m, fsca_0p_sc18m, fsca_ut_krew, fsca_0p_krew,
           fsca_ut_jmz, fsca_0p_jmz, fsca_ut_nrc, fsca_0p_nrc]

meanTemp = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
            meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
            meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label = ['SCW 26MAR2016-UnderTree','SCW 26MAR2016-0pen','SCW 17APR2016-UnderTree','SCW 17APR2016-0pen',
         'SCW 18MAY2016-UnderTree','SCW 18MAY2016-0pen','Krew 2010-UnderTree','Krew 2010-0pen',
         'JMZ 2010-UnderTree','JMZ 2010 - 0pen','NRC 2010-UnderTree','NRC 2010 - 0pen']
color = ['darkred','gold','darkred','gold','darkred','gold',
         'red','orange','darkgreen','lightgreen','navy','lightblue']
marker = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

location = ['lower left','lower left','uper right','lower left','lower left','lower left']
xlimit_fsca_DJF = [(-2.5,-1),(-2.5,-1),(-2.5,-1),(0.5,3),(-9,-3),(-10,-5)]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
    for pnt in range (10):
        plt.scatter(meanTemp[2*i][pnt],fSCA_0u[2*i][pnt], s=25**2, color = color[2*i], marker = marker[i][pnt]) #label = label[2*i],
        plt.scatter(meanTemp[2*i+1][pnt],fSCA_0u[2*i+1][pnt],  s=25**2, color = color[2*i+1], marker = marker[i][pnt]) #label = label[2*i+1],

    plt.ylabel('fSCA%', fontsize=25)
    plt.yticks(fontsize=15)
    plt.xlabel('average temp of DJF (C)', fontsize=25)
    plt.xticks(fontsize=15)
    
    plt.legend([label[2*i], label[2*i+1]], fontsize=20, loc = location[i])
    plt.ylim((0,110))  # adjust the top leaving bottom unchanged
    plt.xlim(xlimit_fsca_DJF[i])

plt.title('fSCA under canopy and in open sites in different temp lapse rate in 4 sites', fontsize=35, y=2.25, x=-1) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_jfd_fsca_all2.png')

#%% ploting DJF temp vs terrain features delta fsca
fSCA_d0u = [fSCA_0u_exp_sc26m, fSCA_0u_shl_sc26m, fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a, 
            fSCA_0u_exp_sc18m, fSCA_0u_shl_sc18m, fSCA_0u_exp_Krew, fSCA_0u_shl_Krew, 
            fSCA_0u_exp_nrc, fSCA_0u_shl_nrc, fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz]
meanTemp2 = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
             meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label = ['SCW 26Mar2016 -- exposed','SCW 26Mar2016 -- sheltered',
         'SCW 17Apr2016 -- exposed','SCW 17Apr2016 -- sheltered',
         'SCW 18May2016 -- exposed','SCW 18May2016 -- sheltered',
         'Krew 2010 -- exposed','Krew 2010 -- sheltered','JMZ 2010 - exposed',
         'JMZ 2010 - sheltered','NRC 2010 - exposed','NRC 2010 - sheltered']

color = ['darkred','gold','darkred','gold','darkred','gold',
         'red','orange','darkgreen','lightgreen','navy','lightblue']

marker = [significant_sc26m2,significant_sc17a2,significant_sc18m2,significant_krew2,
          significant_nrc2,significant_jmz2]

location = ['lower left','lower left','lower left','lower left','lower left','lower left']
xlimit_deltAfsca = [(-2.5,-1),(-2.5,-1),(-2.5,-1),(0,3),(-9,-3),(-10.5,-5)]
arrow_loc = [-2.4, -2.4, -2.4, 0.3, -8.5, -10]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
    for pnt in range (10):
        plt.scatter(meanTemp2[2*i][pnt],fSCA_d0u[2*i][pnt], s=25**2, color = color[2*i], marker = marker[i][pnt]) #label = label[2*i],
        plt.scatter(meanTemp2[2*i+1][pnt],fSCA_d0u[2*i+1][pnt],  s=25**2, color = color[2*i+1], marker = marker[i][pnt]) #label = label[2*i+1],

    plt.ylabel('(fSCA_open - fSCA_under) / fSCA_open', fontsize=20)
    plt.yticks(fontsize=15)
    plt.xlabel('average temp of DJF (C)', fontsize=25)
    plt.xticks(fontsize=15)
    
    plt.legend([label[2*i], label[2*i+1]], fontsize=15, loc = location[i])
    plt.ylim((-1.2,0.8))
    plt.xlim(xlimit_deltAfsca[i])
    
    x = [-11,-8,-5,-2,0,2,4]
    y = [0,0,0,0,0,0,0]
    plt.plot(x,y, color = 'black') #line_temp_ls[2*i],
    
    #plt.arrow(arrow_loc[i], 0.05, 0, 0.5, fc="k", ec="k", head_width=0.1, head_length=0.1)

plt.title('(delta fSCA)/fSCA_op in different temp lapse rate and northness in 4 sites', fontsize=35, y=2.25, x=-1) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrth_deltaFsca_all3.png')

#for i in range (4):
#    plt.subplot(221+i)
#    plt.scatter(meanTemp[2*i],fSCA_0u[2*i], label = label[2*i], s=20**2, color = color[2*i])
#    plt.scatter(meanTemp[2*i+1],fSCA_0u[2*i+1], label = label[2*i+1], s=20**2, color = color[2*i+1])
#
#    plt.ylabel('delta fSCA% / fSCA_open', fontsize=25)
#    plt.yticks(fontsize=15)
#    plt.xlabel('average temp of DJF (C)', fontsize=25)
#    plt.xticks(fontsize=15)
#    plt.legend(fontsize=20, loc = location[i])














