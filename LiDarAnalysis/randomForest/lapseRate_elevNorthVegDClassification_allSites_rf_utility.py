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
from matplotlib import colors
import matplotlib.patches as mpatches
import seaborn as sns
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
# functions
def classificationElevation(points_df_ls,elevationBand,elev_class_inx):#0.25
    
    for indx in range (len(elev_class_inx)):
        points_df_ls[:,9][(points_df_ls[:,2]>=elevationBand[indx])&(points_df_ls[:,2]<elevationBand[indx+1])]=elev_class_inx[indx]
    
    points_df_ls[:,9][(points_df_ls[:,9]<=100)]=110    
    
    return points_df_ls#, pointFile_df
               
def fSCAcalculation_elevClassific(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,9]==101 contains elevation classification

    fSCA_unTr = []
    fSCA_0pen = []
    
    for indx in range (len(elev_class_inx)): 
        #under tree with snow, exposed
        ascii_ut_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])])
        ascii_ut_nSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])])
        #fSCA under tree
        fSCA_ut = 100.*ascii_ut_snow/(ascii_ut_snow+ascii_ut_nSnow)
        fSCA_unTr.append(fSCA_ut)
    
        #0pen with snow
        ascii_0p_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])])
        #0pen no snow
        ascii_0p_nSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])])
        #fSCA 0pen
        fSCA_0p = 100.*ascii_0p_snow/(ascii_0p_nSnow+ascii_0p_snow)
        fSCA_0pen.append(fSCA_0p)
    
    return fSCA_unTr, fSCA_0pen

def fSCAcalculation_elevNorthClassific(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,9]==101 contains elevation classification

    fSCA_unTr_exp = []
    fSCA_unTr_shl = []
    fSCA_0pen_exp = []
    fSCA_0pen_shl = []
    
    for indx in range (len(elev_class_inx)): 
        #under tree with snow, exposed
        ascii_ut_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #under tree no snow, exposed
        ascii_ut_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #under tree with snow, sheltered
        ascii_ut_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
        #under tree no snow, sheltered
        ascii_ut_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])

        #fSCA under tree, exposed
        fSCA_ut_exp = 100.*ascii_ut_snow_exp/(ascii_ut_snow_exp+ascii_ut_noSnow_exp)
        #fSCA under tree, sheltered
        fSCA_ut_shl = 100.*ascii_ut_snow_shl/(ascii_ut_snow_shl+ascii_ut_nSnow_shl)
    
        fSCA_unTr_exp.append(fSCA_ut_exp)
        fSCA_unTr_shl.append(fSCA_ut_shl)
    
        #0pen with snow, exposed
        ascii_0p_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #0pen no snow, exposed
        ascii_0p_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
        #0pen with snow, sheltered
        ascii_0p_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
        #0pen no snow, sheltered
        ascii_0p_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])

        #fSCA 0pen, exposed
        fSCA_0p_exp = 100.*ascii_0p_snow_exp/(ascii_0p_snow_exp+ascii_0p_noSnow_exp)
        #fSCA 0pen, sheltered
        fSCA_0p_shl = 100.*ascii_0p_snow_shl/(ascii_0p_snow_shl+ascii_0p_nSnow_shl)
    
        fSCA_0pen_exp.append(fSCA_0p_exp)
        fSCA_0pen_shl.append(fSCA_0p_shl)
    
    return fSCA_unTr_exp, fSCA_unTr_shl, fSCA_0pen_exp, fSCA_0pen_shl


def fSCAcalculation_elevNorthVegDensClassific(veg_snow_temp_elev_clsf_file,elev_class_inx): #veg_snow_temp_elev_clsf_file[:,9]==101 contains elevation classification

    fSCA_unTr_exp_l = []
    fSCA_unTr_exp_h = []

    fSCA_unTr_shl_l = []
    fSCA_unTr_shl_h = []
    
    for indx in range (len(elev_class_inx)): 
        #under tree with snow, exposed
        ascii_ut_snow_exp_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)])
        ascii_ut_snow_exp_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)])
        #under tree no snow, exposed
        ascii_ut_noSnow_exp_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)])
        ascii_ut_noSnow_exp_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)])
        #under tree with snow, sheltered
        ascii_ut_snow_shl_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)])
        ascii_ut_snow_shl_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)])
        #under tree no snow, sheltered
        ascii_ut_nSnow_shl_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)])
        ascii_ut_nSnow_shl_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)])

        #fSCA under tree, exposed
        fSCA_ut_exp_l = 100.*ascii_ut_snow_exp_l/(ascii_ut_snow_exp_l+ascii_ut_noSnow_exp_l)
        fSCA_ut_exp_h = 100.*ascii_ut_snow_exp_h/(ascii_ut_snow_exp_h+ascii_ut_noSnow_exp_h)

        #fSCA under tree, sheltered
        fSCA_ut_shl_l = 100.*ascii_ut_snow_shl_l/(ascii_ut_snow_shl_l+ascii_ut_nSnow_shl_l)
        fSCA_ut_shl_h = 100.*ascii_ut_snow_shl_h/(ascii_ut_snow_shl_h+ascii_ut_nSnow_shl_h)
    
        fSCA_unTr_exp_l.append(fSCA_ut_exp_l)
        fSCA_unTr_shl_l.append(fSCA_ut_shl_l)
        fSCA_unTr_exp_h.append(fSCA_ut_exp_h)
        fSCA_unTr_shl_h.append(fSCA_ut_shl_h)
    
    return [fSCA_unTr_exp_l, fSCA_unTr_shl_l, fSCA_unTr_exp_h, fSCA_unTr_shl_h]


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
        temp_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])][:,7])
        meanTemp.append(temp_mean)    
        
    return meanTemp

def fSCAcalculation_elevRandom(veg_snow_temp_elev_clsf_file,elev_class_inx,sample_size): #veg_snow_temp_elev_clsf_file[:,9]==101 contains elevation classification

    randomNum = np.random.randint(0,sample_size,size=len(veg_snow_temp_elev_clsf_file[:,0]))
    veg_snow_temp_elev_clsf_file[:,6]=randomNum
    
    fSCA_0p = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2} 
    fSCA_ut = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
    
    for indx in range (len(elev_class_inx)): 
        
        fSCA_0p_rand = np.zeros(sample_size) 
        fSCA_ut_rand = np.zeros(sample_size)
        
        for rand in range (sample_size):
            
            #under tree with snow
            ascii_ut_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow
            ascii_ut_noSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #fSCA under tree
            fSCA_ut_each = 100.*ascii_ut_snow/(ascii_ut_snow+ascii_ut_noSnow)
            
            #0pen with snow
            ascii_0p_snow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen no snow
            ascii_0p_noSnow = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #fSCA 0pen
            fSCA_0p_each = 100.*ascii_0p_snow/(ascii_0p_snow+ascii_0p_noSnow)
            

            fSCA_ut_rand[rand] = fSCA_ut_each
            fSCA_0p_rand[rand] = fSCA_0p_each
    
        fSCA_0p[indx] = fSCA_ut_rand
        fSCA_ut[indx] = fSCA_0p_rand
        
    return fSCA_ut, fSCA_0p


def fSCAcalculation_elevNorthRandom(veg_snow_temp_elev_clsf_file,elev_class_inx,sample_size): #veg_snow_temp_elev_clsf_file[:,9]==101 contains elevation classification

    randomNum = np.random.randint(0,sample_size,size=len(veg_snow_temp_elev_clsf_file[:,0]))
    veg_snow_temp_elev_clsf_file[:,6]=randomNum
    
    fSCA_0u_exp = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2} 
    fSCA_0u_shl = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
    
    for indx in range (len(elev_class_inx)): 
        
        fSCA_0u_exp_rand = np.zeros(sample_size) 
        fSCA_0u_shl_rand = np.zeros(sample_size)
        
        for rand in range (sample_size):
            
            #under tree with snow, exposed
            ascii_ut_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow, exposed
            ascii_ut_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree with snow, sheltered
            ascii_ut_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow, sheltered
            ascii_ut_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            
            #fSCA under tree, exposed
            #fSCA_ut_exp = 100.*ascii_ut_snow_exp/(ascii_ut_snow_exp+ascii_ut_noSnow_exp)
            #fSCA under tree, sheltered
            #fSCA_ut_shl = 100.*ascii_ut_snow_shl/(ascii_ut_snow_shl+ascii_ut_nSnow_shl)
            
            #0pen with snow, exposed
            ascii_0p_snow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen no snow, exposed
            ascii_0p_noSnow_exp = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen with snow, sheltered
            ascii_0p_snow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #0pen no snow, sheltered
            ascii_0p_nSnow_shl = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])

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


def fSCAcalculation_elevNorthVegDensRandom(veg_snow_temp_elev_clsf_file,elev_class_inx,sample_size): #veg_snow_temp_elev_clsf_file[:,9]==101 contains elevation classification

    randomNum = np.random.randint(0,sample_size,size=len(veg_snow_temp_elev_clsf_file[:,0]))
    veg_snow_temp_elev_clsf_file[:,6]=randomNum
    
    fSCA_ut_exp_l = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2} 
    fSCA_ut_shl_l = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
    fSCA_ut_exp_h = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2} 
    fSCA_ut_shl_h = {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2}
    
    for indx in range (len(elev_class_inx)): 
        
        fSCA_ut_exp_l_rand = np.zeros(sample_size) 
        fSCA_ut_shl_l_rand = np.zeros(sample_size)
        fSCA_ut_exp_h_rand = np.zeros(sample_size) 
        fSCA_ut_shl_h_rand = np.zeros(sample_size)
        
        for rand in range (sample_size):
            
            #under tree with snow, exposed
            ascii_ut_snow_exp_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            ascii_ut_snow_exp_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow, exposed
            ascii_ut_noSnow_exp_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            ascii_ut_noSnow_exp_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree with snow, sheltered
            ascii_ut_snow_shl_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            ascii_ut_snow_shl_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            #under tree no snow, sheltered
            ascii_ut_nSnow_shl_l = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]<0.4)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
            ascii_ut_nSnow_shl_h = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,9]==elev_class_inx[indx])&(veg_snow_temp_elev_clsf_file[:,3]>0.1)&(veg_snow_temp_elev_clsf_file[:,8]>0.6)&(veg_snow_temp_elev_clsf_file[:,6]==rand)])
        
            #low density
            #fSCA under canopy, exposed
            if ascii_ut_snow_exp_l+ascii_ut_noSnow_exp_l == 0:
                fSCA_ut_exp_l_each = 0
            else:            
                fSCA_ut_exp_l_each = 100.*ascii_ut_snow_exp_l/(ascii_ut_snow_exp_l+ascii_ut_noSnow_exp_l)
            #fSCA under canopy, sheltered
            if ascii_ut_snow_shl_l+ascii_ut_nSnow_shl_l == 0:
                fSCA_ut_shl_l_each = 0
            else:
                fSCA_ut_shl_l_each = 100.*ascii_ut_snow_shl_l/(ascii_ut_snow_shl_l+ascii_ut_nSnow_shl_l)

            #high density
            #fSCA under canopy, exposed
            if ascii_ut_snow_exp_h+ascii_ut_noSnow_exp_h == 0:
                fSCA_ut_exp_h_each = 0
            else:            
                fSCA_ut_exp_h_each = 100.*ascii_ut_snow_exp_h/(ascii_ut_snow_exp_h+ascii_ut_noSnow_exp_h)
            #fSCA under canopy, sheltered
            if ascii_ut_snow_shl_h+ascii_ut_nSnow_shl_h == 0:
                fSCA_ut_shl_h_each = 0
            else:
                fSCA_ut_shl_h_each = 100.*ascii_ut_snow_shl_h/(ascii_ut_snow_shl_h+ascii_ut_nSnow_shl_h)

            fSCA_ut_exp_l_rand[rand] = fSCA_ut_exp_l_each
            fSCA_ut_shl_l_rand[rand] = fSCA_ut_shl_l_each
            fSCA_ut_exp_h_rand[rand] = fSCA_ut_exp_h_each
            fSCA_ut_shl_h_rand[rand] = fSCA_ut_shl_h_each
    
        fSCA_ut_exp_l[indx] = fSCA_ut_exp_l_rand
        fSCA_ut_shl_l[indx] = fSCA_ut_shl_l_rand
        fSCA_ut_exp_h[indx] = fSCA_ut_exp_h_rand
        fSCA_ut_shl_h[indx] = fSCA_ut_shl_h_rand
        
    return fSCA_ut_exp_l, fSCA_ut_exp_h, fSCA_ut_shl_l, fSCA_ut_shl_h 

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
            significant = 'P'
        else:
            significant = '_'
        
        significant_sign.append(significant)
        
        average_sample_shl.append(np.average(fSCA_0u_shl_samp[iii]))
        average_sample_exp.append(np.average(fSCA_0u_exp_samp[iii]))
        
    return statistic,pvalue,average_sample_exp,average_sample_shl,significant_sign

def plotSnowMap(snow_temp_vegdens_allExtent,ouputPath,T): 
    
    lat_long_nad83 = [np.min(snow_temp_vegdens_allExtent[:,0]),np.max(snow_temp_vegdens_allExtent[:,0]),np.min(snow_temp_vegdens_allExtent[:,1]),np.max(snow_temp_vegdens_allExtent[:,1])]
    ncols = (lat_long_nad83[1]-lat_long_nad83[0]+1).astype(int)
    nrows = (lat_long_nad83[3]-lat_long_nad83[2]+1).astype(int)
    snwMap = np.reshape(snow_temp_vegdens_allExtent[:,5],[ncols,nrows])

    #plotting
    cmap = colors.ListedColormap(['darkgray', 'green', 'wheat', 'white', 'blue'])
    labels = {-99:'excluded',-7:'no snow under trees',-4:'no snow in 0pens',44:'snow in 0pens',77:'snow under trees'}
    colors_cmap = cmap(np.arange(cmap.N))
    cmap = {-99:colors_cmap[0], -7:colors_cmap[1], -4:colors_cmap[2], 44:colors_cmap[3], 77:colors_cmap[4]}

    arrayShow = np.array([[cmap[i] for i in j] for j in snwMap])    
    patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]

    elev_2d = (np.reshape(snow_temp_vegdens_allExtent[:,2],[ncols,nrows])).astype(int)

    plt.figure(figsize=(200,190))
    
    contours = plt.contour(elev_2d, 20, colors='black', linewidths= 6, extent=lat_long_nad83)
    plt.clabel(contours, inline=True, fontsize=100, fmt = '%1.0f')
    
    plt.imshow(arrayShow, extent=lat_long_nad83)#'gist_earth'
    plt.xticks(fontsize=170)
    plt.yticks(fontsize=170)

    plt.legend(handles=patches, loc=2, fontsize=150, borderaxespad=0.)
    plt.title('Snow presence map under forest canopy and open areas in %s' %T, fontsize=250, y=1.015)

    plt.savefig(ouputPath)

def plotSnowMap_northness(snow_temp_vegdens_allExtent,ouputPath,T): 
    snow_temp_vegdens_allExtent[:,10] = -99
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==77)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=771
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==77)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=770
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==44)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=441
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==44)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=440
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-7)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=-71
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-7)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=-70
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-4)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=-41
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-4)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=-40

    
    lat_long_nad83 = [np.min(snow_temp_vegdens_allExtent[:,0]),np.max(snow_temp_vegdens_allExtent[:,0]),np.min(snow_temp_vegdens_allExtent[:,1]),np.max(snow_temp_vegdens_allExtent[:,1])]
    ncols = (lat_long_nad83[1]-lat_long_nad83[0]+1).astype(int)
    nrows = (lat_long_nad83[3]-lat_long_nad83[2]+1).astype(int)
    snwMap = np.reshape(snow_temp_vegdens_allExtent[:,10],[ncols,nrows])

    #plotting
    cmap = colors.ListedColormap(['darkgray', 'lightgreen', 'darkgreen', 'rosybrown', 'saddlebrown', 'white', 'wheat','lightblue','darkblue'])
    
    labels = {-99:'excluded',-71:'no snow under trees--exposed',-70:'no snow under trees--sheltered',
              -41:'no snow in 0pens--exposed',-40:'no snow in 0pens--sheltered',441:'snow in 0pens--exposed',
              440:'snow in 0pens--sheltered',771:'snow under trees--exposed',770:'snow under trees--sheltered'}
    
    colors_cmap = cmap(np.arange(cmap.N))
    cmap = {-99:colors_cmap[0],-71:colors_cmap[1],-70:colors_cmap[2],-41:colors_cmap[3],-40:colors_cmap[4],
            441:colors_cmap[5],440:colors_cmap[6],771:colors_cmap[7],770:colors_cmap[8]}


    arrayShow = np.array([[cmap[i] for i in j] for j in snwMap])    
    patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]
    
    elev_2d = (np.reshape(snow_temp_vegdens_allExtent[:,2],[ncols,nrows])).astype(int)

    plt.figure(figsize=(210,200))
    
    contours = plt.contour(elev_2d, 30, colors='black', linewidths=5, extent=lat_long_nad83)
    plt.clabel(contours, inline=True, fontsize=100, fmt = '%1.0f')
    
    plt.imshow(arrayShow, extent=lat_long_nad83)#'gist_earth'
    plt.xticks(fontsize=180)
    plt.yticks(fontsize=180)

    plt.legend(handles=patches, loc=2, fontsize=120, borderaxespad=0.)
    plt.title('Snow presence map under forest canopy and open areas based on northness in %s' %T, fontsize=250, y=1.015)

    plt.savefig(ouputPath)

def plotSnowMap_lowResol(snow_temp_vegdens_allExtent,ouputPath,T): 
    
    lat_long_nad83 = [np.min(snow_temp_vegdens_allExtent[:,0]),np.max(snow_temp_vegdens_allExtent[:,0]),np.min(snow_temp_vegdens_allExtent[:,1]),np.max(snow_temp_vegdens_allExtent[:,1])]
    ncols = (lat_long_nad83[1]-lat_long_nad83[0]+1).astype(int)
    nrows = (lat_long_nad83[3]-lat_long_nad83[2]+1).astype(int)
    snwMap = np.reshape(snow_temp_vegdens_allExtent[:,5],[ncols,nrows])

    #plotting
    cmap = colors.ListedColormap(['darkgray', 'green', 'wheat', 'white', 'blue'])
    labels = {-99:'excluded',-7:'no snow under trees',-4:'no snow in 0pens',44:'snow in 0pens',77:'snow under trees'}
    colors_cmap = cmap(np.arange(cmap.N))
    cmap = {-99:colors_cmap[0], -7:colors_cmap[1], -4:colors_cmap[2], 44:colors_cmap[3], 77:colors_cmap[4]}

    arrayShow = np.array([[cmap[i] for i in j] for j in snwMap])    
    patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]

    elev_2d = (np.reshape(snow_temp_vegdens_allExtent[:,2],[ncols,nrows])).astype(int)

    plt.figure(figsize=(50,50))
    
    contours = plt.contour(elev_2d, 15, colors='black', linewidths= 4, extent=lat_long_nad83)#
    plt.clabel(contours, inline=True, fontsize=30, fmt = '%1.0f')
    
    plt.imshow(arrayShow, extent=lat_long_nad83)#'gist_earth'
    #lat_long_nad83_krew = [np.min(snow_temp_vegdens_allExtent[:,0]),np.max(snow_temp_vegdens_allExtent[:,0]),
    #                       np.min(snow_temp_vegdens_allExtent[:,1])+1000,np.max(snow_temp_vegdens_allExtent[:,1])-2000]
    #plt.imshow(arrayShow, extent=lat_long_nad83_krew)#'gist_earth'
    plt.xticks(fontsize=60, ha = 'left')#
    plt.yticks(fontsize=60)

    plt.legend(handles=patches, loc=2, fontsize=40, borderaxespad=0.)
    plt.title('Snow presence map under forest canopy and open areas in %s' %T, fontsize=70, y=1.015)

    plt.savefig(ouputPath)

def plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent,ouputPath,T): 
    snow_temp_vegdens_allExtent[:,10] = -99
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==77)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=771
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==77)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=770
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==44)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=441
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==44)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=440
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-7)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=-71
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-7)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=-70
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-4)&(snow_temp_vegdens_allExtent[:,3]<=-0.1)]=-41
    snow_temp_vegdens_allExtent[:,10][(snow_temp_vegdens_allExtent[:,5]==-4)&(snow_temp_vegdens_allExtent[:,3]>=0.1)]=-40

    
    lat_long_nad83 = [np.min(snow_temp_vegdens_allExtent[:,0]),np.max(snow_temp_vegdens_allExtent[:,0]),np.min(snow_temp_vegdens_allExtent[:,1]),np.max(snow_temp_vegdens_allExtent[:,1])]
    ncols = (lat_long_nad83[1]-lat_long_nad83[0]+1).astype(int)
    nrows = (lat_long_nad83[3]-lat_long_nad83[2]+1).astype(int)
    snwMap = np.reshape(snow_temp_vegdens_allExtent[:,10],[ncols,nrows])

    #plotting
    cmap = colors.ListedColormap(['darkgray', 'lightgreen', 'darkgreen', 'rosybrown', 'saddlebrown', 'white', 'wheat','lightblue','darkblue'])
    
    labels = {-99:'excluded',-71:'no snow under trees--exposed',-70:'no snow under trees--sheltered',
              -41:'no snow in 0pens--exposed',-40:'no snow in 0pens--sheltered',441:'snow in 0pens--exposed',
              440:'snow in 0pens--sheltered',771:'snow under trees--exposed',770:'snow under trees--sheltered'}
    
    colors_cmap = cmap(np.arange(cmap.N))
    cmap = {-99:colors_cmap[0],-71:colors_cmap[1],-70:colors_cmap[2],-41:colors_cmap[3],-40:colors_cmap[4],
            441:colors_cmap[5],440:colors_cmap[6],771:colors_cmap[7],770:colors_cmap[8]}


    arrayShow = np.array([[cmap[i] for i in j] for j in snwMap])    
    patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]
    
    elev_2d = (np.reshape(snow_temp_vegdens_allExtent[:,2],[ncols,nrows])).astype(int)

    plt.figure(figsize=(50,50))
    
    contours = plt.contour(elev_2d, 15, colors='black', linewidths=3, extent=lat_long_nad83)
    plt.clabel(contours, inline=True, fontsize=30, fmt = '%1.0f')

#    lat_long_nad83_krew = [np.min(snow_temp_vegdens_allExtent[:,0]),np.max(snow_temp_vegdens_allExtent[:,0]),
#                           np.min(snow_temp_vegdens_allExtent[:,1])+1000,np.max(snow_temp_vegdens_allExtent[:,1])-2000]
#    plt.imshow(arrayShow, extent=lat_long_nad83_krew)#'gist_earth'
    
    plt.imshow(arrayShow, extent=lat_long_nad83)#'gist_earth'
    plt.xticks(fontsize=60, ha = 'left')#
    plt.yticks(fontsize=60)

    plt.legend(handles=patches, loc=2, fontsize=30, borderaxespad=0.)
    plt.title('Snow presence map under forest canopy and open areas based on northness in %s' %T, fontsize=60, y=1.015)

    plt.savefig(ouputPath)

def fsca0p_fscaUt_fscaLimit (fsca_0p_sc26m,fsca_0p_sc17a,fsca_0p_sc18m,fsca_0p_krew,fsca_0p_jmz,fsca_0p_nrc,
                             fsca_ut_sc26m,fsca_ut_sc17a,fsca_ut_sc18m,fsca_ut_krew,fsca_ut_jmz,fsca_ut_nrc,limitFsca1,limitFsca2):
    
    fsca_0p_sc26m_arr = np.array(fsca_0p_sc26m); fsca_0p_sc17a_arr = np.array(fsca_0p_sc17a)
    fsca_0p_sc18m_arr = np.array(fsca_0p_sc18m); fsca_0p_krew_arr = np.array(fsca_0p_krew)
    fsca_0p_jmz_arr = np.array(fsca_0p_jmz); fsca_0p_nrc_arr = np.array(fsca_0p_nrc)

    fsca_ut_sc26m_arr = np.array(fsca_ut_sc26m); fsca_ut_sc17a_arr = np.array(fsca_ut_sc17a)
    fsca_ut_sc18m_arr = np.array(fsca_ut_sc18m); fsca_ut_krew_arr = np.array(fsca_ut_krew)
    fsca_ut_jmz_arr = np.array(fsca_ut_jmz); fsca_ut_nrc_arr = np.array(fsca_ut_nrc)

    #%82 open- under
    fsca_0p_sc26m_lim = np.mean(fsca_0p_sc26m_arr[(fsca_0p_sc26m_arr<=limitFsca2)&(fsca_0p_sc26m_arr>=limitFsca1)])
    fsca_0p_sc17a_lim = np.mean(fsca_0p_sc17a_arr[(fsca_0p_sc17a_arr<=limitFsca2)&(fsca_0p_sc17a_arr>=limitFsca1)])
    fsca_0p_sc18m_lim = np.mean(fsca_0p_sc18m_arr[(fsca_0p_sc18m_arr<=limitFsca2)&(fsca_0p_sc18m_arr>=limitFsca1)])
    fsca_0p_krew_lim = np.mean(fsca_0p_krew_arr[(fsca_0p_krew_arr<=limitFsca2)&(fsca_0p_krew_arr>=limitFsca1)])
    fsca_0p_jmz_lim = np.mean(fsca_0p_jmz_arr[(fsca_0p_jmz_arr<=limitFsca2)&(fsca_0p_jmz_arr>=limitFsca1)])
    fsca_0p_nrc_lim = np.mean(fsca_0p_nrc_arr[(fsca_0p_nrc_arr<=limitFsca2)&(fsca_0p_nrc_arr>=limitFsca1)])

    fsca_ut_sc26m_lim = np.mean(fsca_ut_sc26m_arr[(fsca_ut_sc26m_arr<=limitFsca2)&(fsca_ut_sc26m_arr>=limitFsca1)])
    fsca_ut_sc17a_lim = np.mean(fsca_ut_sc17a_arr[(fsca_ut_sc17a_arr<=limitFsca2)&(fsca_ut_sc17a_arr>=limitFsca1)])
    fsca_ut_sc18m_lim = np.mean(fsca_ut_sc18m_arr[(fsca_ut_sc18m_arr<=limitFsca2)&(fsca_ut_sc18m_arr>=limitFsca1)])
    fsca_ut_krew_lim = np.mean(fsca_ut_krew_arr[(fsca_ut_krew_arr<=limitFsca2)&(fsca_ut_krew_arr>=limitFsca1)])
    fsca_ut_jmz_lim = np.mean(fsca_ut_jmz_arr[(fsca_ut_jmz_arr<=limitFsca2)&(fsca_ut_jmz_arr>=limitFsca1)])
    fsca_ut_nrc_lim = np.mean(fsca_ut_nrc_arr[(fsca_ut_nrc_arr<=limitFsca2)&(fsca_ut_nrc_arr>=limitFsca1)])

    fsca_op = np.array([fsca_0p_sc26m_lim,fsca_0p_sc17a_lim,fsca_0p_sc18m_lim,fsca_0p_krew_lim,fsca_0p_jmz_lim,fsca_0p_nrc_lim])
    fsca_ut = np.array([fsca_ut_sc26m_lim,fsca_ut_sc17a_lim,fsca_ut_sc18m_lim,fsca_ut_krew_lim,fsca_ut_jmz_lim,fsca_ut_nrc_lim])
    fsca_0t = fsca_op-fsca_ut
    
    return fsca_0t

def mean_fsca_Limit (fsca_0p_sc26m,fsca_0p_sc17a,fsca_0p_sc18m,fsca_0p_krew,fsca_0p_jmz,fsca_0p_nrc,
                     limitFsca1,limitFsca2):
    
    fsca_0p_sc26m_arr = np.array(fsca_0p_sc26m); fsca_0p_sc17a_arr = np.array(fsca_0p_sc17a)
    fsca_0p_sc18m_arr = np.array(fsca_0p_sc18m); fsca_0p_krew_arr = np.array(fsca_0p_krew)
    fsca_0p_jmz_arr = np.array(fsca_0p_jmz); fsca_0p_nrc_arr = np.array(fsca_0p_nrc)

    fsca_0p_sc26m_lim = np.mean(fsca_0p_sc26m_arr[(fsca_0p_sc26m_arr<=limitFsca2)&(fsca_0p_sc26m_arr>=limitFsca1)])
    fsca_0p_sc17a_lim = np.mean(fsca_0p_sc17a_arr[(fsca_0p_sc17a_arr<=limitFsca2)&(fsca_0p_sc17a_arr>=limitFsca1)])
    fsca_0p_sc18m_lim = np.mean(fsca_0p_sc18m_arr[(fsca_0p_sc18m_arr<=limitFsca2)&(fsca_0p_sc18m_arr>=limitFsca1)])
    fsca_0p_krew_lim = np.mean(fsca_0p_krew_arr[(fsca_0p_krew_arr<=limitFsca2)&(fsca_0p_krew_arr>=limitFsca1)])
    fsca_0p_jmz_lim = np.mean(fsca_0p_jmz_arr[(fsca_0p_jmz_arr<=limitFsca2)&(fsca_0p_jmz_arr>=limitFsca1)])
    fsca_0p_nrc_lim = np.mean(fsca_0p_nrc_arr[(fsca_0p_nrc_arr<=limitFsca2)&(fsca_0p_nrc_arr>=limitFsca1)])

    fsca_op_lim = np.array([fsca_0p_sc26m_lim,fsca_0p_sc17a_lim,fsca_0p_sc18m_lim,fsca_0p_krew_lim,fsca_0p_jmz_lim,fsca_0p_nrc_lim])
    
    return fsca_op_lim

def lin_reg(X,Y):
    barX=np.mean(X); barY=np.mean(Y)
    XminusbarX=X-barX; YminusbarY=Y-barY
    b1=sum(XminusbarX*YminusbarY)/sum(XminusbarX**2)
    b0=barY-b1*barX
    Yhat=b0+np.multiply(b1,X)
    e_i=Y-Yhat
    sse=np.sum(e_i**2)
    ssr=np.sum((Yhat-barY )**2)
    
    return {'Yhat':Yhat,'b0':b0,'b1':b1,'e_i': e_i,'sse':sse,'ssr':ssr}

def reshapeData(snow_temp_north_vegdens,col): 
    
    lat_long_nad83 = [np.min(snow_temp_north_vegdens[:,0]),np.max(snow_temp_north_vegdens[:,0]),np.min(snow_temp_north_vegdens[:,1]),np.max(snow_temp_north_vegdens[:,1])]
    ncols = (lat_long_nad83[1]-lat_long_nad83[0]+1).astype(int)
    nrows = (lat_long_nad83[3]-lat_long_nad83[2]+1).astype(int)
    snwMap = np.reshape(snow_temp_north_vegdens[:,col],[ncols,nrows])
    
    return snwMap

def readTiff_creatDF(tiffFile):    
    demset = gdal.Open(tiffFile)
    band = demset.GetRasterBand(1)
    lwr = band.ReadAsArray()
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = lwr.shape
    #x1 = x0 + dx * ncols
    #y1 = y0 + dy * nrows
    #extent1=[x0, x1, y1, y0]
    
    latitude =[]
    for x in range (ncols):
        latitude.append(x+x0)
    longitude = []
    for y in range (nrows):
        longitude.append(y0-y)
   
    latitude_rp = (np.tile(latitude, nrows))
    longitude_rp = (np.repeat(longitude, ncols))
    lwr_rp = np.reshape(lwr,(nrows*ncols)).T
    lwr_lat_lon1 = np.vstack([latitude_rp,longitude_rp,lwr_rp]).T

    lwr_df = pd.DataFrame(lwr_lat_lon1,columns=['x','y','lwr'])
    lwr_df.sort_values(by=['x','y'],inplace=True)
    lwr_df.index=np.arange(0,len(lwr_df))
    
    return lwr_df

def readTiff_creatDF_res100(tiffFile,columnName):    
    demset = gdal.Open(tiffFile)
    band = demset.GetRasterBand(1)
    lwr = band.ReadAsArray()
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = lwr.shape
    #x1 = x0 + dx * ncols
    #y1 = y0 + dy * nrows
    #extent1=[x0, x1, y1, y0]
    
    latitude =[]
    for x in range (ncols):
        latitude.append(100*x+x0)
    longitude = []
    for y in range (nrows):
        longitude.append(y0-100*y)
   
    latitude_rp = (np.tile(latitude, nrows))
    longitude_rp = (np.repeat(longitude, ncols))
    lwr_rp = np.reshape(lwr,(nrows*ncols)).T
    lwr_lat_lon1 = np.vstack([latitude_rp,longitude_rp,lwr_rp]).T

    lwr_df = pd.DataFrame(lwr_lat_lon1,columns=['x','y',columnName])
    lwr_df.sort_values(by=['x','y'],inplace=True)
    lwr_df.index=np.arange(0,len(lwr_df))
    
    return lwr_df

def changeResolution0fMultipleTiffsAndCreateDF (path2imagesIn,path2images0ut,columnName):
    # get the list of images
    listOfImgs=os.listdir(path2imagesIn)
    
    # load image using p2i
    fullPath2imgs = []
    for p2i in listOfImgs:
        fullPath2img=os.path.join(path2imagesIn,p2i)
        fullPath2imgs.append(fullPath2img)
        
    ##changing resolution
    for tif in range (len(fullPath2imgs)):
        gdal.Warp(path2images0ut[tif], fullPath2imgs[tif], xRes=100, yRes=100)
    
    df_100 = []
    for tif in range (len(path2images0ut)):
        wr_100 = readTiff_creatDF_res100(path2images0ut[0],columnName)
        df_100.append(wr_100) 
        
    all_wr = []
    for df in range (len(path2images0ut)):
        aaaa = df_100[df]['lwr'].values
        all_wr.append(aaaa)
    
    all_wr_mean = np.mean(all_wr, axis = 0)
    lat_swr = (df_100[0][['x']].values).T
    lon_swr = (df_100[0][['y']].values).T
    latLon_swr = np.vstack([lat_swr,lon_swr,all_wr_mean])
    
    return latLon_swr.T, fullPath2imgs

def creatingMeanRadiationfrom100b100tifFiles (path2tiff_100):
    lwr_df_krew = []    
    for lwr in range (len(path2tiff_100)):
        lwr_df_kr = readTiff_creatDF_res100(path2tiff_100[lwr],'lwr')
        lwr_df_krew.append(lwr_df_kr)
    
    latLon_kr = lwr_df_krew[0][['x','y']] 
    
    lwr_krew = []
    for df in range (len(lwr_df_krew)):
        aaaa = lwr_df_krew[df]['lwr'].values
        lwr_krew.append(aaaa)
    lwr_mean_kr = np.mean(lwr_krew, axis = 0)
    latLon_meanRadiation = np.vstack([latLon_kr['x'].values,latLon_kr['y'].values,lwr_mean_kr]).T
    
    return latLon_meanRadiation

def addClassifier2changeResolutionTo100B100 (lat_long_nad83_krew,snow_temp_vegdens_sL30_elevCls_precip_krew):
    x_krew = np.arange(lat_long_nad83_krew[0], lat_long_nad83_krew[1]+1, 100)
    y_krew = np.arange(lat_long_nad83_krew[2], lat_long_nad83_krew[3]+1, 100)
    xx_krew, yy_krew = np.meshgrid(x_krew, y_krew, sparse=True)

    yindx = np.zeros((len(snow_temp_vegdens_sL30_elevCls_precip_krew),1))
    snow_temp_vegdens_sL30_elevCls_precip_krew1 = np.append(snow_temp_vegdens_sL30_elevCls_precip_krew, yindx, 1)
    snow_temp_vegdens_sL30_elevCls_precip_krew2 = np.append(snow_temp_vegdens_sL30_elevCls_precip_krew1, yindx, 1)

    classifier_krew = np.arange(1000,len(x_krew)*len(y_krew)+1001)

    counter_x = 1
    for clss in range (len(x_krew)-1):
        snow_temp_vegdens_sL30_elevCls_precip_krew2[:,11][(snow_temp_vegdens_sL30_elevCls_precip_krew2[:,0]>=x_krew[clss])&
        (snow_temp_vegdens_sL30_elevCls_precip_krew2[:,0]<=x_krew[clss+1])] = counter_x
        counter_x += 1

    counter_y = 100
    for clss in range (len(y_krew)-1):
        snow_temp_vegdens_sL30_elevCls_precip_krew2[:,12][(snow_temp_vegdens_sL30_elevCls_precip_krew2[:,1]>=y_krew[clss])&
        (snow_temp_vegdens_sL30_elevCls_precip_krew2[:,1]<=y_krew[clss+1])] = counter_y
        counter_y += 100  

    classifier_xy = np.reshape((snow_temp_vegdens_sL30_elevCls_precip_krew2[:,11]+snow_temp_vegdens_sL30_elevCls_precip_krew2[:,12]), 
                               (len(snow_temp_vegdens_sL30_elevCls_precip_krew2),1))

    #adding classifier column for changing resolution                           
    snow_temp_vegdens_sL30_elevCls_precip_cls_krew = np.append(snow_temp_vegdens_sL30_elevCls_precip_krew2, classifier_xy, 1)

    return snow_temp_vegdens_sL30_elevCls_precip_cls_krew

def define_rf_features_fscaLatLonElevNorthTempPrecip_100B100 (classifier_xy,snow_temp_vegdens_sL30_elevCls_krew4):
    classifier_xy_n0duplicate = (pd.DataFrame(classifier_xy)).drop_duplicates(inplace=False)
    classifier_xy_n0duplicate.index = np.arange(0,len(classifier_xy_n0duplicate))

    #calculating fsca in 100 * 100 pixel
    fSCA_cls_krew = []
    for clsf in range (len(classifier_xy_n0duplicate)):#
         snw = len(snow_temp_vegdens_sL30_elevCls_krew4[(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsf])&((snow_temp_vegdens_sL30_elevCls_krew4[:,5]==77)|(snow_temp_vegdens_sL30_elevCls_krew4[:,5]==44))])
         n0snw = len(snow_temp_vegdens_sL30_elevCls_krew4[(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsf])&((snow_temp_vegdens_sL30_elevCls_krew4[:,5]==-7)|(snow_temp_vegdens_sL30_elevCls_krew4[:,5]==-4))])
         
         if snw+n0snw != 0:
             fSCA_ut = float(snw)/(snw+n0snw)
         else: fSCA_ut = 0
       
         fSCA_cls_krew.append(fSCA_ut)

    #calculating fsca under canopy and open in 100 * 100 pixel
    fSCA_cls_ut_krew = []
    fSCA_cls_0p_krew = []
    for cls in range (len(classifier_xy_n0duplicate)):#
         snwUT = len(snow_temp_vegdens_sL30_elevCls_krew4[(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][cls])&(snow_temp_vegdens_sL30_elevCls_krew4[:,5]==77)])
         n0snwUT = len(snow_temp_vegdens_sL30_elevCls_krew4[(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][cls])&(snow_temp_vegdens_sL30_elevCls_krew4[:,5]==-7)])
         snw0p = len(snow_temp_vegdens_sL30_elevCls_krew4[(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][cls])&(snow_temp_vegdens_sL30_elevCls_krew4[:,5]==44)])
         n0snw0p = len(snow_temp_vegdens_sL30_elevCls_krew4[(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][cls])&(snow_temp_vegdens_sL30_elevCls_krew4[:,5]==-4)])
         
         if snwUT+n0snwUT != 0:
             fSCA_ut = float(snwUT)/(snwUT+n0snwUT)
         else: fSCA_ut = 0
         
         if snw0p+n0snw0p != 0:
             fSCA_0p = float(snw0p)/(snw0p+n0snw0p) 
         else: fSCA_0p = 0
              
         fSCA_cls_ut_krew.append(fSCA_ut)
         fSCA_cls_0p_krew.append(fSCA_0p)
    
    #calculating average lat in 100 * 100 pixel
    lat_cls_krew = []
    for clslt in range (len(classifier_xy_n0duplicate)):#
         meanlat = np.min(snow_temp_vegdens_sL30_elevCls_krew4[:,0][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clslt])])
         lat_cls_krew.append(meanlat)
         
    #calculating average lon in 100 * 100 pixel
    lon_cls_krew = []
    for clsln in range (len(classifier_xy_n0duplicate)):#
         meanlon = np.min(snow_temp_vegdens_sL30_elevCls_krew4[:,1][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsln])])
         lon_cls_krew.append(meanlon)
    
    #calculating average elev in 100 * 100 pixel
    elev_cls_krew = []
    for clselv in range (len(classifier_xy_n0duplicate)):#
         meanlv = np.mean(snow_temp_vegdens_sL30_elevCls_krew4[:,2][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clselv])])
         elev_cls_krew.append(meanlv)
    
    #calculating average northness in 100 * 100 pixel
    nrth_cls_krew = []
    for clsN in range (len(classifier_xy_n0duplicate)):#
         meanNrth = np.mean(snow_temp_vegdens_sL30_elevCls_krew4[:,3][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsN])])
         nrth_cls_krew.append(meanNrth)
    
    #calculating average temp in 100 * 100 pixel
    temp_cls_krew = []
    for clsT in range (len(classifier_xy_n0duplicate)):#
         meanTemp = np.mean(snow_temp_vegdens_sL30_elevCls_krew4[:,7][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsT])])
         temp_cls_krew.append(meanTemp)
    
    #calculating average precip in 100 * 100 pixel
    precip_cls_krew = []
    for clsP in range (len(classifier_xy_n0duplicate)):#
         meanPrec = np.mean(snow_temp_vegdens_sL30_elevCls_krew4[:,10][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsP])])
         precip_cls_krew.append(meanPrec)
         
    #calculating average veg density in 100 * 100 pixel
    vegD_cls_krew = []
    for clsV in range (len(classifier_xy_n0duplicate)):#
         meanVD = np.mean(snow_temp_vegdens_sL30_elevCls_krew4[:,8][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsV])])
         vegD_cls_krew.append(meanVD)
    
    #calculating veg density under tree in 100 * 100 pixel
    vegD_ut_cls_krew = []
    for clsVu in range (len(classifier_xy_n0duplicate)):#
         meanVDut = np.mean(snow_temp_vegdens_sL30_elevCls_krew4[:,8][(snow_temp_vegdens_sL30_elevCls_krew4[:,13]==classifier_xy_n0duplicate[0][clsVu]) & ((snow_temp_vegdens_sL30_elevCls_krew4[:,5]>=70) | (snow_temp_vegdens_sL30_elevCls_krew4[:,5]<=-6))])
         vegD_ut_cls_krew.append(meanVDut)

    in_out_0p_rf_krew = pd.DataFrame(np.array([lat_cls_krew,lon_cls_krew,elev_cls_krew,nrth_cls_krew,
                                     temp_cls_krew,precip_cls_krew,vegD_cls_krew, vegD_ut_cls_krew, fSCA_cls_0p_krew,fSCA_cls_ut_krew,fSCA_cls_krew]).T, 
                                     columns = ['x','y','elev','nrth','temp','precip','VegD','vegD_ut','fsca_0p','fsca_ut','fsca'])

    return in_out_0p_rf_krew

def randomForest_MLM (in_0p_rf_krew, out_0p_rf_krew, random_state,n_trees):

    # split train and test; Kfold method
#    kf = KFold(n_splits=5) #spliting by kfold
#    index1, index2, index3, index4, index5= kf.split(in_0p_rf_krew)
#    
#    in_0p_rf_krew_train1 = in_0p_rf_krew.iloc[index1[0]] 
#    in_0p_rf_krew_test1 = in_0p_rf_krew.iloc[index1[1]] 
#    out_0p_rf_krew_train1 = out_0p_rf_krew.iloc[index1[0]] 
#    out_0p_rf_krew_test1 = out_0p_rf_krew.iloc[index1[1]]
    
    #another method to split data
    in_0p_rf_krew_train1, in_0p_rf_krew_test1, out_0p_rf_krew_train1, out_0p_rf_krew_test1 = train_test_split(in_0p_rf_krew, out_0p_rf_krew, test_size=0.3, random_state=random_state, shuffle= True) # 70% training and 30% test
                
    #Create a Gaussian Classifier
    clf1=RandomForestRegressor(n_estimators=n_trees)
    #Train the model using the training sets y_pred=clf.predict(X_test)
    clf1.fit(in_0p_rf_krew_train1,out_0p_rf_krew_train1)
    out_0p_rf_krew_pred1=clf1.predict(in_0p_rf_krew_test1)
    # Display the performance metrics
    #errors = np.mean(abs(y_pred - out_0p_rf_krew_test))
    errors1 = abs(out_0p_rf_krew_pred1 - out_0p_rf_krew_test1)
    mean_abs_error1 = round(np.mean(errors1), 3)
    #mape = np.mean(100 * (errors / out_0p_rf_krew_test))
    #accuracy = 100 - mape
    #print mean_abs_error1
    # Get numerical feature importances
    importances1 = list(clf1.feature_importances_)
    feature_list1 = in_0p_rf_krew.columns
    # List of tuples with variable and importance
    feature_importances1 = [(feature, round(importance, 3)) for feature, importance in zip(feature_list1, importances1)]
    # Sort the feature importances by most important first
    feature_importances2 = sorted(feature_importances1, key = lambda x: x[1], reverse = True)

    return mean_abs_error1,feature_importances1    

def randomForest_sensMod_predict (in_0p_rf_krew, out_0p_rf_krew, random_state, n_trees, input_sensitivit):
  
    in_0p_rf_krew_train1, in_0p_rf_krew_test1, out_0p_rf_krew_train1, out_0p_rf_krew_test1 = train_test_split(in_0p_rf_krew, out_0p_rf_krew, test_size=0.3, random_state=random_state, shuffle= True) # 70% training and 30% test
                
    #Create a Gaussian Classifier
    clf1=RandomForestRegressor(n_estimators=n_trees)
    #Train the model using the training sets y_pred=clf.predict(X_test)
    clf1.fit(in_0p_rf_krew_train1,out_0p_rf_krew_train1)

    out_0p_rf_krew_pred1=clf1.predict(in_0p_rf_krew_test1)

    errors1 = abs(out_0p_rf_krew_pred1 - out_0p_rf_krew_test1)
    mean_abs_error1 = round(np.mean(errors1), 3)
    # Get numerical feature importances
    importances1 = list(clf1.feature_importances_)
    feature_list1 = in_0p_rf_krew.columns
    # List of tuples with variable and importance
    feature_importances1 = [(feature, round(importance, 3)) for feature, importance in zip(feature_list1, importances1)]

    output_senst =clf1.predict(input_sensitivit)

    return [mean_abs_error1,feature_importances1,output_senst]   


def plot_rfm_importance (mean_abs_error1, importances1, pathName, facecolor,name):
    #ploting
    plt.subplots(figsize=(20,15))
    #plt.style.use('fivethirtyeight')
    # list of x locations for plotting
    x_values = list(range(len(importances1.columns)))
    # Make a bar chart
    plt.barh(x_values, importances1.values[0], facecolor = facecolor)#orientation = 'horizontal', 
    # Tick labels for x axis
    feature_list1 = importances1.columns
    plt.xticks(fontsize=35) #rotation='vertical', 
    plt.yticks(x_values, feature_list1, fontsize=35)
    #plt.invert_yaxis()
    # Axis labels and title
    plt.ylabel('Features', fontsize=40); plt.xlabel('Importance', fontsize=40) 
    plt.title('Variable Importances for {} (MAE={})'.format(name,mean_abs_error1), fontsize=35)
    plt.savefig('{}{}.png'.format(pathName,mean_abs_error1))

#rounding coodinates
def roundingCoordinates (in_out_rf_krew):
    for coords in range (len(in_out_rf_krew['x'])):
        if (in_out_rf_krew['x'][coords]/100. - (in_out_rf_krew['x'][coords]/100).astype(int)) != 0.0:
            adjust_num = (in_out_rf_krew['x'][coords]/100. - (in_out_rf_krew['x'][coords]/100).astype(int))
            new_x = in_out_rf_krew['x'][coords]-(adjust_num)*100.
            in_out_rf_krew ['x'][coords] = new_x .astype(int)  
        
        if (in_out_rf_krew['y'][coords]/100. - (in_out_rf_krew['y'][coords]/100).astype(int)) != 0.0:
            adjust_numY = (in_out_rf_krew['y'][coords]/100. - (in_out_rf_krew['y'][coords]/100).astype(int))
            new_y = in_out_rf_krew['y'][coords]-(adjust_numY)*100.
            in_out_rf_krew ['y'][coords] = new_y.astype(int)   
           

def randomForestModel44FscaBins(in_out_rf_krew_lswr_fcsa_drp_df):
    
    in_0pUt_rf_krew_b1 = in_out_rf_krew_lswr_fcsa_drp_df[['temp','precip','lwr','swr','VegD']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==1]
    out_0p_rf_krew_b1 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_0p']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==1])
    out_Ut_rf_krew_b1 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_ut']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==1])
    out_0pUt_rf_krew_b1 = pd.Series(out_0p_rf_krew_b1['fsca_0p'] - out_Ut_rf_krew_b1['fsca_ut'])
    
    in_0pUt_rf_krew_b2 = in_out_rf_krew_lswr_fcsa_drp_df[['temp','precip','lwr','swr','VegD']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==2]
    out_0p_rf_krew_b2 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_0p']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==2])
    out_Ut_rf_krew_b2 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_ut']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==2])
    out_0pUt_rf_krew_b2 = pd.Series(out_0p_rf_krew_b2['fsca_0p'] - out_Ut_rf_krew_b2['fsca_ut'])
    
    in_0pUt_rf_krew_b3 = in_out_rf_krew_lswr_fcsa_drp_df[['temp','precip','lwr','swr','VegD']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==3]
    out_0p_rf_krew_b3 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_0p']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==3])
    out_Ut_rf_krew_b3 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_ut']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==3])
    out_0pUt_rf_krew_b3 = pd.Series(out_0p_rf_krew_b3['fsca_0p'] - out_Ut_rf_krew_b3['fsca_ut'])
    
    in_0pUt_rf_krew_b4 = in_out_rf_krew_lswr_fcsa_drp_df[['temp','precip','lwr','swr','VegD']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==4]
    out_0p_rf_krew_b4 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_0p']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==4])
    out_Ut_rf_krew_b4 = (in_out_rf_krew_lswr_fcsa_drp_df[['fsca_ut']][in_out_rf_krew_lswr_fcsa_drp_df['fsca_classifier']==4])
    out_0pUt_rf_krew_b4 = pd.Series(out_0p_rf_krew_b4['fsca_0p'] - out_Ut_rf_krew_b4['fsca_ut'])
    
    #random forest for bins 1   
    mean_abs_error_0p11,feature_importances_0p11 = randomForest_MLM (in_0pUt_rf_krew_b1, out_0pUt_rf_krew_b1,60,200)
    mean_abs_error_0p12,feature_importances_0p12 = randomForest_MLM (in_0pUt_rf_krew_b1, out_0pUt_rf_krew_b1,60,200)
    mean_abs_error_0p13,feature_importances_0p13 = randomForest_MLM (in_0pUt_rf_krew_b1, out_0pUt_rf_krew_b1,60,200)
    mean_abs_error_0p14,feature_importances_0p14 = randomForest_MLM (in_0pUt_rf_krew_b1, out_0pUt_rf_krew_b1,60,200)
    mean_abs_error_0p15,feature_importances_0p15 = randomForest_MLM (in_0pUt_rf_krew_b1, out_0pUt_rf_krew_b1,60,200)
    
    feature_importances_0p1= [[feature_importances_0p11[i][1] for i in range (len(feature_importances_0p11))],
                             [feature_importances_0p12[i][1] for i in range (len(feature_importances_0p11))],
                             [feature_importances_0p13[i][1] for i in range (len(feature_importances_0p11))],
                             [feature_importances_0p14[i][1] for i in range (len(feature_importances_0p11))],
                             [feature_importances_0p15[i][1] for i in range (len(feature_importances_0p11))]]
    feature_importances_0p_mean1 = pd.DataFrame([np.mean([feature_importances_0p1 [i][0] for i in range (5)]),
                                   np.mean([feature_importances_0p1 [i][1] for i in range (5)]),
                                   np.mean([feature_importances_0p1 [i][2] for i in range (5)]),
                                   np.mean([feature_importances_0p1 [i][3] for i in range (5)]),
                                   np.mean([feature_importances_0p1 [i][4] for i in range (5)])]).T
    
    feature_importances_0p_mean1.columns = in_0pUt_rf_krew_b1.columns
    mean_abs_error_0p1 = (mean_abs_error_0p11+mean_abs_error_0p12+mean_abs_error_0p13+mean_abs_error_0p14+mean_abs_error_0p15)/5.

    #random forest for bins 2   
    mean_abs_error_0p21,feature_importances_0p21 = randomForest_MLM (in_0pUt_rf_krew_b2, out_0pUt_rf_krew_b2,60,200)
    mean_abs_error_0p22,feature_importances_0p22 = randomForest_MLM (in_0pUt_rf_krew_b2, out_0pUt_rf_krew_b2,60,200)
    mean_abs_error_0p23,feature_importances_0p23 = randomForest_MLM (in_0pUt_rf_krew_b2, out_0pUt_rf_krew_b2,60,200)
    mean_abs_error_0p24,feature_importances_0p24 = randomForest_MLM (in_0pUt_rf_krew_b2, out_0pUt_rf_krew_b2,60,200)
    mean_abs_error_0p25,feature_importances_0p25 = randomForest_MLM (in_0pUt_rf_krew_b2, out_0pUt_rf_krew_b2,60,200)
    
    feature_importances_0p2= [[feature_importances_0p21[i][1] for i in range (len(feature_importances_0p21))],
                             [feature_importances_0p22[i][1] for i in range (len(feature_importances_0p21))],
                             [feature_importances_0p23[i][1] for i in range (len(feature_importances_0p21))],
                             [feature_importances_0p24[i][1] for i in range (len(feature_importances_0p21))],
                             [feature_importances_0p25[i][1] for i in range (len(feature_importances_0p21))]]
    feature_importances_0p_mean2 = pd.DataFrame([np.mean([feature_importances_0p2 [i][0] for i in range (5)]),
                                   np.mean([feature_importances_0p2 [i][1] for i in range (5)]),
                                   np.mean([feature_importances_0p2 [i][2] for i in range (5)]),
                                   np.mean([feature_importances_0p2 [i][3] for i in range (5)]),
                                   np.mean([feature_importances_0p2 [i][4] for i in range (5)])]).T
    
    feature_importances_0p_mean2.columns = in_0pUt_rf_krew_b2.columns
    mean_abs_error_0p2 = (mean_abs_error_0p21+mean_abs_error_0p22+mean_abs_error_0p23+mean_abs_error_0p24+mean_abs_error_0p25)/5.

    #random forest for bins 3   
    mean_abs_error_0p31,feature_importances_0p31 = randomForest_MLM (in_0pUt_rf_krew_b3, out_0pUt_rf_krew_b3,50,200)
    mean_abs_error_0p32,feature_importances_0p32 = randomForest_MLM (in_0pUt_rf_krew_b3, out_0pUt_rf_krew_b3,50,200)
    mean_abs_error_0p33,feature_importances_0p33 = randomForest_MLM (in_0pUt_rf_krew_b3, out_0pUt_rf_krew_b3,50,200)
    mean_abs_error_0p34,feature_importances_0p34 = randomForest_MLM (in_0pUt_rf_krew_b3, out_0pUt_rf_krew_b3,50,200)
    mean_abs_error_0p35,feature_importances_0p35 = randomForest_MLM (in_0pUt_rf_krew_b3, out_0pUt_rf_krew_b3,50,200)
    
    feature_importances_0p3= [[feature_importances_0p31[i][1] for i in range (len(feature_importances_0p31))],
                             [feature_importances_0p32[i][1] for i in range (len(feature_importances_0p31))],
                             [feature_importances_0p33[i][1] for i in range (len(feature_importances_0p31))],
                             [feature_importances_0p34[i][1] for i in range (len(feature_importances_0p31))],
                             [feature_importances_0p35[i][1] for i in range (len(feature_importances_0p31))]]
    feature_importances_0p_mean3 = pd.DataFrame([np.mean([feature_importances_0p3 [i][0] for i in range (5)]),
                                   np.mean([feature_importances_0p3 [i][1] for i in range (5)]),
                                   np.mean([feature_importances_0p3 [i][2] for i in range (5)]),
                                   np.mean([feature_importances_0p3 [i][3] for i in range (5)]),
                                   np.mean([feature_importances_0p3 [i][4] for i in range (5)])]).T
    
    feature_importances_0p_mean3.columns = in_0pUt_rf_krew_b3.columns
    mean_abs_error_0p3 = (mean_abs_error_0p31+mean_abs_error_0p32+mean_abs_error_0p33+mean_abs_error_0p34+mean_abs_error_0p35)/5.

    #random forest for bins 4   
    mean_abs_error_0p41,feature_importances_0p41 = randomForest_MLM (in_0pUt_rf_krew_b4, out_0pUt_rf_krew_b4,50,200)
    mean_abs_error_0p42,feature_importances_0p42 = randomForest_MLM (in_0pUt_rf_krew_b4, out_0pUt_rf_krew_b4,50,200)
    mean_abs_error_0p43,feature_importances_0p43 = randomForest_MLM (in_0pUt_rf_krew_b4, out_0pUt_rf_krew_b4,50,200)
    mean_abs_error_0p44,feature_importances_0p44 = randomForest_MLM (in_0pUt_rf_krew_b4, out_0pUt_rf_krew_b4,50,200)
    mean_abs_error_0p45,feature_importances_0p45 = randomForest_MLM (in_0pUt_rf_krew_b4, out_0pUt_rf_krew_b4,50,200)
    
    feature_importances_0p4= [[feature_importances_0p41[i][1] for i in range (len(feature_importances_0p41))],
                             [feature_importances_0p42[i][1] for i in range (len(feature_importances_0p41))],
                             [feature_importances_0p43[i][1] for i in range (len(feature_importances_0p41))],
                             [feature_importances_0p44[i][1] for i in range (len(feature_importances_0p41))],
                             [feature_importances_0p45[i][1] for i in range (len(feature_importances_0p41))]]
    feature_importances_0p_mean4 = pd.DataFrame([np.mean([feature_importances_0p4 [i][0] for i in range (5)]),
                                   np.mean([feature_importances_0p4 [i][1] for i in range (5)]),
                                   np.mean([feature_importances_0p4 [i][2] for i in range (5)]),
                                   np.mean([feature_importances_0p4 [i][3] for i in range (5)]),
                                   np.mean([feature_importances_0p4 [i][4] for i in range (5)])]).T
    
    feature_importances_0p_mean4.columns = in_0pUt_rf_krew_b4.columns
    mean_abs_error_0p4 = (mean_abs_error_0p41+mean_abs_error_0p42+mean_abs_error_0p43+mean_abs_error_0p44+mean_abs_error_0p45)/5.

    return [[feature_importances_0p_mean1,mean_abs_error_0p1],[feature_importances_0p_mean2,mean_abs_error_0p2],[feature_importances_0p_mean3,mean_abs_error_0p3],[feature_importances_0p_mean4,mean_abs_error_0p4]]


def randomForestModel4Fsca0pen (in_0p_rf_krew,out_0p_rf_krew,random_state,n_trees):
        
    mean_abs_error_0p1,feature_importances_0p1 = randomForest_MLM (in_0p_rf_krew, out_0p_rf_krew,random_state,n_trees)
    mean_abs_error_0p2,feature_importances_0p2 = randomForest_MLM (in_0p_rf_krew, out_0p_rf_krew,random_state,n_trees)
    mean_abs_error_0p3,feature_importances_0p3 = randomForest_MLM (in_0p_rf_krew, out_0p_rf_krew,random_state,n_trees)
    mean_abs_error_0p4,feature_importances_0p4 = randomForest_MLM (in_0p_rf_krew, out_0p_rf_krew,random_state,n_trees)
    mean_abs_error_0p5,feature_importances_0p5 = randomForest_MLM (in_0p_rf_krew, out_0p_rf_krew,random_state,n_trees)
    
    feature_importances_0p= [[feature_importances_0p1[i][1] for i in range (len(feature_importances_0p1))],
                             [feature_importances_0p2[i][1] for i in range (len(feature_importances_0p1))],
                             [feature_importances_0p3[i][1] for i in range (len(feature_importances_0p1))],
                             [feature_importances_0p4[i][1] for i in range (len(feature_importances_0p1))],
                             [feature_importances_0p5[i][1] for i in range (len(feature_importances_0p1))]]
    feature_importances_0p_mean = pd.DataFrame([np.mean([feature_importances_0p [i][0] for i in range (5)]),
                                  np.mean([feature_importances_0p [i][1] for i in range (5)]),
                                  np.mean([feature_importances_0p [i][2] for i in range (5)]),
                                  np.mean([feature_importances_0p [i][3] for i in range (5)])]).T
    feature_importances_0p_mean.columns = in_0p_rf_krew.columns
    mean_abs_error_0p = (mean_abs_error_0p1+mean_abs_error_0p2+mean_abs_error_0p3+mean_abs_error_0p4+mean_abs_error_0p5)/5.
    
    return [feature_importances_0p_mean,mean_abs_error_0p]

def randomForestModel4FscaUnderTree0r40p_ut(in_ut_rf_krew,out_ut_rf_krew,random_state,n_trees):
    
    mean_abs_error_ut1,feature_importances_ut1 = randomForest_MLM (in_ut_rf_krew, out_ut_rf_krew,random_state,n_trees)
    mean_abs_error_ut2,feature_importances_ut2 = randomForest_MLM (in_ut_rf_krew, out_ut_rf_krew,random_state,n_trees)
    mean_abs_error_ut3,feature_importances_ut3 = randomForest_MLM (in_ut_rf_krew, out_ut_rf_krew,random_state,n_trees)
    mean_abs_error_ut4,feature_importances_ut4 = randomForest_MLM (in_ut_rf_krew, out_ut_rf_krew,random_state,n_trees)
    mean_abs_error_ut5,feature_importances_ut5 = randomForest_MLM (in_ut_rf_krew, out_ut_rf_krew,random_state,n_trees)
    
    feature_importances_ut= [[feature_importances_ut1[i][1] for i in range (len(feature_importances_ut1))],
                             [feature_importances_ut2[i][1] for i in range (len(feature_importances_ut1))],
                             [feature_importances_ut3[i][1] for i in range (len(feature_importances_ut1))],
                             [feature_importances_ut4[i][1] for i in range (len(feature_importances_ut1))],
                             [feature_importances_ut5[i][1] for i in range (len(feature_importances_ut1))]]
    feature_importances_ut_mean = pd.DataFrame([np.mean([feature_importances_ut [i][0] for i in range (5)]),
                                  np.mean([feature_importances_ut [i][1] for i in range (5)]),
                                  np.mean([feature_importances_ut [i][2] for i in range (5)]),
                                  np.mean([feature_importances_ut [i][3] for i in range (5)]),
                                  np.mean([feature_importances_ut [i][4] for i in range (5)])]).T
    feature_importances_ut_mean.columns = in_ut_rf_krew.columns
    mean_abs_error_ut = (mean_abs_error_ut1+mean_abs_error_ut2+mean_abs_error_ut3+mean_abs_error_ut4+mean_abs_error_ut5)/5.

    return [feature_importances_ut_mean,mean_abs_error_ut]

def fscaPlottingByColorMap4tempNrth(in_out_rf_krew_lswr_ls_drp_df,out_0p_rf_krew,cmap,cbarLabel,savepath):
    plt.subplots(figsize=(20,15))
    #norm= matplotlib.colors.Normalize(vmin=0,vmax=1)

    plt.scatter(in_out_rf_krew_lswr_ls_drp_df['temp'], in_out_rf_krew_lswr_ls_drp_df['nrth'], 
                c=out_0p_rf_krew, cmap=cmap, s = 20**2, linewidth=0.5)#spcb = 
    
    vmin = np.min(out_0p_rf_krew)
    vmax = np.max(out_0p_rf_krew)
    mid = (vmin+vmax)/2

    cbar= plt.colorbar(ticks=[vmin, mid, vmax])#spcb,-0.6, -0.2, 0.2, 0.6
    cbar.set_label(cbarLabel, fontsize=25, labelpad=+1)#
    cbar.ax.set_yticklabels([format(vmin, '.2f'), format(mid, '.2f'), format(vmax, '.2f')],fontsize=25) #,'0.2','0.4','0.6','0.8'
    
    plt.xticks(fontsize=35) #rotation='vertical', 
    plt.yticks(fontsize=35)
    plt.ylabel('northness', fontsize=40); plt.xlabel('temperature (C)', fontsize=40) 
    plt.title('Change in {} based on temp&nrth'.format(cbarLabel), fontsize=40)
    plt.savefig(savepath)
    
def fscaPlottingByColorMap4swrLwr(in_out_rf_krew_lswr_ls_drp_df,out_0p_rf_krew,cmap,cbarLabel,savepath):
    plt.subplots(figsize=(20,15))
    plt.scatter(in_out_rf_krew_lswr_ls_drp_df['swr'], in_out_rf_krew_lswr_ls_drp_df['lwr'], 
                c=out_0p_rf_krew, cmap=cmap, s = 20**2, linewidth=0.5)
    
    vmin = np.min(out_0p_rf_krew)
    vmax = np.max(out_0p_rf_krew)
    mid = (vmin+vmax)/2

    cbar= plt.colorbar(ticks=[vmin, mid, vmax])#spcb,-0.6, -0.2, 0.2, 0.6
    cbar.set_label(cbarLabel, fontsize=25, labelpad=+1)#
    cbar.ax.set_yticklabels([format(vmin, '.2f'), format(mid, '.2f'), format(vmax, '.2f')],fontsize=25) #,'0.2','0.4','0.6','0.8'
    
    plt.xticks(fontsize=35) #rotation='vertical', 
    plt.yticks(fontsize=35)
    plt.ylabel('lwr', fontsize=40); plt.xlabel('swr', fontsize=40) 
    plt.title('change in {} based on swr & lwr'.format(cbarLabel), fontsize=40)
    plt.savefig(savepath)

def plotrandomForestimportancesfor4bins(feature_importances_error_0putBin,importance_bins_krew,color_krew,siteName,savePath):   
    plt.subplots(figsize=(20,15))
    y_values = np.arange(0,len(feature_importances_error_0putBin[0][0].columns))
    indx = -1.6
    bar_width = 0.2
    for bins in range(len(importance_bins_krew)):
        plt.barh(y_values+indx*bar_width, importance_bins_krew[bins], bar_width, facecolor = color_krew[bins])#, alpha=0.2,orientation = 'horizontal', 
        indx += 1
        
    # Tick labels for x axis
    feature_list1 = feature_importances_error_0putBin[0][0].columns
    plt.xticks(fontsize=35) #rotation='vertical', 
    plt.yticks(y_values, feature_list1, fontsize=35)
    
    plt.ylabel('Features', fontsize=40); plt.xlabel('Importance', fontsize=40) 
    
    meanAbsError = [feature_importances_error_0putBin[i][1] for i in range (4)]
    fSCA_bin=['<0.3','0.3-0.55','0.55-0.8','>0.8']
    legend = ['fSCA_bin {}; MAE = {}'.format(fSCA_bin[i],meanAbsError[i]) for i in range (4)]
    plt.legend(legend,fontsize = 22)#, loc = 'upper center'
    plt.title('Importance of features in 4 fSCA bins in {}'.format(siteName), fontsize=40)#, position=(0.04, 0.88), ha='left'
    
    plt.savefig(savePath)
randomForest_sensMod_predict
def plotSensitivitySWRvegDensity(sens_vegSwr_krew_up,sens_vegSwr_krew_low,lent1,lent2,savePath1,color_up,color_low,siteName):

    vegDens_sens = np.arange(0.01*lent1,lent2*0.01,0.01)

    plt.subplots(figsize=(30,20))
    plt.plot(vegDens_sens,sens_vegSwr_krew_up[2], c = color_up, mec = 'black', mfc = color_up, 
             marker = '^', linewidth=7, markersize=20)
    plt.plot(vegDens_sens,sens_vegSwr_krew_low[2],c = color_low, mec = 'black', mfc = color_low, 
             marker = 'v', linewidth=7, markersize=20)
    plt.xticks(fontsize=35) #rotation='vertical', 
    plt.yticks(fontsize=35)
    plt.ylabel('Predicted fSCAopen-fSCAunder in {}'.format(siteName), fontsize=40)
    plt.xlabel('Vegetation density', fontsize=40) 
    legend = ['Predicted fSCAopen-fSCAunder for 80% of SWR in {}'.format(siteName),'Predicted fSCAopen-fSCAunder for 20% of SWR in {}'.format(siteName)]
    plt.legend(legend,fontsize = 34)#, loc = 'upper center'
    plt.title('Predicted fSCAopen-fSCAunder for 2 percentiles of SWR in {}'.format(siteName), fontsize=40, y = 1.04)#, position=(0.04, 0.88), ha='left'
    plt.savefig(savePath1)


#def plotSensitivitySWRvegDensity(sens_vegSwr_krew_up,sens_vegSwr_krew_low,savePath1,savePath2,color_up,color_low,siteName,cmap):
#
#    vegDens_sens = np.arange(0.01,1.01,0.01)
#
#    plt.subplots(figsize=(30,20))
#    plt.plot(vegDens_sens,sens_vegSwr_krew_up[2], c = color_up, mec = 'black', mfc = color_up, 
#             marker = '^', linewidth=7, markersize=20)
#    plt.plot(vegDens_sens,sens_vegSwr_krew_low[2],c = color_low, mec = 'black', mfc = color_low, 
#             marker = 'v', linewidth=7, markersize=20)
#    plt.xticks(fontsize=35) #rotation='vertical', 
#    plt.yticks(fontsize=35)
#    plt.ylabel('Predicted fSCAopen-fSCAunder in {}'.format(siteName), fontsize=40)
#    plt.xlabel('Vegetation density', fontsize=40) 
#    legend = ['Predicted fSCAopen-fSCAunder for 80% of SWR in {}'.format(siteName),'Predicted fSCAopen-fSCAunder for 20% of SWR in {}'.format(siteName)]
#    plt.legend(legend,fontsize = 34)#, loc = 'upper center'
#    plt.title('Predicted fSCAopen-fSCAunder for 2 percentiles of SWR in {}'.format(siteName), fontsize=40, y = 1.04)#, position=(0.04, 0.88), ha='left'
#    plt.savefig(savePath1)
                
    
#    fig, ax = plt.subplots(figsize=(30,25))
#    plt.scatter(sens_vegSwr_krew_low[2], sens_vegSwr_krew_up[2], 
#                c=vegDens_sens, cmap=cmap, s = 27**2, linewidth=0.5)#spcb = 
#    vmin = 0.01
#    vmax = 1
#    mid = (vmin+vmax)/2
#    cbar= plt.colorbar(ticks=[vmin, mid, vmax])#spcb,-0.6, -0.2, 0.2, 0.6
#    cbar.set_label('Vegetation density (VD)', fontsize=30, labelpad=+1)#
#    cbar.ax.set_yticklabels([format(vmin, '.2f'), format(mid, '.2f'), format(vmax, '.2f')],fontsize=25) #,'0.2','0.4','0.6','0.8'
#    plt.xticks(np.arange(round(min(min(sens_vegSwr_krew_low[2]),min(sens_vegSwr_krew_up[2]))-0.02,2),
#                round(max(max(sens_vegSwr_krew_low[2]),max(sens_vegSwr_krew_up[2]))+0.02,2),0.04),fontsize=40) #rotation='vertical', 
#    plt.yticks(np.arange(round(min(min(sens_vegSwr_krew_low[2]),min(sens_vegSwr_krew_up[2]))-0.02,2),
#                round(max(max(sens_vegSwr_krew_low[2]),max(sens_vegSwr_krew_up[2]))+0.02,2),0.04),fontsize=40)
#    plt.ylabel('Predicted fSCAopen-fSCAunder for 80% of SWR in {}'.format(siteName), fontsize=40)
#    plt.xlabel('Predicted fSCAopen-fSCAunder for 20% of SWR in {}'.format(siteName), fontsize=40) 
#    ax.xaxis.set_label_coords(0.5, -0.07)
#    plt.title('Predicted fSCAopen-fSCAunder for 80%  vs. 20% of SWR for different VD in {}'.format(siteName), 
#              fontsize=32, y =1.03)#, position=(0.04, 0.88), ha='left'
#    
#    plt.savefig(savePath2)

        
# ploting features with lines   
#plt.subplots(figsize=(20,15))
#for countFeat in range (len(feature_importances_error_0putBin)): 
#    plt.plot(feature_importances_error_0putBin[countFeat][0].values[0], marker='s', markersize=20, linewidth=3.5)
#
#meanAbsError = [feature_importances_error_0putBin[i][1] for i in range (4)]
#legend = ['bin{}; MAE = {}'.format(i+1,meanAbsError[i]) for i in range (4)]
#plt.legend(legend,fontsize = 40)#, loc = 'upper center'
#plt.title('Importance of features in 4 fscaBins in Krew', fontsize=40)#, position=(0.04, 0.88), ha='left'
#
#ax = np.reshape(np.arange(0,len(feature_importances_error_0putBin[0][0].columns)),(1,5))
#sa_xticks = feature_importances_error_0putBin[0][0].columns
#plt.xticks(ax[0], sa_xticks, fontsize=40)#, rotation=25 
#plt.yticks(fontsize=40)
#plt.xlabel('Features', fontsize=40)
#plt.ylabel('Importance', fontsize=40)   
    
    
    
    
    
    