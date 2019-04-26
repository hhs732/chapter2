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

#%% load data nrc 2010
#calculating topo_dimension from cut northness to plot snow on
snow_temp_vegdens_allExtent_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_nrc.npy')
ouputPath_nrc = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_contour_nrc.png'
T_nrc = 'NRC'
plotSnowMap(snow_temp_vegdens_allExtent_nrc,ouputPath_nrc,T_nrc)

ouputPath_nrt_nrc = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_nrth_contour_nrc.png'
plotSnowMap_northness(snow_temp_vegdens_allExtent_nrc,ouputPath_nrt_nrc,T_nrc)

#ascii_grid_veg_indexnrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_nrc.npy')
#veg_density_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/veg_density_bl.npy')

#lat_long_nad83_nrc = [np.min(snow_temp_vegdens_allExtent_nrc[:,0]),np.max(snow_temp_vegdens_allExtent_nrc[:,0]),np.min(snow_temp_vegdens_allExtent_nrc[:,1]),np.max(snow_temp_vegdens_allExtent_nrc[:,1])]
#lat_long_degre_nrc = [39.978183, 40.025319, -105.58556840332032, -105.52620245393197]

##y = -0.007512x + 14.168475  # R² = 0.969091
#elev_nrc = ascii_grid_veg_indexnrc[:,2]
#tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475 
#
#veg_snow_temp_density_nrc = np.vstack([ascii_grid_veg_indexnrc[:,0],ascii_grid_veg_indexnrc[:,1].astype(int),ascii_grid_veg_indexnrc[:,2],
#                                       ascii_grid_veg_indexnrc[:,3],ascii_grid_veg_indexnrc[:,4],ascii_grid_veg_indexnrc[:,5],
#                                       ascii_grid_veg_indexnrc[:,6],tempDJF_gm_nrc,veg_density_nrc[:,2],tempDJF_gm_nrc,tempDJF_gm_nrc]).T
# slope less than 30
#snow_temp_vegdens_index_sL30_nrc = veg_snow_temp_density_nrc[(veg_snow_temp_density_nrc[:,5]>-90)&(veg_snow_temp_density_nrc[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_nrc',veg_snow_temp_density_nrc)

snow_temp_vegdens_index_sL30_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_nrc.npy')
aaaa_test1 = snow_temp_vegdens_index_sL30_nrc[0:19710,:]

elev_nrc = snow_temp_vegdens_index_sL30_nrc[:,2]
tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475

snow_temp_vegdens_index_sL30_nrc_df = pd.DataFrame(snow_temp_vegdens_index_sL30_nrc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
elevationBand_nrc = np.arange(min(snow_temp_vegdens_index_sL30_nrc_df['z']),max(snow_temp_vegdens_index_sL30_nrc_df['z']),65)
snow_temp_vegdens_sL30_elevCls_nrc = classificationElevation(snow_temp_vegdens_index_sL30_nrc,elevationBand_nrc,elev_class_inx)
aaaa_test2 = snow_temp_vegdens_sL30_elevCls_nrc[0:710,:]

meanTemp_elevClass_nrc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)   

#sampling
sample_size = 100
fsca_ut_nrc, fsca_0p_nrc = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)
fsca_ut_samp_nrc, fsca_0p_samp_nrc = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx,sample_size)
statistic_nrc1,pvalue_nrc1,average_sample_exp_nrc1,average_sample_shl_nrc1,significant_nrc1 = wilcoxonTest(fsca_ut_samp_nrc, fsca_0p_samp_nrc)

# under canopy and 0pen fsca // exposed vs. sheltered
fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_nrc, fSCA_0u_shl_nrc = fsca0pn_fscaUnTr(fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc, elev_class_inx)
#sampling
fSCA_0u_exp_nrc_samp, fSCA_0u_shl_nrc_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx,sample_size)
statistic_nrc2,pvalue_nrc2,average_sample_exp_nrc2,average_sample_shl_nrc2,significant_nrc2 = wilcoxonTest(fSCA_0u_exp_nrc_samp, fSCA_0u_shl_nrc_samp)

# December
#y = -0.00777x + 17.81629
tempDec_nrc = -0.00777*elevationBand_nrc + 17.81629

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_nrc = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)
fSCA_ut_exp_l_sam_nrc, fSCA_ut_exp_h_sam_nrc, fSCA_ut_shl_l_sam_nrc, fSCA_ut_shl_h_sam_nrc = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx,sample_size)
statistic_ex_nrc,pvalue_ex_nrc,average_sample_l_ex_nrc,average_sample_h_ex_nrc,significant_ex_nrc = wilcoxonTest(fSCA_ut_exp_l_sam_nrc, fSCA_ut_exp_h_sam_nrc)
statistic_sh_nrc,pvalue_sh_nrc,average_sample_l_sh_nrc,average_sample_h_sh_nrc,significant_sh_nrc = wilcoxonTest(fSCA_ut_shl_l_sam_nrc, fSCA_ut_shl_h_sam_nrc)
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/significant_nrc',[significant_nrc1,significant_nrc2,significant_ex_nrc,significant_sh_nrc])

#%% load data Jemez 2010
#ascii_grid_veg_indexJmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_jmz.npy')
#veg_density_jmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/veg_density_jm.npy')
#
#elev_jmz = ascii_grid_veg_indexJmz[:,2]
#tempDJF_gm_jmz = -0.0046*elev_jmz + 7.7058 
#
#veg_snow_temp_density_jmz = np.vstack([ascii_grid_veg_indexJmz[:,0],ascii_grid_veg_indexJmz[:,1].astype(int),ascii_grid_veg_indexJmz[:,2],
#                                       ascii_grid_veg_indexJmz[:,3],ascii_grid_veg_indexJmz[:,4],ascii_grid_veg_indexJmz[:,5],
#                                       ascii_grid_veg_indexJmz[:,6],tempDJF_gm_jmz,veg_density_jmz[:,2],tempDJF_gm_jmz,tempDJF_gm_jmz]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_jmz',veg_snow_temp_density_jmz)

## slope less than 30
#snow_temp_vegdens_index_sL30_jmz = veg_snow_temp_density_jmz[(veg_snow_temp_density_jmz[:,5]>-50)&(veg_snow_temp_density_jmz[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_jmz',snow_temp_vegdens_index_sL30_jmz)

snow_temp_vegdens_index_sL30_jmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_jmz.npy')
aaaa_test1 = snow_temp_vegdens_index_sL30_jmz[0:710,:]

elev_jmz = snow_temp_vegdens_index_sL30_jmz[:,2]
tempDJF_gm_jmz = -0.0046*elev_jmz + 7.7058 

#calculating topo_dimension from cut northness
lat_long_nad83_Jmz = [np.min(snow_temp_vegdens_index_sL30_jmz[:,0]),np.max(snow_temp_vegdens_index_sL30_jmz[:,0]),np.min(snow_temp_vegdens_index_sL30_jmz[:,1]),np.max(snow_temp_vegdens_index_sL30_jmz[:,1])]
lat_long_degre_Jmz = [35.863504, 35.918451, -106.60600014573875, -106.54061446340079]

snow_temp_vegdens_index_sL30_jmz_df = pd.DataFrame(snow_temp_vegdens_index_sL30_jmz, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_jmz = np.arange(min(snow_temp_vegdens_index_sL30_jmz_df['z']),max(snow_temp_vegdens_index_sL30_jmz_df['z']),90)
snow_temp_vegdens_sL30_elevCls_jmz = classificationElevation(snow_temp_vegdens_index_sL30_jmz,elevationBand_jmz,elev_class_inx)

meanTemp_elevClass_jmz = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_jmz, fsca_0p_jmz = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)
fsca_ut_samp_jmz, fsca_0p_samp_jmz = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx,sample_size)
#statistical test
statistic_jmz1,pvalue_jmz1,average_sample_exp_jmz1,average_sample_shl_jmz1,significant_jmz1 = wilcoxonTest(fsca_ut_samp_jmz, fsca_0p_samp_jmz)

# under canopy and 0pen fsca exposed vs. open
fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz = fsca0pn_fscaUnTr(fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz, elev_class_inx)
#statistical test
fSCA_0u_exp_Jmz_samp, fSCA_0u_shl_Jmz_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx,sample_size)
statistic_jmz2,pvalue_jmz2,average_sample_exp_jmz2,average_sample_shl_jmz2,significant_jmz2 = wilcoxonTest(fSCA_0u_exp_Jmz_samp, fSCA_0u_shl_Jmz_samp)

# y = -0.0043x + 5.6707
tempDec_jmz = -0.0048*elevationBand_jmz + 8.6707

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_jmz = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)
fSCA_ut_exp_l_sam_jmz, fSCA_ut_exp_h_sam_jmz, fSCA_ut_shl_l_sam_jmz, fSCA_ut_shl_h_sam_jmz = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx,sample_size)
statistic_ex_jmz,pvalue_ex_jmz,average_sample_l_ex_jmz,average_sample_h_ex_jmz,significant_ex_jmz = wilcoxonTest(fSCA_ut_exp_l_sam_jmz, fSCA_ut_exp_h_sam_jmz)
statistic_sh_jmz,pvalue_sh_jmz,average_sample_l_sh_jmz,average_sample_h_sh_jmz,significant_sh_jmz = wilcoxonTest(fSCA_ut_shl_l_sam_jmz, fSCA_ut_shl_h_sam_jmz)
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/significant_jmz',[significant_jmz1,significant_jmz2,significant_ex_jmz,significant_sh_jmz])

#calculating topo_dimension from cut northness to plot snow on
snow_temp_vegdens_allExtent_jmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_jmz.npy')
ouputPath_jmz = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_contour_jmz.png'
T_jmz = 'JRB'
plotSnowMap(snow_temp_vegdens_allExtent_jmz,ouputPath_jmz,T_jmz)

lat_long_nad83_jmz = [np.min(snow_temp_vegdens_allExtent_jmz[:,0]),np.max(snow_temp_vegdens_allExtent_jmz[:,0]),np.min(snow_temp_vegdens_allExtent_jmz[:,1]),np.max(snow_temp_vegdens_allExtent_jmz[:,1])]

ouputPath_nrt_jmz = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_nrth_contour_jmz.png'
plotSnowMap_northness(snow_temp_vegdens_allExtent_jmz,ouputPath_nrt_jmz,T_jmz)
#%% load data sagehen creek 26 March 2016
ascii_grid_veg_index26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc26m.npy')
#calculating topo_dimension from cut northness
veg_density_sc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/veg_density_sc.npy')

elev_sc26m = ascii_grid_veg_index26m[:,2]
tempDJF_sc26m = -0.0016 * elev_sc26m + 1.802

veg_snow_temp_density_sc26m = np.vstack([ascii_grid_veg_index26m[:,0],ascii_grid_veg_index26m[:,1].astype(int),ascii_grid_veg_index26m[:,2],
                                         ascii_grid_veg_index26m[:,3],ascii_grid_veg_index26m[:,4],ascii_grid_veg_index26m[:,5],
                                         ascii_grid_veg_index26m[:,6],tempDJF_sc26m,veg_density_sc[:,2],tempDJF_sc26m,tempDJF_sc26m]).T
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sc26m',veg_snow_temp_density_sc26m)
# slope less than 30
snow_temp_vegdens_index_sL30_sc26m = veg_snow_temp_density_sc26m[(veg_snow_temp_density_sc26m[:,5]>-50)&(veg_snow_temp_density_sc26m[:,4]<30)]
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_sc26m',snow_temp_vegdens_index_sL30_sc26m)

snow_temp_vegdens_index_sL30_sc26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_sc26m.npy')
aaaa_test1 = snow_temp_vegdens_index_sL30_sc26m[0:710,:]

elev_sc26m = snow_temp_vegdens_index_sL30_sc26m[:,2]
tempDJF_sc26m = -0.0016 * elev_sc26m + 1.802

snow_temp_vegdens_index_sL30_sc26m_df = pd.DataFrame(snow_temp_vegdens_index_sL30_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_sc = np.arange(min(snow_temp_vegdens_index_sL30_sc26m_df['z']),max(snow_temp_vegdens_index_sL30_sc26m_df['z']),67)
snow_temp_vegdens_sL30_elevCls_sc26m = classificationElevation(snow_temp_vegdens_index_sL30_sc26m,elevationBand_sc,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc26m, fsca_0p_sc26m = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)
fsca_ut_samp_sc26m, fsca_0p_samp_sc26m = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx,sample_size)
#statistical test
statistic_sc26m1,pvalue_sc26m1,average_sample_exp_sc26m1,average_sample_shl_sc26m1,significant_sc26m1 = wilcoxonTest(fsca_ut_samp_sc26m, fsca_0p_samp_sc26m)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc26m, fSCA_ut_shl_sc26m, fSCA_0p_exp_sc26m, fSCA_0p_shl_sc26m = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc26m, fSCA_0u_shl_sc26m = fsca0pn_fscaUnTr(fSCA_ut_exp_sc26m, fSCA_ut_shl_sc26m, fSCA_0p_exp_sc26m, fSCA_0p_shl_sc26m, elev_class_inx)
#statistical test
fSCA_0u_exp_sc26m_samp, fSCA_0u_shl_sc26m_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx,sample_size)
statistic_sc26m2,pvalue_sc26m2,average_sample_exp_sc26m2,average_sample_shl_sc26m2,significant_sc26m2 = wilcoxonTest(fSCA_0u_exp_sc26m_samp, fSCA_0u_shl_sc26m_samp)

#y = y = -0.0013047442x + 2.0480940252
tempDec_sc26m = -0.0013047442 * elevationBand_sc + + 2.0480940252

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_sc26m = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)
fSCA_ut_exp_l_sam_sc26m, fSCA_ut_exp_h_sam_sc26m, fSCA_ut_shl_l_sam_sc26m, fSCA_ut_shl_h_sam_sc26m = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx,sample_size)
statistic_ex_sc26m,pvalue_ex_sc26m,average_sample_l_ex_sc26m,average_sample_h_ex_sc26m,significant_ex_sc26m = wilcoxonTest(fSCA_ut_exp_l_sam_sc26m, fSCA_ut_exp_h_sam_sc26m)
statistic_sh_sc26m,pvalue_sh_sc26m,average_sample_l_sh_sc26m,average_sample_h_sh_sc26m,significant_sh_sc26m = wilcoxonTest(fSCA_ut_shl_l_sam_sc26m, fSCA_ut_shl_h_sam_sc26m)
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/significant_sc26m',[significant_sc26m1,significant_sc26m2,significant_ex_sc26m,significant_sh_sc26m])

#calculating topo_dimension from cut northness to plot snow on
snow_temp_vegdens_allExtent_sc26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sc26m.npy')
ouputPath_sc26m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_contour_sc26m.png'
T_sc26m = 'SCW_26March'
plotSnowMap(snow_temp_vegdens_allExtent_sc26m,ouputPath_sc26m,T_sc26m)

lat_long_nad83_sc26m = [np.min(snow_temp_vegdens_allExtent_sc26m[:,0]),np.max(snow_temp_vegdens_allExtent_sc26m[:,0]),np.min(snow_temp_vegdens_allExtent_sc26m[:,1]),np.max(snow_temp_vegdens_allExtent_sc26m[:,1])]

ouputPath_nrt_sc26m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_nrth_contour_sc26m.png'
plotSnowMap_northness(snow_temp_vegdens_allExtent_sc26m,ouputPath_nrt_sc26m,T_sc26m)

#%% load data sagehen creek 17 April 2016
#ascii_grid_veg_index17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc17A.npy')
##calculating topo_dimension from cut northness
#veg_density_sc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/veg_density_sc.npy')
#
#elev_sc17a = ascii_grid_veg_index17a[:,2]
#tempDJF_sc17a = -0.0016 * elev_sc17a + 1.802
#
#veg_snow_temp_density_sc17a = np.vstack([ascii_grid_veg_index17a[:,0],ascii_grid_veg_index17a[:,1].astype(int),ascii_grid_veg_index17a[:,2],
#                                         ascii_grid_veg_index17a[:,3],ascii_grid_veg_index17a[:,4],ascii_grid_veg_index17a[:,5],
#                                         ascii_grid_veg_index17a[:,6],tempDJF_sc17a,veg_density_sc[:,2],tempDJF_sc17a,tempDJF_sc17a]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sc17a',veg_snow_temp_density_sc17a)

## slope less than 30
#snow_temp_vegdens_index_sL30_sc17a = veg_snow_temp_density_sc17a[(veg_snow_temp_density_sc17a[:,5]>-50)&(veg_snow_temp_density_sc17a[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_sc17a',snow_temp_vegdens_index_sL30_sc17a)

snow_temp_vegdens_index_sL30_sc17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_sc17a.npy')
aaaa_test2 = snow_temp_vegdens_index_sL30_sc17a[0:710,:]

elev_sc17a = snow_temp_vegdens_index_sL30_sc17a[:,2]
tempDJF_sc17a = -0.0016 * elev_sc17a + 1.802

snow_temp_vegdens_index_sL30_sc17a_df = pd.DataFrame(snow_temp_vegdens_index_sL30_sc17a, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])
aaaa_test = snow_temp_vegdens_index_sL30_sc17a[700780:710057,:]

elevationBand_sc17a = np.arange(min(snow_temp_vegdens_index_sL30_sc17a_df['z']),max(snow_temp_vegdens_index_sL30_sc17a_df['z']),67)
snow_temp_vegdens_sL30_elevCls_sc17a = classificationElevation(snow_temp_vegdens_index_sL30_sc17a,elevationBand_sc17a,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc17a, fsca_0p_sc17a = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)
fsca_ut_samp_sc17a, fsca_0p_samp_sc17a = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx,sample_size)
#statistical test
statistic_sc17a1,pvalue_sc17a1,average_sample_exp_sc17a1,average_sample_shl_sc17a1,significant_sc17a1 = wilcoxonTest(fsca_ut_samp_sc17a, fsca_0p_samp_sc17a)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a = fsca0pn_fscaUnTr(fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a, elev_class_inx)
#statistical test
fSCA_0u_exp_sc17a_samp, fSCA_0u_shl_sc17a_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx,sample_size)
statistic_sc17a2,pvalue_sc17a2,average_sample_exp_sc17a2,average_sample_shl_sc17a2,significant_sc17a2 = wilcoxonTest(fSCA_0u_exp_sc17a_samp, fSCA_0u_shl_sc17a_samp)

#y = -0.0013x + 2.7982
tempDec_sc17a = -0.0013 * elevationBand_sc17a + 2.7982

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_sc17a = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)
fSCA_ut_exp_l_sam_sc17a, fSCA_ut_exp_h_sam_sc17a, fSCA_ut_shl_l_sam_sc17a, fSCA_ut_shl_h_sam_sc17a = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx,sample_size)
statistic_ex_sc17a,pvalue_ex_sc17a,average_sample_l_ex_sc17a,average_sample_h_ex_sc17a,significant_ex_sc17a = wilcoxonTest(fSCA_ut_exp_l_sam_sc17a, fSCA_ut_exp_h_sam_sc17a)
statistic_sh_sc17a,pvalue_sh_sc17a,average_sample_l_sh_sc17a,average_sample_h_sh_sc17a,significant_sh_sc17a = wilcoxonTest(fSCA_ut_shl_l_sam_sc17a, fSCA_ut_shl_h_sam_sc17a)
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/significant_sc17a',[significant_sc17a1,significant_sc17a2,significant_ex_sc17a,significant_sh_sc17a])

snow_temp_vegdens_allExtent_sc17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sc17a.npy')
ouputPath_sc17a = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_contour_sc17a.png'
T_sc17a = 'SCW_17April'
plotSnowMap(snow_temp_vegdens_allExtent_sc17a,ouputPath_sc17a,T_sc17a)

lat_long_nad83_sc17a = [np.min(snow_temp_vegdens_allExtent_sc17a[:,0]),np.max(snow_temp_vegdens_allExtent_sc17a[:,0]),np.min(snow_temp_vegdens_allExtent_sc17a[:,1]),np.max(snow_temp_vegdens_allExtent_sc17a[:,1])]

ouputPath_nrt_sc17a = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_nrth_contour_sc17a.png'
plotSnowMap_northness(snow_temp_vegdens_allExtent_sc17a,ouputPath_nrt_sc17a,T_sc17a)

#%% load data sagehen creek 18 May 2016
#ascii_grid_veg_index18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc24May.npy')
##calculating topo_dimension from cut northness
#
#veg_density_sc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/veg_density_sc.npy')
#
#elev_sc18m = ascii_grid_veg_index18m[:,2]
#tempDJF_sc18m = -0.0016 * elev_sc18m + 1.802
#
#veg_snow_temp_density_sc18m = np.vstack([ascii_grid_veg_index18m[:,0],ascii_grid_veg_index18m[:,1].astype(int),ascii_grid_veg_index18m[:,2],
#                                         ascii_grid_veg_index18m[:,3],ascii_grid_veg_index18m[:,4],ascii_grid_veg_index18m[:,5],
#                                         ascii_grid_veg_index18m[:,6],tempDJF_sc18m,veg_density_sc[:,2],tempDJF_sc18m,tempDJF_sc18m]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sc18m',veg_snow_temp_density_sc18m)

## slope less than 30
#snow_temp_vegdens_index_sL30_sc18m = veg_snow_temp_density_sc18m[(veg_snow_temp_density_sc18m[:,5]>-50)&(veg_snow_temp_density_sc18m[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_sc18m',snow_temp_vegdens_index_sL30_sc18m)

snow_temp_vegdens_index_sL30_sc18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_sc18m.npy')
aaaa_test3 = snow_temp_vegdens_index_sL30_sc18m[0:710,:]

elev_sc18m = snow_temp_vegdens_index_sL30_sc18m[:,2]
tempDJF_sc18m = -0.0016 * elev_sc18m + 1.802

snow_temp_vegdens_index_sL30_sc18m_df = pd.DataFrame(snow_temp_vegdens_index_sL30_sc18m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_sc18m = np.arange(min(snow_temp_vegdens_index_sL30_sc18m_df['z']),max(snow_temp_vegdens_index_sL30_sc18m_df['z']),67)
snow_temp_vegdens_sL30_elevCls_sc18m = classificationElevation(snow_temp_vegdens_index_sL30_sc18m,elevationBand_sc18m,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc18m, fsca_0p_sc18m = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)
fsca_ut_samp_sc18m, fsca_0p_samp_sc18m = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx,sample_size)
#statistical test
statistic_sc18m1,pvalue_sc18m1,average_sample_exp_sc18m1,average_sample_shl_sc18m1,significant_sc18m1 = wilcoxonTest(fsca_ut_samp_sc18m, fsca_0p_samp_sc18m)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc18m, fSCA_ut_shl_sc18m, fSCA_0p_exp_sc18m, fSCA_0p_shl_sc18m = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc18m, fSCA_0u_shl_sc18m = fsca0pn_fscaUnTr(fSCA_ut_exp_sc18m, fSCA_ut_shl_sc18m, fSCA_0p_exp_sc18m, fSCA_0p_shl_sc18m, elev_class_inx)
#statistical test
fSCA_0u_exp_sc18m_samp, fSCA_0u_shl_sc18m_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx,sample_size)
statistic_sc18m2,pvalue_sc18m2,average_sample_exp_sc18m2,average_sample_shl_sc18m2,significant_sc18m2 = wilcoxonTest(fSCA_0u_exp_sc18m_samp, fSCA_0u_shl_sc18m_samp)

#y = -0.0020x + 5.2481
tempDec_sc18m = -0.002 * elevationBand_sc18m + 5.2481

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_sc18m = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)
fSCA_ut_exp_l_sam_sc18m, fSCA_ut_exp_h_sam_sc18m, fSCA_ut_shl_l_sam_sc18m, fSCA_ut_shl_h_sam_sc18m = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx,sample_size)
statistic_ex_sc18m,pvalue_ex_sc18m,average_sample_l_ex_sc18m,average_sample_h_ex_sc18m,significant_ex_sc18m = wilcoxonTest(fSCA_ut_exp_l_sam_sc18m, fSCA_ut_exp_h_sam_sc18m)
statistic_sh_sc18m,pvalue_sh_sc18m,average_sample_l_sh_sc18m,average_sample_h_sh_sc18m,significant_sh_sc18m = wilcoxonTest(fSCA_ut_shl_l_sam_sc18m, fSCA_ut_shl_h_sam_sc18m)
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/significant_sc18m',[significant_sc18m1,significant_sc18m2,significant_ex_sc18m,significant_sh_sc18m])

snow_temp_vegdens_allExtent_sc18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sc18m.npy')
ouputPath_sc18m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_contour_sc18m.png'
T_sc18m = 'SCW_18May'
plotSnowMap(snow_temp_vegdens_allExtent_sc18m,ouputPath_sc18m,T_sc18m)

lat_long_nad83_sc18m = [np.min(snow_temp_vegdens_allExtent_sc18m[:,0]),np.max(snow_temp_vegdens_allExtent_sc18m[:,0]),np.min(snow_temp_vegdens_allExtent_sc18m[:,1]),np.max(snow_temp_vegdens_allExtent_sc18m[:,1])]

ouputPath_nrt_sc18m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_nrth_contour_sc18m.png'
plotSnowMap_northness(snow_temp_vegdens_allExtent_sc18m,ouputPath_nrt_sc18m,T_sc18m)

#%% load data KREW 2010
#ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
#
#veg_density_kr = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/veg_density_kr.npy')
#
#elev_krew = ascii_grid_veg_indexKrew[:,2]
#tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339 #R² = 0.5074
#
#veg_snow_temp_density_krew = np.vstack([ascii_grid_veg_indexKrew[:,0],ascii_grid_veg_indexKrew[:,1].astype(int),ascii_grid_veg_indexKrew[:,2],
#                                        ascii_grid_veg_indexKrew[:,3],ascii_grid_veg_indexKrew[:,4],ascii_grid_veg_indexKrew[:,5],
#                                        ascii_grid_veg_indexKrew[:,6],tempDJF_gm_krew,veg_density_kr[:,2],tempDJF_gm_krew,tempDJF_gm_krew]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_krew',veg_snow_temp_density_krew)

## slope less than 30
#snow_temp_vegdens_index_sL30_krew = veg_snow_temp_density_krew[(veg_snow_temp_density_krew[:,5]>-50)&(veg_snow_temp_density_krew[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_krew',snow_temp_vegdens_index_sL30_krew)
snow_temp_vegdens_index_sL30_krew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_krew.npy')
aaaa_test1 = snow_temp_vegdens_index_sL30_krew[0:710,:]

#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(snow_temp_vegdens_index_sL30_krew[:,0]),np.max(snow_temp_vegdens_index_sL30_krew[:,0]),np.min(snow_temp_vegdens_index_sL30_krew[:,1]),np.max(snow_temp_vegdens_index_sL30_krew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

elev_krew = snow_temp_vegdens_index_sL30_krew[:,2]
tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339 #R² = 0.5074

snow_temp_vegdens_index_sL30_krew_df = pd.DataFrame(snow_temp_vegdens_index_sL30_krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_krew = np.arange(min(snow_temp_vegdens_index_sL30_krew_df['z']),max(snow_temp_vegdens_index_sL30_krew_df['z']),67)
snow_temp_vegdens_sL30_elevCls_krew = classificationElevation(snow_temp_vegdens_index_sL30_krew,elevationBand_krew,elev_class_inx)

meanTemp_elevClass_krew = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_krew, fsca_0p_krew = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)
fsca_ut_samp_krew, fsca_0p_samp_krew = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx,sample_size)
#statistical test
statistic_krew1,pvalue_krew1,average_sample_exp_krew1,average_sample_shl_krew1,significant_krew1 = wilcoxonTest(fsca_ut_samp_krew, fsca_0p_samp_krew)

# under canopy and 0pen fsca----exposed vs. sheltered
fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_Krew, fSCA_0u_shl_Krew = fsca0pn_fscaUnTr(fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew, elev_class_inx)
#statistical test
sample_size2 = 100      
fSCA_0u_exp_Krew_samp, fSCA_0u_shl_Krew_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx,sample_size2)
statistic_Krew2,pvalue_Krew2,average_sample_exp_Krew2,average_sample_shl_Krew2,significant_krew2 = wilcoxonTest(fSCA_0u_exp_Krew_samp, fSCA_0u_shl_Krew_samp)

#y = -0.003300x + 7.287576
tempDec_krew = -0.0033*elevationBand_krew + 7.287576 

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_krew = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)
fSCA_ut_exp_l_sam_krew, fSCA_ut_exp_h_sam_krew, fSCA_ut_shl_l_sam_krew, fSCA_ut_shl_h_sam_krew = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx,sample_size)
statistic_ex_krew,pvalue_ex_krew,average_sample_l_ex_krew,average_sample_h_ex_krew,significant_ex_krew = wilcoxonTest(fSCA_ut_exp_l_sam_krew, fSCA_ut_exp_h_sam_krew)
statistic_sh_krew,pvalue_sh_krew,average_sample_l_sh_krew,average_sample_h_sh_krew,significant_sh_krew = wilcoxonTest(fSCA_ut_shl_l_sam_krew, fSCA_ut_shl_h_sam_krew)
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/significant_krew',[significant_krew1,significant_krew2,significant_ex_krew,significant_sh_krew])

snow_temp_vegdens_allExtent_krew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_krew.npy')
ouputPath_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_contour_krew.png'
T_krew = 'KREW'
plotSnowMap(snow_temp_vegdens_allExtent_krew,ouputPath_krew,T_krew)

lat_long_nad83_krew = [np.min(snow_temp_vegdens_allExtent_krew[:,0]),np.max(snow_temp_vegdens_allExtent_krew[:,0]),np.min(snow_temp_vegdens_allExtent_krew[:,1]),np.max(snow_temp_vegdens_allExtent_krew[:,1])]

ouputPath_nrt_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snowMap_nrth_contour_krew.png'
plotSnowMap_northness(snow_temp_vegdens_allExtent_krew,ouputPath_nrt_krew,T_krew)

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

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_dec_fsca_all4.png')
#%% ploting fsca vs. DJF temp lapse rate
fSCA_0u = [fsca_ut_sc26m, fsca_0p_sc26m, fsca_ut_sc17a, fsca_0p_sc17a,
           fsca_ut_sc18m, fsca_0p_sc18m, fsca_ut_krew, fsca_0p_krew,
           fsca_ut_jmz, fsca_0p_jmz, fsca_ut_nrc, fsca_0p_nrc]

meanTemp = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
            meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
            meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label0 = ['SCW 26MAR2016-UnderTree','SCW 26MAR2016-0pen','SCW 17APR2016-UnderTree','SCW 17APR2016-0pen',
         'SCW 18MAY2016-UnderTree','SCW 18MAY2016-0pen','Krew 2010-UnderTree','Krew 2010-0pen',
         'JMZ 2010-UnderTree','JMZ 2010 - 0pen','NRC 2010-UnderTree','NRC 2010 - 0pen']
color0 = ['plum','purple','plum','purple','plum','purple',
          'gold','darkred','lightgreen','darkgreen','deepskyblue','navy']#'olive',
markerS = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

marker0 = ['^','o','^','o','^','o','^','o','^','o','^','o']

location0 = ['lower left','lower left','uper right','lower left','lower left','lower left']
xlimit_fsca_DJF = [(-2.4,-1.2),(-2.4,-1.2),(-2.4,-1.2),(0.4,2.8),(-8.3,-3.5),(-10.5,-5)]

plt.subplots(figsize=(60,30)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
    for pnt in range (10):
        if fSCA_0u[2*i][pnt]>fSCA_0u[2*i+1][pnt]:
            plt.scatter(meanTemp[2*i][pnt],fSCA_0u[2*i][pnt]+7, s=25**2, color = 'k', marker = markerS[i][pnt]) #label = label[2*i],
        else:
            plt.scatter(meanTemp[2*i+1][pnt],fSCA_0u[2*i+1][pnt]+7,  s=25**2, color = 'k', marker = markerS[i][pnt]) #label = label[2*i+1],

    plt.plot(meanTemp[2*i],fSCA_0u[2*i], color = color0[2*i], linewidth=4, marker = marker0[2*i], markersize=40) #, 
    plt.plot(meanTemp[2*i+1],fSCA_0u[2*i+1], color = color0[2*i+1], linewidth=4, marker = marker0[2*i+1], markersize=40) #

    plt.ylabel('fSCA%', fontsize=50)
    plt.yticks(fontsize=40)
    plt.xlabel('average temp of DJF (C)', fontsize=40)
    plt.xticks(fontsize=40)
    
    plt.legend([label0[2*i], label0[2*i+1]], fontsize=35, loc = location0[i])
    plt.ylim((0,110))  # adjust the top leaving bottom unchanged
    plt.xlim(xlimit_fsca_DJF[i])

plt.title('fSCA under canopy and in open sites in different temp lapse rate in 4 sites', fontsize=70, y=2.25, x=-0.8) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_jfd_fsca_all6.png')

#%% ploting DJF temp vs terrain features delta fsca
fSCA_d0u = [fSCA_0u_exp_sc26m, fSCA_0u_shl_sc26m, fSCA_0u_exp_Krew, fSCA_0u_shl_Krew, 
            fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a, fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz, 
            fSCA_0u_exp_sc18m, fSCA_0u_shl_sc18m, fSCA_0u_exp_nrc, fSCA_0u_shl_nrc]
meanTemp2 = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label = ['SCW 26Mar2016 - exposed','SCW 26Mar2016 - sheltered',
         'JMZ 2010 - exposed','JMZ 2010 - sheltered',         
         'SCW 17Apr2016 - exposed','SCW 17Apr2016 - sheltered',
         'Krew 2010 - exposed','Krew 2010 - sheltered',
         'SCW 18May2016 - exposed','SCW 18May2016 - sheltered',
         'NRC 2010 - exposed','NRC 2010 - sheltered']

color1 = ['plum','purple','gold','darkred','plum','purple','lightgreen','darkgreen','plum','purple','deepskyblue','navy']#'olive',
#marker = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

marker1 = ['o','^','o','^','o','^','o','^','o','^','o','^']

markerS2 = [significant_sc26m2,significant_krew2,significant_sc17a2,significant_jmz2,significant_sc18m2,significant_nrc2]

location = ['lower left','lower left','lower left','lower left','lower left','lower left']
xlimit_deltAfsca = [(-2.4,-1.2),(0.4,2.8),(-2.4,-1.2),(-8.3,-3.5),(-2.4,-1.2),(-10.5,-5)]
arrow_loc = [-2.4, -2.4, -2.4, 0.3, -8.5, -10]

plt.subplots(figsize=(40,60)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(321+i)
    
#    for pnt in range (10):
#        plt.scatter(meanTemp2[2*i][pnt],fSCA_d0u[2*i][pnt], s=50**2, color = color1[2*i], marker = marker1[2*i]) #label = label[2*i],
#        plt.scatter(meanTemp2[2*i+1][pnt],fSCA_d0u[2*i+1][pnt],  s=50**2, color = color1[2*i+1], marker = marker1[2*i+1]) #label = label[2*i+1],
    for pnt in range (10):
        if fSCA_d0u[2*i][pnt]>fSCA_d0u[2*i+1][pnt]:
            plt.scatter(meanTemp2[2*i][pnt],fSCA_d0u[2*i][pnt]+0.1, s=25**2, color = 'k', marker = markerS2[i][pnt]) #label = label[2*i],
        else:
            plt.scatter(meanTemp2[2*i+1][pnt],fSCA_d0u[2*i+1][pnt]+0.1,  s=25**2, color = 'k', marker = markerS2[i][pnt]) #label = label[2*i+1],

    plt.plot(meanTemp2[2*i],fSCA_d0u[2*i], color = color1[2*i], linewidth=4, marker = marker1[2*i], markersize=40) #, 
    plt.plot(meanTemp2[2*i+1],fSCA_d0u[2*i+1], color = color1[2*i+1], linewidth=4, marker = marker1[2*i+1], markersize=40) #

    plt.ylabel('(fSCA_open - fSCA_under) / fSCA_open', fontsize=45)
    plt.yticks(fontsize=30)
    plt.xlabel('average temp of DJF (C)', fontsize=45)
    plt.xticks(fontsize=30)
    
    plt.legend([label[2*i], label[2*i+1]], fontsize=40, loc = location[i])
    plt.ylim((-1.2,0.8))
    plt.xlim(xlimit_deltAfsca[i])
    
    x = [-11,-8,-5,-2,0,2,4]
    y = [0,0,0,0,0,0,0]
    plt.plot(x,y, color = 'black') #line_temp_ls[2*i],
    
    #plt.arrow(arrow_loc[i], 0.05, 0, 0.5, fc="k", ec="k", head_width=0.1, head_length=0.1)

plt.title('(fSCA_open - fSCA_under)/fSCA_op based on northness in 4 sites', fontsize=60, y=3.45, x=-0.1) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrth_deltaFsca_all8.png')

#%% ploting DJF temp vs terrain features delta fsca
#fSCA_uT_rad_vegDens_krew = [fSCA_unTr_exp_l, fSCA_unTr_shl_l, fSCA_unTr_exp_h, fSCA_unTr_shl_h]

fSCA_0uV = [fSCA_uT_rad_vegDens_sc26m[0], fSCA_uT_rad_vegDens_sc26m[1], fSCA_uT_rad_vegDens_sc26m[2], fSCA_uT_rad_vegDens_sc26m[3], 
            fSCA_uT_rad_vegDens_krew[0], fSCA_uT_rad_vegDens_krew[1], fSCA_uT_rad_vegDens_krew[2], fSCA_uT_rad_vegDens_krew[3],
            fSCA_uT_rad_vegDens_sc17a[0], fSCA_uT_rad_vegDens_sc17a[1], fSCA_uT_rad_vegDens_sc17a[2], fSCA_uT_rad_vegDens_sc17a[3],
            fSCA_uT_rad_vegDens_jmz[0], fSCA_uT_rad_vegDens_jmz[1], fSCA_uT_rad_vegDens_jmz[2], fSCA_uT_rad_vegDens_jmz[3], 
            fSCA_uT_rad_vegDens_sc18m[0], fSCA_uT_rad_vegDens_sc18m[1], fSCA_uT_rad_vegDens_sc18m[2], fSCA_uT_rad_vegDens_sc18m[3],
            fSCA_uT_rad_vegDens_nrc[0], fSCA_uT_rad_vegDens_nrc[1], fSCA_uT_rad_vegDens_nrc[2], fSCA_uT_rad_vegDens_nrc[3]]

meanTemp2 = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_krew,meanTemp_elevClass_krew,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_nrc,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label = ['SC26Mar2016-exposed-lowVD','SC26Mar2016-shelter-lowVD',
         'SC26Mar2016-exposed-highVD','SC26Mar2016-shelter-highVD',
         'Krew2010-exposed-lowVD','Krew2010-shelter-lowVD',
         'Krew2010-exposed-highVD','Krew2010-shelter-highVD',
         'SC17Apr2016-exp-lowVD','SC17Apr2016-shel-lowVD',
         'SC17Apr2016-exp-highVD','SC17Apr2016-shel-highVD',
         'JMZ2010-exposed-lowVD','JMZ2010-shelter-lowVD',
         'JMZ2010-exposed-highVD','JMZ2010-shelter-highVD',
         'SC18May2016-exposed-lowVD','SC18May2016-shelter-lowVD',
         'SC18May2016-exposed-highVD','SC18May2016-shelter-highVD',
         'NRC2010-exposed-lowVD','NRC2010-shelter-lowVD',
         'NRC2010-exposed-highVD','NRC2010-shelter-highVD']

color = ['orchid','purple','plum','darkorchid',
         'red','darkred','gold','orange',
         'orchid','purple','plum','darkorchid',
         'lightgreen','darkgreen','olive','green',
         'orchid','purple','plum','darkorchid',
         'deepskyblue','navy','lightblue','blue']

marker = ['o','s','^','*',
          'o','s','^','*',
          'o','s','^','*',
          'o','s','^','*',
          'o','s','^','*',
          'o','s','^','*']

#markerS3 = [significant_ex_sc26m,significant_ex_sc17a,significant_ex_sc18m,significant_ex_krew,significant_ex_nrc,significant_ex_jmz]
#markerS4 = [significant_sh_sc26m,significant_sh_sc17a,significant_sh_sc18m,significant_sh_krew,significant_sh_nrc,significant_sh_jmz]


location = ['lower left','lower left','lower left','lower left','uper right','lower left']
xlimit_deltAfsca = [(-2.4,-1.2),(0.4,2.8),(-2.4,-1.2),(-8.3,-3.5),(-2.4,-1.2),(-10.5,-5)]

plt.subplots(figsize=(40,60)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(321+i)
    
#    for pnt in range (10):
#        plt.scatter(meanTemp2[4*i][pnt],fSCA_d0u[4*i][pnt], s=30**2, color = color[4*i], marker = marker[4*i]) #label = label[2*i],
#        plt.scatter(meanTemp2[4*i+1][pnt],fSCA_d0u[4*i+1][pnt],  s=25**2, color = color[4*i+1], marker = marker[4*i+1]) #label = label[2*i+1],
#        plt.scatter(meanTemp2[4*i+2][pnt],fSCA_d0u[4*i+2][pnt],  s=35**2, color = color[4*i+2], marker = marker[4*i+2]) #label = label[2*i+1],
#        plt.scatter(meanTemp2[4*i+3][pnt],fSCA_d0u[4*i+3][pnt],  s=40**2, color = color[4*i+3], marker = marker[4*i+3]) #label = label[2*i+1],

#    for pnt in range (10):
#        if fSCA_0uV[4*i][pnt]>fSCA_0uV[4*i+2][pnt]:
#            plt.scatter(meanTemp2[4*i][pnt],fSCA_0uV[4*i][pnt]+7, s=25**2, color = 'k', marker = markerS3[i][pnt]) #label = label[2*i],
#        else:
#            plt.scatter(meanTemp2[4*i+2][pnt],fSCA_0uV[4*i+2][pnt]+7,  s=25**2, color = 'k', marker = markerS3[i][pnt]) #label = label[2*i+1],

    plt.plot(meanTemp2[4*i],fSCA_0uV[4*i], color = color[4*i], linewidth=4, marker = marker[4*i], markersize=30) #, 
    plt.plot(meanTemp2[4*i+1],fSCA_0uV[4*i+1], color = color[4*i+1], linewidth=4, marker = marker[4*i+1], markersize=25) #
    plt.plot(meanTemp2[4*i+2],fSCA_0uV[4*i+2], color = color[4*i+2], linewidth=4, marker = marker[4*i+2], markersize=35) #
    plt.plot(meanTemp2[4*i+3],fSCA_0uV[4*i+3], color = color[4*i+3], linewidth=4, marker = marker[4*i+3], markersize=40) #

    plt.ylabel('fSCA%', fontsize=50)
    plt.yticks(fontsize=50)
    plt.xlabel('average temp of DJF (C)', fontsize=50)
    plt.xticks(fontsize=40)
    
    plt.legend([label[4*i], label[4*i+1], label[4*i+2], label[4*i+3]], fontsize=30, loc = location[i])
    plt.ylim((0,110))
    plt.xlim(xlimit_deltAfsca[i])
    

plt.title('fSCA classification based on northness and  vegetation density (VD)', fontsize=60, y=3.45, x=-0.1) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_nrth_vegDens_all8.png')





#%%
snow_temp_vegdens_index_sL30_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/snow_temp_vegdens_index_sL30_nrc.npy')

lat_nrc =[]
for dx in range (ncols):
    lat_nrc.append(450000+dx)

lon_nrc = []
for dy in range (nrows):
    lon_nrc.append(4425500+dy)
   
lat_rp_nrc = (np.repeat(lat_nrc, nrows)).astype(int)
lon_rp_nrc = (np.tile(lon_nrc, ncols)).astype(int)
minX = np.min(lat_rp_nrc)
minY = np.min(lon_rp_nrc)

snwPres_nrc = np.vstack([snow_temp_vegdens_index_sL30_nrc[:,0],snow_temp_vegdens_index_sL30_nrc[:,1].astype(int),
                         snow_temp_vegdens_index_sL30_nrc[:,5]]).T
snwPres_df_nrc = pd.DataFrame(snwPres_nrc, columns = ['x','y','z'])
indx_DF2 = np.vstack([snwPres_df_nrc['x']-minX,snwPres_df_nrc['y']-minY]).T
indx_nrc2 = indx_DF2[:,1]+ncols*indx_DF2[:,0]
snwPres_df_nrc.index = indx_nrc2.astype(int)

snwDist_nrc = np.vstack([lat_rp_nrc,lon_rp_nrc,lon_rp_nrc]).T
snwDist_nrc[:,2]=-99
snwDist_df_nrc = pd.DataFrame(snwDist_nrc, columns = ['x','y','z'])

N=len(snwDist_df_nrc)
snowMap = np.array([0]*N)
sm=snwPres_df_nrc['z']
print("please wait ...")
print("total length: %s" %N)
nn=0
indx_snwPres = list(snwPres_df_nrc.index)
snowMap[indx_snwPres[:]] = sm[indx_snwPres[:]]

snwMap_test = np.reshape(snowMap,[ncols,nrows])

plt.figure(figsize=(30,20))
cmap = colors.ListedColormap(['darkgray', 'green', 'wheat', 'white', 'blue'])
bounds=[-100,-80,-5,0,50,80]
norm = colors.BoundaryNorm(bounds, cmap.N)

implot = plt.imshow(snwMap, cmap=cmap, norm=norm, extent=lat_long_nad83_nrc)#'gist_earth'
plt.colorbar(orientation='vertical')








