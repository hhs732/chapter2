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
def cutPart0fMapSort(points_df,minX,maxX,minY,maxY):
    points_df_int = points_df[(points_df['x'] >= minX) & (points_df['x'] <= maxX)]
    points_df_int_sp0 = points_df_int[(points_df_int['y'] >= minY) & (points_df_int['y'] <= maxY)]
    #points_df_int_sp = points_df_int_sp0.values
    points_df_sp_intg = pd.concat([points_df_int_sp0[['x','y']].astype(int),points_df_int_sp0['z']], axis=1)
    points_df_sp_intg.sort_values(by=['x','y'],inplace=True)
    points_df_sp_intg.index=np.arange(0,len(points_df_sp_intg))
    
    return points_df_sp_intg

def cutPart0fASCIISort(points_df,minX,maxX,minY,maxY):
    points_df_int = points_df[(points_df['x'] >= minX) & (points_df['x'] <= maxX)]
    points_df_int_sp0 = points_df_int[(points_df_int['y'] >= minY) & (points_df_int['y'] <= maxY)]
    points_df_sp_intg = pd.concat([points_df_int_sp0[['x','y']].astype(int),points_df_int_sp0['z'],points_df_int_sp0['North'],points_df_int_sp0['slope']], axis=1)
    points_df_sp_intg.sort_values(by=['x','y'],inplace=True)
    points_df_sp_intg.index=np.arange(0,len(points_df_sp_intg))
    
    return points_df_sp_intg

# classification based on lat and long in 100x100 pixels
#veg_snow_temp_sc_df_tst = pd.DataFrame(veg_snow_temp_sc[0:1000,:], columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
#lat_long_df_sc_tst = lat_long_df_sc[['x','y']][0:200]
#cluster100_sc = []
#for row in veg_snow_temp_sc_df.itertuples():
#    indexFinl = np.where(((lat_long_df_sc['x']-row.x<=100) & (lat_long_df_sc['y']-row.y<=100) & (lat_long_df_sc['x']-row.x>=0) & (lat_long_df_sc['y']-row.y>=0)))[0]#
#    if (np.size(indexFinl))>0:
#        cluster100_sc.append([lat_long_df_sc['x'][indexFinl[0]],lat_long_df_sc['y'][indexFinl[0]]])

def classificationGridCells(points_df_ls,gridcell):#0.25
 
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[0][0][0])&(points_df_ls[:,0]<gridcell[0][2][0])&
                      (points_df_ls[:,1]>=gridcell[0][0][1])&(points_df_ls[:,1]<gridcell[0][1][1])]=11
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[1][0][0])&(points_df_ls[:,0]<gridcell[1][2][0])&
                      (points_df_ls[:,1]>=gridcell[1][0][1])&(points_df_ls[:,1]<gridcell[1][1][1])]=12
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[2][0][0])&(points_df_ls[:,0]<gridcell[2][2][0])&
                      (points_df_ls[:,1]>=gridcell[2][0][1])&(points_df_ls[:,1]<gridcell[2][1][1])]=13
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[3][0][0])&(points_df_ls[:,0]<gridcell[3][2][0])&
                      (points_df_ls[:,1]>=gridcell[3][0][1])&(points_df_ls[:,1]<gridcell[3][1][1])]=21
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[4][0][0])&(points_df_ls[:,0]<gridcell[4][2][0])&
                      (points_df_ls[:,1]>=gridcell[4][0][1])&(points_df_ls[:,1]<gridcell[4][1][1])]=22
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[5][0][0])&(points_df_ls[:,0]<gridcell[5][2][0])&
                      (points_df_ls[:,1]>=gridcell[5][0][1])&(points_df_ls[:,1]<gridcell[5][1][1])]=23
    
    pointFile_df = pd.DataFrame(points_df_ls,columns=['x','y','z','north','slope','vegClass','gridClass'])
    
    return pointFile_df.values#, pointFile_df

def classificationElevation(points_df_ls,elevationBand):#0.25
 
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[0])&(points_df_ls[:,2]<elevationBand[1])]=101
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[1])&(points_df_ls[:,2]<elevationBand[2])]=102
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[2])&(points_df_ls[:,2]<elevationBand[3])]=103
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[3])&(points_df_ls[:,2]<elevationBand[4])]=104
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[4])&(points_df_ls[:,2]<elevationBand[5])]=105
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[5])&(points_df_ls[:,2]<elevationBand[6])]=106
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[6])&(points_df_ls[:,2]<elevationBand[7])]=107
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[7])&(points_df_ls[:,2]<elevationBand[8])]=108
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[8])&(points_df_ls[:,2]<elevationBand[9])]=109
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[9])&(points_df_ls[:,2]<elevationBand[10])]=110
    
    points_df_ls[:,8][(points_df_ls[:,8]<=100)]=110
                
    return points_df_ls#, pointFile_df

def fSCAcalculation_elevClassific(veg_snow_temp_elev_clsf_file): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification
    
    ascii_ut_snow101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==101)])
    ascii_ut_snow102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==102)])
    ascii_ut_snow103 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==103)])
    ascii_ut_snow104 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==104)])
    ascii_ut_snow105 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==105)])
    ascii_ut_snow106 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==106)])
    ascii_ut_snow107 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==107)])
    ascii_ut_snow108 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==108)])
    ascii_ut_snow109 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==109)])
    ascii_ut_snow110 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==110)])

    ascii_ut_nsnow101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==101)])
    ascii_ut_nsnow102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==102)])
    ascii_ut_nsnow103 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==103)])
    ascii_ut_nsnow104 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==104)])
    ascii_ut_nsnow105 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==105)])
    ascii_ut_nsnow106 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==106)])
    ascii_ut_nsnow107 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==107)])
    ascii_ut_nsnow108 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==108)])
    ascii_ut_nsnow109 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==109)])
    ascii_ut_nsnow110 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==110)])

    ascii_fSCA_ut101 = 100.*ascii_ut_snow101/(ascii_ut_snow101+ascii_ut_nsnow101)
    ascii_fSCA_ut102 = 100.*ascii_ut_snow102/(ascii_ut_snow102+ascii_ut_nsnow102)
    ascii_fSCA_ut103 = 100.*ascii_ut_snow103/(ascii_ut_snow103+ascii_ut_nsnow103)
    ascii_fSCA_ut104 = 100.*ascii_ut_snow104/(ascii_ut_snow104+ascii_ut_nsnow104)
    ascii_fSCA_ut105 = 100.*ascii_ut_snow105/(ascii_ut_snow105+ascii_ut_nsnow105)
    ascii_fSCA_ut106 = 100.*ascii_ut_snow106/(ascii_ut_snow106+ascii_ut_nsnow106)
    ascii_fSCA_ut107 = 100.*ascii_ut_snow107/(ascii_ut_snow107+ascii_ut_nsnow107)
    ascii_fSCA_ut108 = 100.*ascii_ut_snow108/(ascii_ut_snow108+ascii_ut_nsnow108)
    ascii_fSCA_ut109 = 100.*ascii_ut_snow109/(ascii_ut_snow109+ascii_ut_nsnow109)
    ascii_fSCA_ut110 = 100.*ascii_ut_snow110/(ascii_ut_snow110+ascii_ut_nsnow110)
    
    fSCA_unTr = [ascii_fSCA_ut101, ascii_fSCA_ut102, ascii_fSCA_ut103, ascii_fSCA_ut104, 
                 ascii_fSCA_ut105, ascii_fSCA_ut106, ascii_fSCA_ut107, ascii_fSCA_ut108, 
                 ascii_fSCA_ut109, ascii_fSCA_ut110]
    
    
    ascii_0p_snow101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==101)])
    ascii_0p_snow102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==102)])
    ascii_0p_snow103 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==103)])
    ascii_0p_snow104 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==104)])
    ascii_0p_snow105 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==105)])
    ascii_0p_snow106 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==106)])
    ascii_0p_snow107 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==107)])
    ascii_0p_snow108 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==108)])
    ascii_0p_snow109 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==109)])
    ascii_0p_snow110 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==110)])

    ascii_0p_nsnow101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==101)])
    ascii_0p_nsnow102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==102)])
    ascii_0p_nsnow103 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==103)])
    ascii_0p_nsnow104 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==104)])
    ascii_0p_nsnow105 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==105)])
    ascii_0p_nsnow106 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==106)])
    ascii_0p_nsnow107 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==107)])
    ascii_0p_nsnow108 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==108)])
    ascii_0p_nsnow109 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==109)])
    ascii_0p_nsnow110 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==110)])

    ascii_fSCA_0p101 = 100.*ascii_0p_snow101/(ascii_0p_snow101+ascii_0p_nsnow101)
    ascii_fSCA_0p102 = 100.*ascii_0p_snow102/(ascii_0p_snow102+ascii_0p_nsnow102)
    ascii_fSCA_0p103 = 100.*ascii_0p_snow103/(ascii_0p_snow103+ascii_0p_nsnow103)
    ascii_fSCA_0p104 = 100.*ascii_0p_snow104/(ascii_0p_snow104+ascii_0p_nsnow104)
    ascii_fSCA_0p105 = 100.*ascii_0p_snow105/(ascii_0p_snow105+ascii_0p_nsnow105)
    ascii_fSCA_0p106 = 100.*ascii_0p_snow106/(ascii_0p_snow106+ascii_0p_nsnow106)
    ascii_fSCA_0p107 = 100.*ascii_0p_snow107/(ascii_0p_snow107+ascii_0p_nsnow107)
    ascii_fSCA_0p108 = 100.*ascii_0p_snow108/(ascii_0p_snow108+ascii_0p_nsnow108)
    ascii_fSCA_0p109 = 100.*ascii_0p_snow109/(ascii_0p_snow109+ascii_0p_nsnow109)
    ascii_fSCA_0p110 = 100.*ascii_0p_snow110/(ascii_0p_snow110+ascii_0p_nsnow110)
    
    fSCA_0pn = [ascii_fSCA_0p101, ascii_fSCA_0p102, ascii_fSCA_0p103, ascii_fSCA_0p104, 
                 ascii_fSCA_0p105, ascii_fSCA_0p106, ascii_fSCA_0p107, ascii_fSCA_0p108, 
                 ascii_fSCA_0p109, ascii_fSCA_0p110]
    
    return fSCA_unTr, fSCA_0pn

def fsca0pn_fscaUnTr(fSCA_ut,fSCA_0p):
    fSCA_0u0 = fSCA_0p[0] - fSCA_ut[0]
    fSCA_0u1 = fSCA_0p[1] - fSCA_ut[1]
    fSCA_0u2 = fSCA_0p[2] - fSCA_ut[2]
    fSCA_0u3 = fSCA_0p[3] - fSCA_ut[3]
    fSCA_0u4 = fSCA_0p[4] - fSCA_ut[4]
    fSCA_0u5 = fSCA_0p[5] - fSCA_ut[5]
    fSCA_0u6 = fSCA_0p[6] - fSCA_ut[6]
    fSCA_0u7 = fSCA_0p[7] - fSCA_ut[7]
    fSCA_0u8 = fSCA_0p[8] - fSCA_ut[8]
    fSCA_0u9 = fSCA_0p[9] - fSCA_ut[9]
        
    fSCA_0p_ut = [fSCA_0u0, fSCA_0u1, fSCA_0u2, fSCA_0u3, fSCA_0u4, fSCA_0u5, fSCA_0u6, fSCA_0u7, fSCA_0u8, fSCA_0u9]
    
    return fSCA_0p_ut

def calculationMeanTempforEachElevClassification(veg_snow_temp_elev_clsf_file):
    temp101_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==101)][:,7])
    temp102_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==102)][:,7])
    temp103_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==103)][:,7])
    temp104_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==104)][:,7])
    temp105_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==105)][:,7])
    temp106_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==106)][:,7])
    temp107_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==107)][:,7])
    temp108_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==108)][:,7])
    temp109_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==109)][:,7])
    temp110_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==110)][:,7])
    
    meanTemp = [temp101_mean, temp102_mean, temp103_mean, temp104_mean, temp105_mean, temp106_mean, temp107_mean, temp108_mean, temp109_mean, temp110_mean]
    
    return meanTemp
#%% lat long manipulation
tempMin2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2016.nc')
tempMax2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2016.nc')
tempMin2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2015.nc')
tempMax2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2015.nc')
tempMin2016_info = tempMin2016.variables

lat_sc_deg = tempMin2016.variables['lat'][344:346] #0.02083333333333215
#Xmin = 728227.344477	  Ymin = 4363840.587794	  Xmax = 741979.782394	  Ymax = 4371767.484312
extent_sc = [730000,741000,4362300,4372000]
edgPoints_sc_lat = [730000,733666.67,737333.33,741000]
edgPoints_sc_long = [4362300,4367150,4372000]

product_edg2 = list(itertools.product(edgPoints_sc_lat, edgPoints_sc_long))
sc_edge_m = []
for tupl2 in product_edg2:
    tupl_ls2 = list(tupl2)
    sc_edge_m.append(tupl_ls2)

gridCells_m = [[sc_edge_m[0],sc_edge_m[1],sc_edge_m[3],sc_edge_m[4]],
               [sc_edge_m[3],sc_edge_m[4],sc_edge_m[6],sc_edge_m[7]],
               [sc_edge_m[6],sc_edge_m[7],sc_edge_m[9],sc_edge_m[10]],
               [sc_edge_m[1],sc_edge_m[2],sc_edge_m[4],sc_edge_m[5]],
               [sc_edge_m[4],sc_edge_m[5],sc_edge_m[7],sc_edge_m[8]],
               [sc_edge_m[7],sc_edge_m[8],sc_edge_m[10],sc_edge_m[11]]]   

#%% elevation gridMet file
elevation_gridMet = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/metdata_elevationdata.nc')
lat_flip = np.flip(elevation_gridMet.variables['lat'][:], axis = 0)
lat_flip_sc = lat_flip [344:346]
long_flip_sc = elevation_gridMet.variables['lon'] [107:110]

elevation_gridMet_flip = np.flip(elevation_gridMet.variables['elevation'][:], axis = 0)
elevation_gridMet_sc = elevation_gridMet_flip[344:346, 107:110]
elevation_sc_ls = [2228.7, 2098.8, 1985.3, 2244.2, 2091.2, 1973.2]

#%% temperature
#tempMin = tempMin2016.variables['air_temperature'][:]
tempMin_sc15d = tempMin2015.variables['air_temperature'][335:366, 344:346, 107:110]-273.15
tempMin_sc16j = tempMin2016.variables['air_temperature'][0:31, 344:346, 107:110]-273.15
tempMin_sc16f = tempMin2016.variables['air_temperature'][31:60, 344:346, 107:110]-273.15

tempMax_sc15d = tempMax2015.variables['air_temperature'][335:366, 344:346, 107:110]-273.15
tempMax_sc16j = tempMax2016.variables['air_temperature'][0:31, 344:346, 107:110]-273.15
tempMax_sc16f = tempMax2016.variables['air_temperature'][31:60, 344:346, 107:110]-273.15

tempMin_sc15dm = np.mean(tempMin_sc15d,0)
tempMin_sc16jm = np.mean(tempMin_sc16j,0)
tempMin_sc16fm = np.mean(tempMin_sc16f,0)
tempMin_sc_djf = (tempMin_sc15dm + tempMin_sc16jm + tempMin_sc16fm)/3

tempMax_sc15dm = np.mean(tempMax_sc15d,0)
tempMax_sc16jm = np.mean(tempMax_sc16j,0)
tempMax_sc16fm = np.mean(tempMax_sc16f,0)
tempMax_sc_djf = (tempMax_sc15dm + tempMax_sc16jm + tempMax_sc16fm)/3

tempMean_sc_djf = (tempMin_sc_djf + tempMax_sc_djf)/2
meanTemp_djf_ls = [-0.725272,-0.482033,-0.790013,-0.729759,-0.80674,-1.1651]
#%% trendline-------fit
z_sc = np.polyfit(elevation_sc_ls,meanTemp_djf_ls, 1)
p_sc = np.poly1d(z_sc)
trendline_sc = p_sc(elevation_sc_ls)
slope_sc, intercept_sc, r_value_sc, p_value_sc, std_err_sc = stats.linregress(elevation_sc_ls,meanTemp_djf_ls)

#plt.scatter(elevation_sc_ls,meanTemp_djf_ls)

#%% load data sagehen creek 26 march 2016
ascii_index_pr26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc26m.npy')
#calculating topo_dimension from cut northness

ascii_grid_veg_index26m = classificationGridCells(ascii_index_pr26m,gridCells_m)
# slope less than 30
ascii_grid_veg_index_ls26m = ascii_grid_veg_index26m[(ascii_grid_veg_index26m[:,5]>-90)&(ascii_grid_veg_index26m[:,4]<30)]

tempDJF_sc = p_sc(ascii_grid_veg_index_ls26m[:,2])

veg_snow_temp_sc26m = np.vstack([ascii_grid_veg_index_ls26m[:,0],ascii_grid_veg_index_ls26m[:,1].astype(int),ascii_grid_veg_index_ls26m[:,2],
                              ascii_grid_veg_index_ls26m[:,3],ascii_grid_veg_index_ls26m[:,4],ascii_grid_veg_index_ls26m[:,5],
                              ascii_grid_veg_index_ls26m[:,6],tempDJF_sc,tempDJF_sc]).T

veg_snow_temp_sc_df26m = pd.DataFrame(veg_snow_temp_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc26m[700780:710057,:]

elevationBand_sc = np.arange(min(veg_snow_temp_sc_df26m['z']),max(veg_snow_temp_sc_df26m['z']),70)
veg_snow_temp_sc_clsf26m = classificationElevation(veg_snow_temp_sc26m,elevationBand_sc)

meanTemp_elevClass_sc26m = calculationMeanTempforEachElevClassification(veg_snow_temp_sc_clsf26m)   
#%% under canopy and 0pen fsca
fSCA_ut_sc26m, fSCA_0p_sc26m = fSCAcalculation_elevClassific(veg_snow_temp_sc_clsf26m)

#fsca open -minus- fsca under canopy
fSCA_0u_sc26m = fsca0pn_fscaUnTr(fSCA_ut_sc26m,fSCA_0p_sc26m)

#%% load data sagehen creek 17 April 2016
ascii_index_pr26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc17A.npy')
#calculating topo_dimension from cut northness

ascii_grid_veg_index26m = classificationGridCells(ascii_index_pr26m,gridCells_m)
# slope less than 30
ascii_grid_veg_index_ls26m = ascii_grid_veg_index26m[(ascii_grid_veg_index26m[:,5]>-90)&(ascii_grid_veg_index26m[:,4]<30)]

tempDJF_sc = p_sc(ascii_grid_veg_index_ls26m[:,2])

veg_snow_temp_sc26m = np.vstack([ascii_grid_veg_index_ls26m[:,0],ascii_grid_veg_index_ls26m[:,1].astype(int),ascii_grid_veg_index_ls26m[:,2],
                              ascii_grid_veg_index_ls26m[:,3],ascii_grid_veg_index_ls26m[:,4],ascii_grid_veg_index_ls26m[:,5],
                              ascii_grid_veg_index_ls26m[:,6],tempDJF_sc,tempDJF_sc]).T

veg_snow_temp_sc_df26m = pd.DataFrame(veg_snow_temp_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc26m[700780:710057,:]

elevationBand_sc = np.arange(min(veg_snow_temp_sc_df26m['z']),max(veg_snow_temp_sc_df26m['z']),70)
veg_snow_temp_sc_clsf26m = classificationElevation(veg_snow_temp_sc26m,elevationBand_sc)

meanTemp_elevClass_sc26m = calculationMeanTempforEachElevClassification(veg_snow_temp_sc_clsf26m)   
#%% under canopy and 0pen fsca
fSCA_ut_sc26m, fSCA_0p_sc26m = fSCAcalculation_elevClassific(veg_snow_temp_sc_clsf26m)

#fsca open -minus- fsca under canopy
fSCA_0u_sc26m = fsca0pn_fscaUnTr(fSCA_ut_sc26m,fSCA_0p_sc26m)





#%% ploting DJF temp
#tempFsca_sc26M = np.polyfit(meanTemp_elevClass_sc,fSCA_0u_sc, 1)
#fun_tempFsca_sc26M = np.poly1d(tempFsca_sc26M)
#trendline_tempFsca_sc26M = fun_tempFsca_sc26M(fSCA_0u_sc)
#slope_tempFsca26M, intercept_tempFsca26M, r_value_tempFsca26M, p_value_tempFsca26M, std_err_tempFsca26M = stats.linregress(meanTemp_elevClass_sc, fSCA_0u_sc)

fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_elevClass_sc,fSCA_0u_sc, label = 'Sagehen Creek 26March2016', s=20**2)#r2=%.6f'%(r_value_tempFsca26M**2)

plt.ylabel('delta fSCA% (0pen-underCanopy)', fontsize=35)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (c)', fontsize=35)
plt.xticks(fontsize=30)

plt.legend(fontsize=30, loc = 'upper right')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_scM26.png')











