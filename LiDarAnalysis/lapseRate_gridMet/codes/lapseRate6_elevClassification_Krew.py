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
tempMin2010 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2010.nc')
tempMax2010 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2010.nc')
tempMin2009 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2009.nc')
tempMax2009 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2009.nc')
tempMin2010_info = tempMin2010.variables

lat_krew_deg = tempMin2010.variables['lat'][286:290] #0.02083333333333215
long_krew_deg = tempMin2010.variables['lon'][133:136] #0.02083333333333215
#Xmin = 728227.344477	  Ymin = 4363840.587794	  Xmax = 741979.782394	  Ymax = 4371767.

extent_krew = [303500,306000,4100400,4106000]
edgPoints_lat = [303500,304125,304750,305375,306000]
edgPoints_long = [4100400,4102266,4104132,4106000]

product_edg2 = list(itertools.product(edgPoints_lat, edgPoints_long))
edge_m = []
for tupl2 in product_edg2:
    tupl_ls2 = list(tupl2)
    edge_m.append(tupl_ls2)

gridCells_m = [[edge_m[0],edge_m[1],edge_m[4],edge_m[5]],#needs more work
               [edge_m[2],edge_m[3],edge_m[6],edge_m[7]],#needs more work
               [edge_m[4],edge_m[5],edge_m[8],edge_m[9]],#needs more work
               [edge_m[6],edge_m[7],edge_m[10],edge_m[11]],#needs more work
               [edge_m[8],edge_m[9],edge_m[12],edge_m[13]],#needs more work
               [edge_m[10],edge_m[11],edge_m[14],edge_m[15]],#needs more work
               [edge_m[12],edge_m[13],edge_m[16],edge_m[17]],#needs more work
               [edge_m[14],edge_m[15],edge_m[18],edge_m[19]]] #needs more work

#%% elevation gridMet file
elevation_gridMet = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/metdata_elevationdata.nc')
lat_flip = np.flip(elevation_gridMet.variables['lat'][:], axis = 0)
lat_flip_local = lat_flip [286:290]
long_flip_local = elevation_gridMet.variables['lon'] [133:136]

elevation_gridMet_flip = np.flip(elevation_gridMet.variables['elevation'][:], axis = 0)
elevation_gridMet_local = elevation_gridMet_flip[286:290, 133:136]
elevation_ls = [864.4,1310,1424,1194.2	,1727.8,1778.8,1615.6,1926.8,1936.7,1977.8,2215,2187.5]

#%% temperature
#tempMin = tempMin2010.variables['air_temperature'][:]
tempMin_09d = tempMin2009.variables['air_temperature'][335:366, 286:290, 133:136]-273.15
tempMin_10j = tempMin2010.variables['air_temperature'][0:31, 286:290, 133:136]-273.15
tempMin_10f = tempMin2010.variables['air_temperature'][31:59, 286:290, 133:136]-273.15

tempMax_09d = tempMax2009.variables['air_temperature'][335:366, 286:290, 133:136]-273.15
tempMax_10j = tempMax2010.variables['air_temperature'][0:31, 286:290, 133:136]-273.15
tempMax_10f = tempMax2010.variables['air_temperature'][31:59, 286:290, 133:136]-273.15

tempMin_09dm = np.mean(tempMin_09d,0)
tempMin_10jm = np.mean(tempMin_10j,0)
tempMin_10fm = np.mean(tempMin_10f,0)
tempMin_djf = (tempMin_09dm + tempMin_10jm + tempMin_10fm)/3

tempMax_09dm = np.mean(tempMax_09d,0)
tempMax_10jm = np.mean(tempMax_10j,0)
tempMax_10fm = np.mean(tempMax_10f,0)
tempMax_djf = (tempMax_09dm + tempMax_10jm + tempMax_10fm)/3

tempMean_djf = (tempMin_djf + tempMax_djf)/2
meanTemp_djf_ls = [6.47281,4.73633,3.97773,4.87242,2.3286,1.74218,2.76521,1.15866,0.308014,1.0037,-0.286167,-0.693362]
#%% trendline-------fit
z_sc = np.polyfit(elevation_ls,meanTemp_djf_ls, 1)
p_sc = np.poly1d(z_sc)
trendline_sc = p_sc(elevation_ls)
slope, intercept, r_value, p_value, std_err = stats.linregress(elevation_ls,meanTemp_djf_ls)

plt.scatter(elevation_ls,meanTemp_djf_ls)

#%% load data sagehen creek 26 march 2010
ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(ascii_grid_veg_indexKrew[:,0]),np.max(ascii_grid_veg_indexKrew[:,0]),np.min(ascii_grid_veg_indexKrew[:,1]),np.max(ascii_grid_veg_indexKrew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

# slope less than 30
ascii_grid_veg_index_lsKrew = ascii_grid_veg_indexKrew[(ascii_grid_veg_indexKrew[:,5]>-90)&(ascii_grid_veg_indexKrew[:,4]<30)]
aaaa_test = ascii_grid_veg_index_lsKrew[700780:710057,:]

tempDJF_gm = p_sc(ascii_grid_veg_index_lsKrew[:,2])

veg_snow_temp_Krew = np.vstack([ascii_grid_veg_index_lsKrew[:,0],ascii_grid_veg_index_lsKrew[:,1].astype(int),ascii_grid_veg_index_lsKrew[:,2],
                                ascii_grid_veg_index_lsKrew[:,3],ascii_grid_veg_index_lsKrew[:,4],ascii_grid_veg_index_lsKrew[:,5],
                                ascii_grid_veg_index_lsKrew[:,6],tempDJF_gm,tempDJF_gm]).T

veg_snow_temp_dfKrew = pd.DataFrame(veg_snow_temp_Krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand = np.arange(min(veg_snow_temp_dfKrew['z']),max(veg_snow_temp_dfKrew['z']),72)
veg_snow_temp_clsfKrew = classificationElevation(veg_snow_temp_Krew,elevationBand)

meanTemp_elevClass = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfKrew)   
# under canopy and 0pen fsca
fSCA_ut_Krew, fSCA_0p_Krew = fSCAcalculation_elevClassific(veg_snow_temp_clsfKrew)

#fsca open -minus- fsca under canopy
fSCA_0u_Krew = fsca0pn_fscaUnTr(fSCA_ut_Krew,fSCA_0p_Krew)

#%% ploting DJF temp
fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_elevClass,fSCA_0u_Krew, label = 'Krew 2010', s=20**2, color = 'orange')#r2=%.6f'%(r_value_tempFscaKrew**2)
plt.plot(meanTemp_elevClass,fSCA_0u_Krew, color = 'orange')


plt.ylabel('delta fSCA% (0pen-underCanopy)', fontsize=35)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (c)', fontsize=35)
plt.xticks(fontsize=30)

plt.legend(fontsize=20, loc = 'upper right')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_krew.png')











