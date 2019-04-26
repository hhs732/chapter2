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
#%% lat long manipulation
tempMin2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2016.nc')
tempMax2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2016.nc')
tempMin2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2015.nc')
tempMax2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2015.nc')
tempMin2016_info = tempMin2016.variables

lat_sc_deg = tempMin2016.variables['lat'][344:346] #0.02083333333333215
#Xmin = 728227.344477	  Ymin = 4363840.587794	  Xmax = 741979.782394	  Ymax = 4371767.484312

lat_sc_minEdg = lat_sc_deg-0.02083333333333215
lat_sc_maxEdg = lat_sc_deg+0.02083333333333215
lat_sc_edg = np.unique(np.append(lat_sc_minEdg,lat_sc_maxEdg))

long_sc_deg = tempMin2016.variables['lon'][107:110]#0.020833333333335702
long_sc_minEdg = long_sc_deg-0.020833333333335702
long_sc_maxEdg = long_sc_deg+0.020833333333335702
long_sc_edg = np.unique(np.append(long_sc_minEdg,long_sc_maxEdg))

product_edg = list(itertools.product(lat_sc_edg, long_sc_edg))
sc_edge = []
for tupl in product_edg:
    tupl_ls = list(tupl)
    sc_edge.append(tupl_ls)

product_coords = list(itertools.product(lat_sc_deg, long_sc_deg))
sc_coords = []
for tupl2 in product_coords:
    tupl_ls2 = list(tupl2)
    sc_coords.append(tupl_ls2)

gridCells = [[sc_edge[0],sc_edge[1],sc_edge[4],sc_edge[5]],
             [sc_edge[1],sc_edge[2],sc_edge[5],sc_edge[6]],
             [sc_edge[2],sc_edge[3],sc_edge[6],sc_edge[7]],
             [sc_edge[4],sc_edge[5],sc_edge[8],sc_edge[9]],
             [sc_edge[5],sc_edge[6],sc_edge[9],sc_edge[10]],
             [sc_edge[6],sc_edge[7],sc_edge[10],sc_edge[11]]] 
gridCellAuto = []
for indx1 in range(len(sc_coords)):
    oneSqr = []   
    for indx2 in range(len(sc_edge)):
        if((sc_edge[indx2][0]-sc_coords[indx1][0])<0.0416667 and (sc_edge[indx2][0]-sc_coords[indx1][0])>=-0.0416667 and (sc_edge[indx2][1]-sc_coords[indx1][1])<0.0416667 and (sc_edge[indx2][1]-sc_coords[indx1][1])>=-0.0416667):
            oneSqr.append(sc_edge[indx2])
    gridCellAuto.append(oneSqr)

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

ncolms_sc = ((extent_sc[1]-extent_sc[0])/100)+1
latitude_sc =[]
for ltx in range (ncolms_sc):
    latitude_sc.append(extent_sc[0]+ltx*100)

nrows_sc = ((extent_sc[3]-extent_sc[2])/100)+1
longitude_sc = []
for ylg in range (nrows_sc):
    longitude_sc.append(extent_sc[2]+ylg*100)

latitude_rp_sc = pd.DataFrame((np.tile(latitude_sc, nrows_sc)).astype(int),columns=['x'])
longitude_rp_sc = pd.DataFrame((np.repeat(longitude_sc, ncolms_sc)).astype(int),columns=['y'])
lat_long_df_sc = pd.concat([latitude_rp_sc,longitude_rp_sc],axis=1)
     
#%% 
elevation_gridMet = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/metdata_elevationdata.nc')
elevation_gridMet_flip = np.flip(elevation_gridMet.variables['elevation'][:], axis = 0)
elevation_gridMet_sc = elevation_gridMet_flip[344:346, 107:110]
elevation_sc_ls = [1985.3, 2098.8, 2228.7, 1973.2, 2091.2	, 2244.2	]

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
meanTemp_djf_ls = [-0.790013, -0.482033, -0.725272, -1.1651, -0.80674, -0.729759]
#%% trendline------------------fit
z_sc = np.polyfit(elevation_sc_ls,meanTemp_djf_ls, 1)
p_sc = np.poly1d(z_sc)
trendline_sc = p_sc(elevation_sc_ls)
slope_sc, intercept_sc, r_value_sc, p_value_sc, std_err_sc = stats.linregress(elevation_sc_ls,meanTemp_djf_ls)

#plt.scatter(elevation_sc_ls,meanTemp_djf_ls)

#%% load data
ascii_index_pr = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_sc26m.npy')
#calculating topo_dimension from cut northness

ascii_grid_veg_index = classificationGridCells(ascii_index_pr,gridCells_m)
# slope less than 30
ascii_grid_veg_index_ls = ascii_grid_veg_index[(ascii_grid_veg_index[:,5]>-90)&(ascii_grid_veg_index[:,4]<30)]

tempDJF_sc = p_sc(ascii_grid_veg_index_ls[:,2])

veg_snow_temp_sc = np.vstack([ascii_grid_veg_index_ls[:,0],ascii_grid_veg_index_ls[:,1].astype(int),ascii_grid_veg_index_ls[:,2],
                              ascii_grid_veg_index_ls[:,3],ascii_grid_veg_index_ls[:,4],ascii_grid_veg_index_ls[:,5],
                              ascii_grid_veg_index_ls[:,6],tempDJF_sc,tempDJF_sc]).T

veg_snow_temp_sc_df = pd.DataFrame(veg_snow_temp_sc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
aaaa_test = veg_snow_temp_sc[700780:710057,:]

#veg_snow_temp_sc_df_tst = pd.DataFrame(veg_snow_temp_sc[0:1000,:], columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
#lat_long_df_sc_tst = lat_long_df_sc[['x','y']][0:200]
cluster100_sc = []
for row in veg_snow_temp_sc_df.itertuples():
    indexFinl = np.where(((lat_long_df_sc['x']-row.x<=100) & (lat_long_df_sc['y']-row.y<=100) & (lat_long_df_sc['x']-row.x>=0) & (lat_long_df_sc['y']-row.y>=0)))[0]#
    if (np.size(indexFinl))>0:
        cluster100_sc.append([lat_long_df_sc['x'][indexFinl[0]],lat_long_df_sc['y'][indexFinl[0]]])
    
#%% under canopy fsca
ascii_ut_snow11 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==77)&(ascii_grid_veg_index_ls[:,6]==11)])
ascii_ut_snow12 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==77)&(ascii_grid_veg_index_ls[:,6]==12)])
ascii_ut_snow13 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==77)&(ascii_grid_veg_index_ls[:,6]==13)])
ascii_ut_snow21 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==77)&(ascii_grid_veg_index_ls[:,6]==21)])
ascii_ut_snow22 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==77)&(ascii_grid_veg_index_ls[:,6]==22)])
ascii_ut_snow23 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==77)&(ascii_grid_veg_index_ls[:,6]==23)])

ascii_ut_nsnow11 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-7)&(ascii_grid_veg_index_ls[:,6]==11)])
ascii_ut_nsnow12 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-7)&(ascii_grid_veg_index_ls[:,6]==12)])
ascii_ut_nsnow13 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-7)&(ascii_grid_veg_index_ls[:,6]==13)])
ascii_ut_nsnow21 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-7)&(ascii_grid_veg_index_ls[:,6]==21)])
ascii_ut_nsnow22 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-7)&(ascii_grid_veg_index_ls[:,6]==22)])
ascii_ut_nsnow23 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-7)&(ascii_grid_veg_index_ls[:,6]==23)])

ascii_fSCA_ut11 = 100.*ascii_ut_snow11/(ascii_ut_snow11+ascii_ut_nsnow11)
ascii_fSCA_ut12 = 100.*ascii_ut_snow12/(ascii_ut_snow12+ascii_ut_nsnow12)
ascii_fSCA_ut13 = 100.*ascii_ut_snow13/(ascii_ut_snow13+ascii_ut_nsnow13)
ascii_fSCA_ut21 = 100.*ascii_ut_snow21/(ascii_ut_snow21+ascii_ut_nsnow21)
ascii_fSCA_ut22 = 100.*ascii_ut_snow22/(ascii_ut_snow22+ascii_ut_nsnow22)
ascii_fSCA_ut23 = 100.*ascii_ut_snow23/(ascii_ut_snow23+ascii_ut_nsnow23)
#%% open fsca
ascii_0p_snow11 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==44)&(ascii_grid_veg_index_ls[:,6]==11)])
ascii_0p_snow12 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==44)&(ascii_grid_veg_index_ls[:,6]==12)])
ascii_0p_snow13 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==44)&(ascii_grid_veg_index_ls[:,6]==13)])
ascii_0p_snow21 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==44)&(ascii_grid_veg_index_ls[:,6]==21)])
ascii_0p_snow22 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==44)&(ascii_grid_veg_index_ls[:,6]==22)])
ascii_0p_snow23 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==44)&(ascii_grid_veg_index_ls[:,6]==23)])

ascii_0p_nsnow11 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-4)&(ascii_grid_veg_index_ls[:,6]==11)])
ascii_0p_nsnow12 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-4)&(ascii_grid_veg_index_ls[:,6]==12)])
ascii_0p_nsnow13 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-4)&(ascii_grid_veg_index_ls[:,6]==13)])
ascii_0p_nsnow21 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-4)&(ascii_grid_veg_index_ls[:,6]==21)])
ascii_0p_nsnow22 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-4)&(ascii_grid_veg_index_ls[:,6]==22)])
ascii_0p_nsnow23 = len(ascii_grid_veg_index_ls[(ascii_grid_veg_index_ls[:,5]==-4)&(ascii_grid_veg_index_ls[:,6]==23)])

ascii_fSCA_0p11 = 100.*ascii_0p_snow11/(ascii_0p_snow11+ascii_0p_nsnow11)
ascii_fSCA_0p12 = 100.*ascii_0p_snow12/(ascii_0p_snow12+ascii_0p_nsnow12)
ascii_fSCA_0p13 = 100.*ascii_0p_snow13/(ascii_0p_snow13+ascii_0p_nsnow13)
ascii_fSCA_0p21 = 100.*ascii_0p_snow21/(ascii_0p_snow21+ascii_0p_nsnow21)
ascii_fSCA_0p22 = 100.*ascii_0p_snow22/(ascii_0p_snow22+ascii_0p_nsnow22)
ascii_fSCA_0p23 = 100.*ascii_0p_snow23/(ascii_0p_snow23+ascii_0p_nsnow23)
#%% fsca open -- fsca under canopy
fSCA_0u11 = ascii_fSCA_0p11 - ascii_fSCA_ut11
fSCA_0u12 = ascii_fSCA_0p12 - ascii_fSCA_ut12
fSCA_0u13 = ascii_fSCA_0p13 - ascii_fSCA_ut13
fSCA_0u21 = ascii_fSCA_0p21 - ascii_fSCA_ut21
fSCA_0u22 = ascii_fSCA_0p22 - ascii_fSCA_ut22
fSCA_0u23 = ascii_fSCA_0p23 - ascii_fSCA_ut23

fSCA_0u = [fSCA_0u11, fSCA_0u12, fSCA_0u13, fSCA_0u21, fSCA_0u22, fSCA_0u23]

#%% ploting DJF temp
z = np.polyfit(meanTemp_djf_ls,fSCA_0u, 1)
p = np.poly1d(z)
trendline = p(meanTemp_djf_ls)
slope, intercept, r_value, p_value, std_err = stats.linregress(meanTemp_djf_ls, fSCA_0u)

fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_djf_ls,fSCA_0u, label = 'Sagehen Creek 26March2016  r2=%.6f'%(r_value**2), s=20**2)

plt.plot(meanTemp_djf_ls,trendline,"r--")
plt.ylabel('delta fSCA% (0pen-underCanopy)', fontsize=30)
plt.yticks(fontsize=20)
plt.xlabel('average temp of DJF (c)', fontsize=30)
plt.xticks(fontsize=20)

plt.legend(fontsize=25, loc = 'upper left')

# the line equation:
print "y=%.6fx+(%.6f)"%(z[0],z[1])
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/temp_lsr_scM26.png')

#%% swr
swr2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/srad_2015.nc')
swr2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/srad_2016.nc')
swr_info = swr2015.variables

swr_sc15d = swr2015.variables['surface_downwelling_shortwave_flux_in_air'][335:366, 344:346, 107:110]-273.15
swr_sc16j = swr2016.variables['surface_downwelling_shortwave_flux_in_air'][0:31, 344:346, 107:110]-273.15
swr_sc16f = swr2016.variables['surface_downwelling_shortwave_flux_in_air'][31:60, 344:346, 107:110]-273.15

swr_sc15dm = -1*np.mean(swr_sc15d,0)
swr_sc16jm = -1*np.mean(swr_sc16j,0)
swr_sc16fm = -1*np.mean(swr_sc16f,0)
swr_sc_djf = (swr_sc15dm + swr_sc16jm + swr_sc16fm)/3
swr_sc_djf_ls = [160.46, 159.635, 158.853, 159.385, 157.935, 156.551]
#%% ploting
zswr = np.polyfit(swr_sc_djf_ls,fSCA_0u, 1)
pswr = np.poly1d(zswr)
trendline_swr = pswr(swr_sc_djf_ls)
slope_swr, intercept_swr, r_value_swr, p_value_swr, std_err_swr = stats.linregress(swr_sc_djf_ls, fSCA_0u)

fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(swr_sc_djf_ls,fSCA_0u, label = 'Sagehen Creek 26March2016  r2=%.6f'%(r_value_swr**2), s=20**2, color = 'red')

plt.plot(swr_sc_djf_ls,trendline_swr,"k-")
plt.ylabel('delta fSCA% (0pen-underCanopy)', fontsize=30)
plt.yticks(fontsize=20)
plt.xlabel('average swr in DJF (w/m2)', fontsize=30)
plt.xticks(fontsize=20)

plt.legend(fontsize=25, loc = 'upper right')

# the line equation:
print "y=%.6fx+(%.6f)"%(z[0],z[1])
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/swr_lsr_sc26M.png')
#%% relationship between swr and temp in each pixel and each month
tempMean_sc_d00 = (tempMin_sc15d[:,0,0] + tempMax_sc15d[:,0,0])/2
tempMean_sc_d01 = (tempMin_sc15d[:,0,1] + tempMax_sc15d[:,0,1])/2
tempMean_sc_d02 = (tempMin_sc15d[:,0,2] + tempMax_sc15d[:,0,2])/2
tempMean_sc_d10 = (tempMin_sc15d[:,1,0] + tempMax_sc15d[:,1,0])/2
tempMean_sc_d11 = (tempMin_sc15d[:,1,1] + tempMax_sc15d[:,1,1])/2
tempMean_sc_d12 = (tempMin_sc15d[:,1,2] + tempMax_sc15d[:,1,2])/2

tempMean_sc_j00 = (tempMin_sc16j[:,0,0] + tempMax_sc16j[:,0,0])/2
tempMean_sc_j01 = (tempMin_sc16j[:,0,1] + tempMax_sc16j[:,0,1])/2
tempMean_sc_j02 = (tempMin_sc16j[:,0,2] + tempMax_sc16j[:,0,2])/2
tempMean_sc_j10 = (tempMin_sc16j[:,1,0] + tempMax_sc16j[:,1,0])/2
tempMean_sc_j11 = (tempMin_sc16j[:,1,1] + tempMax_sc16j[:,1,1])/2
tempMean_sc_j12 = (tempMin_sc16j[:,1,2] + tempMax_sc16j[:,1,2])/2

tempMean_sc_f00 = (tempMin_sc16f[:,0,0] + tempMax_sc16f[:,0,0])/2
tempMean_sc_f01 = (tempMin_sc16f[:,0,1] + tempMax_sc16f[:,0,1])/2
tempMean_sc_f02 = (tempMin_sc16f[:,0,2] + tempMax_sc16f[:,0,2])/2
tempMean_sc_f10 = (tempMin_sc16f[:,1,0] + tempMax_sc16f[:,1,0])/2
tempMean_sc_f11 = (tempMin_sc16f[:,1,1] + tempMax_sc16f[:,1,1])/2
tempMean_sc_f12 = (tempMin_sc16f[:,1,2] + tempMax_sc16f[:,1,2])/2

swr_sc_d00 = -1*swr_sc15d[:,0,0]
swr_sc_d01 = -1*swr_sc15d[:,0,1]
swr_sc_d02 = -1*swr_sc15d[:,0,2]
swr_sc_d10 = -1*swr_sc15d[:,1,0] 
swr_sc_d11 = -1*swr_sc15d[:,1,1]
swr_sc_d12 = -1*swr_sc15d[:,1,2]

swr_sc_j00 = -1*swr_sc16j[:,0,0]
swr_sc_j01 = -1*swr_sc16j[:,0,1]
swr_sc_j02 = -1*swr_sc16j[:,0,2]
swr_sc_j10 = -1*swr_sc16j[:,1,0]
swr_sc_j11 = -1*swr_sc16j[:,1,1]
swr_sc_j12 = -1*swr_sc16j[:,1,2]

swr_sc_f00 = -1*swr_sc16f[:,0,0]
swr_sc_f01 = -1*swr_sc16f[:,0,1]
swr_sc_f02 = -1*swr_sc16f[:,0,2]
swr_sc_f10 = -1*swr_sc16f[:,1,0]
swr_sc_f11 = -1*swr_sc16f[:,1,1]
swr_sc_f12 = -1*swr_sc16f[:,1,2]
#%%finding trendline properties
#ztsj00 = np.polyfit(tempMean_sc_j00,swr_sc_j00, 1)
#ptsj00 = np.poly1d(ztsj00)
#trendline_tsj00 = ptsj00(tempMean_sc_j00)

tempMean_sc = [tempMean_sc_d00,tempMean_sc_d01,tempMean_sc_d02,
               tempMean_sc_d10,tempMean_sc_d11,tempMean_sc_d12]
swr_sc = [swr_sc_d00,swr_sc_d01,swr_sc_d02,
          swr_sc_d10,swr_sc_d11,swr_sc_d12]
gridcells_num = [0.0, 0.1, 0.2, 1.0, 1.1, 1.2]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))
plt.title('Sagehen Creek temp-swr_dec', fontsize=35)

for i in range (6):
    slope_tsd, intercept_tsd, r_value_tsd, p_value_tsd, std_err_tsd = stats.linregress(tempMean_sc[i], swr_sc[i])

    plt.subplot(321+i)
    plt.scatter(tempMean_sc[i],swr_sc[i], label = 'dec%.1f  r2=%.6f'%(gridcells_num[i],r_value_tsd**2), s=10**2, color = 'red')

    plt.ylabel('SWR (w/m2)', fontsize=20)
    plt.yticks(fontsize=10)
    plt.xlabel('daily Dec. 2016 temp in %.1f (c)'%(gridcells_num[i]), fontsize=20)
    plt.xticks(fontsize=10)
    plt.legend(fontsize=20, loc = 'upper right')

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/temp_swr_dec.png')

# the line equation:
#print "y=%.6fx+(%.6f)"%(z[0],z[1])
#%%
tempMean_sc = [tempMean_sc_j00,tempMean_sc_j01,tempMean_sc_j02,
               tempMean_sc_j10,tempMean_sc_j11,tempMean_sc_j12]
swr_sc = [swr_sc_j00,swr_sc_j01,swr_sc_j02,
          swr_sc_j10,swr_sc_j11,swr_sc_j12]
gridcells_num = [0.0, 0.1, 0.2, 1.0, 1.1, 1.2]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))
plt.title('Sagehen Creek temp-swr_jan', fontsize=35)
for i in range (6):
    slope_tsj, intercept_tsj, r_value_tsj, p_value_tsj, stj_err_tsj = stats.linregress(tempMean_sc[i], swr_sc[i])

    plt.subplot(321+i)
    plt.scatter(tempMean_sc[i],swr_sc[i], label = 'dec%.1f  r2=%.6f'%(gridcells_num[i],r_value_tsd**2), s=10**2, color = 'red')

    plt.ylabel('SWR (w/m2)', fontsize=20)
    plt.yticks(fontsize=10)
    plt.xlabel('daily Jan. 2016 temp in %.1f (c)'%(gridcells_num[i]), fontsize=20)
    plt.xticks(fontsize=10)
    plt.legend(fontsize=20, loc = 'upper right')

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/temp_swr_jan.png')
#%%

tempMean_sc = [tempMean_sc_d00,tempMean_sc_d01,tempMean_sc_d02,
               tempMean_sc_d10,tempMean_sc_d11,tempMean_sc_d12]
swr_sc = [swr_sc_d00,swr_sc_d01,swr_sc_d02,
          swr_sc_d10,swr_sc_d11,swr_sc_d12]
gridcells_num = [0.0, 0.1, 0.2, 1.0, 1.1, 1.2]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))
plt.title('Sagehen Creek temp-swr_dec', fontsize=35)
for i in range (6):
    slope_tsd, intercept_tsd, r_value_tsd, p_value_tsd, std_err_tsd = stats.linregress(tempMean_sc[i], swr_sc[i])

    plt.subplot(321+i)
    plt.scatter(tempMean_sc[i],swr_sc[i], label = 'dec%.1f  r2=%.6f'%(gridcells_num[i],r_value_tsd**2), s=10**2, color = 'red')

    plt.ylabel('SWR (w/m2)', fontsize=20)
    plt.yticks(fontsize=10)
    plt.xlabel('daily Feb. 2016 temp in %.1f (c)'%(gridcells_num[i]), fontsize=20)
    plt.xticks(fontsize=10)
    plt.legend(fontsize=20, loc = 'upper right')

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/temp_swr_feb.png')

#%%
