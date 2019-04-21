import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import os
import rasterio
#%% functions

def readPlotDEM(filename,elevationMissNo):#,pathNameextent2 = [indx0=4364802,indx1,col0=730002,col1]
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows
    extent1=[x0, x1, y1, y0]
    
    plt.figure(figsize=(30,20))
    plt.imshow(elevation, cmap='gist_earth', extent=extent1)
    #plt.savefig(pathName)
    
    return elevation

def creatingCentroidGroundpointsFromDem(tiffFilename,elevationMissNo):#,pathNameforDemImage):
    demset = gdal.Open(tiffFilename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    elevation[elevation > 10000.] = elevationMissNo

    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape
    #x1 = x0 + dx * ncols
    #y1 = y0 + dy * nrows
    #extent1=[x0, x1, y1, y0]
    
    latitude =[]
    for x in range (ncols):
        latitude.append(x+x0)
    longitude = []
    for y in range (nrows):
        longitude.append(y0-y)
   
    latitude_rp = (np.tile(latitude, nrows)).astype(int)
    longitude_rp = (np.repeat(longitude, ncols)).astype(int)
    elevation_rp = np.reshape(elevation,(nrows*ncols)).T
    dem_groundPoints = np.vstack([latitude_rp,longitude_rp,elevation_rp]).T
    
    dem_groundPoints_df0 = pd.DataFrame(dem_groundPoints,columns=['x','y','z'])
    dem_groundPoints_df = pd.concat([dem_groundPoints_df0[['x','y']].astype(int),dem_groundPoints_df0['z']], axis=1)
    dem_groundPoints_df.sort_values(by=['x','y'],inplace=True)
    dem_groundPoints_df.index=np.arange(0,len(dem_groundPoints_df))
    
    return dem_groundPoints_df

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

def calculatingSlopeAspectNorthness(demFilename,elevationMissNoS0f,slopeFileName,aspectFileName):
    gdal.DEMProcessing(slopeFileName, demFilename, 'slope')
    with rasterio.open(slopeFileName) as datasets1:
        slope=datasets1.read(1) #slope in degree
    slope_rad = slope * 0.0174533 #slope in radian
    slope_sin = np.sin(slope_rad)

    gdal.DEMProcessing(aspectFileName, demFilename, 'aspect')
    with rasterio.open(aspectFileName) as dataseta:
        aspect=dataseta.read(1)
    aspect_rad = aspect * 0.0174533 #radian
    aspect_cos = np.cos(aspect_rad)
    
    northness=aspect_cos*slope_sin
    
    #elevation,nrows,ncols,lat,Long
    demset = gdal.Open(filenameIn)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNoS0f
    elevation[elevation < 0] = elevationMissNoS0f
    elevation[elevation > 10000.] = elevationMissNoS0f

    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape

    latitude =[]
    for x in range (ncols):
        latitude.append(x+x0)
    longitude = []
    for y in range (nrows):
        longitude.append(y0-y)
   
    latitude_rp = (np.tile(latitude, nrows)).astype(int)
    longitude_rp = (np.repeat(longitude, ncols)).astype(int)
    
    northness_col = np.reshape(northness,(nrows*ncols)).T
    elevation_col = np.reshape(elevation,(nrows*ncols)).T
    slope_col = np.reshape(slope,(nrows*ncols)).T
    elevNrth_gp = np.vstack([latitude_rp,longitude_rp,northness_col,elevation_col]).T
    index_ls = np.vstack([latitude_rp,longitude_rp,northness_col,elevation_col,slope_col,elevation_col]).T
    return elevNrth_gp, index_ls

def classificationGridCells(points_df_ls,gridcell):#0.25
 
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[0][0][0])&(points_df_ls[:,0]<gridcell[0][2][0])&
                      (points_df_ls[:,1]<=gridcell[0][0][1])&(points_df_ls[:,0]<gridcell[0][1][1])]=11
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[1][0][0])&(points_df_ls[:,0]<gridcell[1][2][0])&
                      (points_df_ls[:,1]<=gridcell[1][0][1])&(points_df_ls[:,0]<gridcell[1][1][1])]=12
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[2][0][0])&(points_df_ls[:,0]<gridcell[2][2][0])&
                      (points_df_ls[:,1]<=gridcell[2][0][1])&(points_df_ls[:,0]<gridcell[2][1][1])]=13
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[3][0][0])&(points_df_ls[:,0]<gridcell[3][2][0])&
                      (points_df_ls[:,1]<=gridcell[3][0][1])&(points_df_ls[:,0]<gridcell[3][1][1])]=21
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[4][0][0])&(points_df_ls[:,0]<gridcell[4][2][0])&
                      (points_df_ls[:,1]<=gridcell[4][0][1])&(points_df_ls[:,0]<gridcell[4][1][1])]=22
                
    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[5][0][0])&(points_df_ls[:,0]<gridcell[5][2][0])&
                      (points_df_ls[:,1]<=gridcell[5][0][1])&(points_df_ls[:,0]<gridcell[5][1][1])]=23
    
    pointFile_df = pd.DataFrame(points_df_ls,columns=['x','y','z','north','slope','vegClass','gridClass'])
    
    return pointFile_df.values#, pointFile_df
#%% temperature
tempMin2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2016.nc')
tempMax2016 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2016.nc')
tempMin2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmn_2015.nc')
tempMax2015 = Dataset('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tmmx_2015.nc')
tempMin2016_info = tempMin2016.variables
inSpatialRef = tempMin2016['crs'].spatial_ref
sr = osr.SpatialReference(str(inSpatialRef))
res = sr.AutoIdentifyEPSG()
srid = sr.GetAuthorityCode(None)

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
#%% lat long manipulation
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

extent1 = [730000,741000,4362300,4372000]
edgPoints_sc_lat = [730000,733666.67,737333.33,741000]
edgPoints_sc_long = [4362300,4367150,4372000]

product_edg2 = list(itertools.product(edgPoints_sc_lat, edgPoints_sc_long))
sc_edge_m = []
for tupl2 in product_edg2:
    tupl_ls2 = list(tupl2)
    sc_edge_m.append(tupl_ls2)

gridCells_m = [[sc_edge_m[0],sc_edge_m[1],sc_edge_m[4],sc_edge_m[5]],
               [sc_edge_m[1],sc_edge_m[2],sc_edge_m[5],sc_edge_m[6]],
               [sc_edge_m[2],sc_edge_m[3],sc_edge_m[6],sc_edge_m[7]],
               [sc_edge_m[4],sc_edge_m[5],sc_edge_m[8],sc_edge_m[9]],
               [sc_edge_m[5],sc_edge_m[6],sc_edge_m[9],sc_edge_m[10]],
               [sc_edge_m[6],sc_edge_m[7],sc_edge_m[10],sc_edge_m[11]]]   

#%% dem snow off (veg) 
elevationMissNoS0f = 1700.
#E:/chapter2/lidar_analysis/fusion/sagehen_allcount/
extent2 = [731000,738000,4365000,4372000]

filenameIn = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/fusion/dem2014snow0ff.tif"
#pathNameOut = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/fusion/dem2014snow0ff_G2.png"
elevationVeg = readPlotDEM(filenameIn,elevationMissNoS0f)#,pathNameOut
#demset = gdal.Open(filenameIn)
#inLayer = demset.GetRasterBand(1)
#inLayer.GetSpatialRef()
# loading projection
#sr = osr.SpatialReference(str(inSpatialRef))
dem_groundPointsVeg_df = creatingCentroidGroundpointsFromDem(filenameIn,elevationMissNoS0f)#,pathNameS)
dem_groundPointsVeg_df_cut = cutPart0fMapSort(dem_groundPointsVeg_df,extent2[0],extent2[1],extent2[2],extent2[3])

elevNrth_gp, index_df = calculatingSlopeAspectNorthness(filenameIn,0,
                                                        'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/fusion/slopeSc.tif',
                                                        'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/fusion/aspectSc.tif')
index_df[:,4][(index_df[:,4]==-9999.)]=31
index_colNm_df = pd.DataFrame(index_df, columns = ['x','y','North','z','slope','z'])
index_df_cut = cutPart0fASCIISort(index_colNm_df,extent2[0],extent2[1],extent2[2],extent2[3])
#%% ascii file snow 0ff
#allcount ==== allcovercount; totalall ==== totalcount
#loading veg ascii file (count)
veg_count_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_veg_count.asc"
elevationVeg_count = readPlotDEM(veg_count_path,0)#,veg_count_path_out,2000,7000,1000,11000
veg_count_df = creatingCentroidGroundpointsFromDem(veg_count_path,0)
veg_count_df_cut = cutPart0fMapSort(veg_count_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading veg ascii file (totalCount)
veg_totalCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_veg_totalCount.asc"
elevationVeg_totalCount = readPlotDEM(veg_totalCount_path,0)#,veg_totalCount_path_out,2000,7000,1000,11000)#
veg_totalCount_df = creatingCentroidGroundpointsFromDem(veg_totalCount_path,0)
veg_totalCount_df_cut = cutPart0fMapSort(veg_totalCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading veg ascii file (allCoverCount)
veg_allCoverCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_veg_allCoverCount.asc"
elevationVeg_allCoverCount = readPlotDEM(veg_allCoverCount_path,0)#,veg_allCoverCount_path_out,2000,7000,1000,11000)#
veg_allCoverCount_df = creatingCentroidGroundpointsFromDem(veg_allCoverCount_path,0)
veg_allCoverCount_df_cut = cutPart0fMapSort(veg_allCoverCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

#%% ascii files snow 0n 23May2016
#loading snow ascii file (count)
snow_count_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_count26Mar.asc"
#snow_count_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_count17Apr.asc"
#snow_count_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_count18May.asc"

elevationsnow_count = readPlotDEM(snow_count_path,0)#,snow_count_path_out,1534,5500,1593,6508
snow_count_df = creatingCentroidGroundpointsFromDem(snow_count_path,0)
snow_count_df_cut = cutPart0fMapSort(snow_count_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading snow ascii file (totalCount)
snow_totalCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_totalCount26Mar.asc"
#snow_totalCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_totalCount17Apr.asc"
#snow_totalCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_totalCount18May.asc"

elevationsnow_totalCount = readPlotDEM(snow_totalCount_path,0)#,snow_totalCount_path_out,1534,5500,1593,6508
snow_totalCount_df = creatingCentroidGroundpointsFromDem(snow_totalCount_path,0)
snow_totalCount_df_cut = cutPart0fMapSort(snow_totalCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading snow ascii file (allCoverCount)
snow_allCoverCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_allCoverCount26Mar.asc"
#snow_allCoverCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_allCoverCount17Apr.asc"
#snow_allCoverCount_path = "E:/chapter2/lidar_analysis/fusion/sagehen_allcount/sagehen_snow_allCoverCount18May.asc"

elevationsnow_allCoverCount = readPlotDEM(snow_allCoverCount_path,0)#,snow_allCoverCount_path_out,1534,5500,1593,6508
snow_allCoverCount_df = creatingCentroidGroundpointsFromDem(snow_allCoverCount_path,0)
snow_allCoverCount_df_cut = cutPart0fMapSort(snow_allCoverCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

#%%cutting some part of maps
dem_Veg_intg = dem_groundPointsVeg_df_cut[:]
index_intg = index_df_cut[:]

#cutting veg maps
veg_count_intg = veg_count_df_cut[:]
veg_totalCount_intg = veg_totalCount_df_cut[:]
veg_allCoverCount_intg = veg_allCoverCount_df_cut[:]

#cutting snow maps
snow_count_intg = snow_count_df_cut[:]
snow_totalCount_intg = snow_totalCount_df_cut[:]
snow_allCoverCount_intg = snow_allCoverCount_df_cut[:]
#%%Classification snow 0ff #veg classification
#veg_totalCount_intg=all returns %%% veg_count_intg = above 0.15m %%% veg_allCoverCount_intg = above 2m 
veg_ascii_index = np.vstack([dem_Veg_intg['x'],dem_Veg_intg['y'],dem_Veg_intg['z'],
                             index_intg['North'],index_intg['slope'],
                             veg_totalCount_intg['z'],veg_allCoverCount_intg['z'],veg_count_intg['z'],
                             dem_Veg_intg['x']]).T
# Tree with no low branch
veg_ascii_index[:,8][(veg_ascii_index[:,6])>0 & (veg_ascii_index[:,7]==veg_ascii_index[:,6])]=2
# Tree with low branch
veg_ascii_index[:,8][((veg_ascii_index[:,6])>0) & (veg_ascii_index[:,6]<veg_ascii_index[:,7])]=-2
# open with low branch
veg_ascii_index[:,8][((veg_ascii_index[:,6])==0) & (veg_ascii_index[:,7]>0)]=-1
# open with no low branch
veg_ascii_index[:,8][((veg_ascii_index[:,5]> 0) & (veg_ascii_index[:,6]==0) & (veg_ascii_index[:,7]==0))]=1
#undefined pixels
veg_ascii_index[:,8][((veg_ascii_index[:,8])>200)]= -99

aaaa_test = veg_ascii_index[14780:15780,:]
#%%Classification snow 0n #snow classification
#snow_totalCount_intg=all returns %%% snow_count_intg = above 0.15m %%% snow_allCoverCount_intg = above 2m 
snow_ascii_index = np.vstack([dem_Veg_intg['x'],dem_Veg_intg['y'],dem_Veg_intg['z'],
                             index_intg['North'],index_intg['slope'],
                             snow_totalCount_intg['z'],snow_allCoverCount_intg['z'],snow_count_intg['z'],
                             dem_Veg_intg['x']]).T

# N0 snow on the ground
snow_ascii_index[:,8][(snow_ascii_index[:,5]>snow_ascii_index[:,7]) & (snow_ascii_index[:,6]==snow_ascii_index[:,7])]=0 #?????
# snow on the Tree 
snow_ascii_index[:,8][(snow_ascii_index[:,6]==snow_ascii_index[:,5]) & (snow_ascii_index[:,7]==snow_ascii_index[:,6])]=-9
# snow on the ground
snow_ascii_index[:,8][(snow_ascii_index[:,5]>0) & (snow_ascii_index[:,6]<snow_ascii_index[:,7])]=10
# N0 retrurn
snow_ascii_index[:,8][(snow_ascii_index[:,5]==0)]=-99

aaaa_test = snow_ascii_index[34780:35780,:]

#%% final snow-veg classification
ascii_index = np.vstack([dem_Veg_intg['x'],dem_Veg_intg['y'],dem_Veg_intg['z'],
                         index_intg['North'],index_intg['slope'],
                         veg_ascii_index[:,8],snow_ascii_index[:,8],dem_Veg_intg['x']]).T
#snow under tree N0 low branch 77
ascii_index[:,7][(ascii_index[:,5]==2) & (ascii_index[:,6]==10)]=77
#N0 snow under tree -7
ascii_index[:,7][(ascii_index[:,5]==2) & (ascii_index[:,6]==0)]=-7
#snow in 0pen N0 low branch 44
ascii_index[:,7][(ascii_index[:,5]==1) & (ascii_index[:,6]==10)]=44
#N0 snow under tree -7
ascii_index[:,7][(ascii_index[:,5]==1) & (ascii_index[:,6]==0)]=-4
   
ascii_index[:,7][(ascii_index[:,7]>1000)]=-99

aaaa_test = ascii_index[34780:35870,:]

#%% terrain properties
ascii_index_pr = np.vstack([ascii_index[:,0],ascii_index[:,1],ascii_index[:,2],ascii_index[:,3],ascii_index[:,4],ascii_index[:,7],dem_Veg_intg['x']]).T
np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr', ascii_index_pr)

#%% load data
ascii_index_pr = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr.npy')
#calculating topo_dimension from cut northness

ascii_grid_veg_index = classificationGridCells(ascii_index_pr,gridCells_m)
# slope less than 30
ascii_grid_veg_index_ls = ascii_grid_veg_index[(ascii_grid_veg_index[:,5]>-90)&(ascii_grid_veg_index[:,4]<30)]

aaaa_test = ascii_grid_veg_index[604780:605780,:]
#%%fSCA no topography
ascii_ut_snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==77])
ascii_ut_n0snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==-7])
ascii_fSCA_ut = 100.*ascii_ut_snow/(ascii_ut_snow+ascii_ut_n0snow)

ascii_0p_snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==44])
ascii_0p_n0snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==-4])
ascii_fSCA_0p = 100.*ascii_0p_snow/(ascii_0p_snow+ascii_0p_n0snow)
#%%
ascii_ut_snow_le = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==11)])
ascii_ut_snow_lf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==21)])
ascii_ut_snow_ls = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==31)])
ascii_ut_snow_me = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==12)])
ascii_ut_snow_mf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==22)])
ascii_ut_snow_ms = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==32)])



























