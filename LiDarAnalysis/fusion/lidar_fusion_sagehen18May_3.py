import laspy as ls
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import csv
import os
import rasterio

#%% functions

def readPlotDEM(filename,elevationMissNo,pathName):#,deletRows
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
    extent=[x0, x1, y1, y0]
    
    plt.figure(figsize=(30,20))
    plt.imshow(elevation, cmap='gist_earth', extent=extent)
    plt.savefig(pathName)
    
    return elevation, extent

def readCropPlotASCII(filename,elevationMissNo,pathName,indxMin,indxMax,colMin,colMax):#,deletRows
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    
    #removing first n rows
    elevation_df = pd.DataFrame(elevation)
    elevation_df.drop(elevation_df.index[:indxMin], inplace=True) #veg2000,snow1534
    elevation_df.drop(elevation_df.index[indxMax:], inplace=True) #7000,5500
    elevation_df.drop(elevation_df.columns[:colMin], axis = 1, inplace = True) #1000,1593
    elevation_df.drop(elevation_df.columns[colMax:], axis = 1, inplace = True) #0,6508
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation_df.values.shape
    X0 = x0+colMin
    Y0 = y0-indxMin
    x1 = X0 + dx * ncols
    y1 = Y0 + dy * nrows
    extent=[X0, x1, y1, Y0]
    
    plt.figure(figsize=(30,20))
    plt.imshow(elevation_df.values, cmap='gist_earth', extent=extent)
    plt.savefig(pathName)
    
    return elevation_df.values

def creatingCentroidGroundpointsFromDem(tiffFilename,elevationMissNo):#,pathNameforDemImage):
    demset = gdal.Open(tiffFilename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    elevation[elevation > 10000.] = elevationMissNo

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
    elevation_rp = np.reshape(elevation,(nrows*ncols)).T
    dem_groundPoints = np.vstack([latitude_rp,longitude_rp,elevation_rp]).T
    
    dem_groundPoints_df0 = pd.DataFrame(dem_groundPoints,columns=['x','y','z'])
    dem_groundPoints_df = pd.concat([dem_groundPoints_df0[['x','y']].astype(int),dem_groundPoints_df0['z']], axis=1)
    dem_groundPoints_df.sort_values(by=['x','y'],inplace=True)
    dem_groundPoints_df.index=np.arange(0,len(dem_groundPoints_df))
    
    return dem_groundPoints,dem_groundPoints_df

def creatingLatituteLongitudeFromDem(tiffFilename,elevationMissNo):#,pathNameforDemImage):
    demset = gdal.Open(tiffFilename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 10000.] = elevationMissNo

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
    
    return latitude_rp, longitude_rp, nrows, ncols

def creatingCentroidGroundpointsFromASCII(tiffFilename,elevationMissNo,indxMin,indxMax,colMin,colMax):#,pathNameforDemImage):
    demset = gdal.Open(tiffFilename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    elevation[elevation > 10000.] = elevationMissNo
    
    elevation_df = pd.DataFrame(elevation)
    elevation_df.drop(elevation_df.index[:indxMin], inplace=True) #veg2000,snow1534
    elevation_df.drop(elevation_df.index[indxMax:], inplace=True) #7000,5500
    elevation_df.drop(elevation_df.columns[:colMin], axis = 1, inplace = True) #1000,1593
    elevation_df.drop(elevation_df.columns[colMax:], axis = 1, inplace = True) #0,6508
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    X0 = x0+colMin
    Y0 = y0-indxMin
    
    nrows, ncols = elevation_df.values.shape
        
    latitude =[]
    for x in range (ncols):
        latitude.append(x+X0)
    longitude = []
    for y in range (nrows):
        longitude.append(Y0-y)
    
    latitude_rp = (np.tile(latitude, nrows)).astype(int)
    longitude_rp = (np.repeat(longitude, ncols)).astype(int)
    elevation_rp = np.reshape(elevation_df.values,(nrows*ncols)).T
    ascii_groundPoints = np.vstack([latitude_rp,longitude_rp,elevation_rp]).T
    
    ascii_groundPoints_df0 = pd.DataFrame(ascii_groundPoints,columns=['x','y','z'])
    ascii_groundPoints_df = pd.concat([ascii_groundPoints_df0[['x','y']].astype(int),ascii_groundPoints_df0['z']], axis=1)
    ascii_groundPoints_df.sort_values(by=['x','y'],inplace=True)
    ascii_groundPoints_df.index=np.arange(0,len(ascii_groundPoints_df))
    
    return ascii_groundPoints,ascii_groundPoints_df

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

def cutPart0fMap(arry,minX,maxX,minY,maxY):
    arry_df = pd.DataFrame(arry)
    points_df_int = arry_df[(arry_df[0] >= minX) & (arry_df[0] <= maxX)]
    points_df_int_sp0 = points_df_int[(points_df_int[1] >= minY) & (points_df_int[1] <= maxY)]
    points_df_sp_intg = pd.concat([points_df_int_sp0], axis=1)
    points_df_sp_intg.index=np.arange(0,len(points_df_sp_intg))
    
    return points_df_sp_intg

def calculatingSlopeAspectNorthness(demFilename,slopeFileName,aspectFileName,elevation,nrows,ncols,lat,Long):
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
    
    northness_col = np.reshape(northness,(nrows*ncols)).T
    elevation_col = np.reshape(elevation,(nrows*ncols)).T
    slope_col = np.reshape(slope,(nrows*ncols)).T
    elevNrth_gp = np.vstack([lat,Long,northness_col,elevation_col]).T
    index_ls = np.vstack([lat,Long,northness_col,elevation_col,slope_col,elevation_col]).T
    return elevNrth_gp, index_ls

#classification with northeness and elevation index
def classificationNorthnessElevASCII(index_ls,nrth1,nrth2,elev1,elev2):#0.25
    #1 <-0.1 exposed; 1 lowElev
    index_ls[:,6][(index_ls[:,3]<nrth1)&(index_ls[:,2]<elev1)]=11
    index_ls[:,6][(index_ls[:,3]>nrth2)&(index_ls[:,2]<elev1)]=31
    index_ls[:,6][(index_ls[:,3]<=nrth2)&(index_ls[:,3]>=nrth1)&(index_ls[:,2]<elev1)]=21
    index_ls[:,6][(index_ls[:,3]<nrth1)&(index_ls[:,2]>elev2)]=13
    index_ls[:,6][(index_ls[:,3]>nrth2)&(index_ls[:,2]>elev2)]=33
    index_ls[:,6][(index_ls[:,3]<=nrth2)&(index_ls[:,3]>=nrth1)&(index_ls[:,2]>elev2)]=23
    index_ls[:,6][(index_ls[:,3]<nrth1)&(index_ls[:,2]>=elev1)&(index_ls[:,2]<=elev2)]=12
    index_ls[:,6][(index_ls[:,3]>nrth2)&(index_ls[:,2]>=elev1)&(index_ls[:,2]<=elev2)]=32
    index_ls[:,6][(index_ls[:,3]<=nrth2)&(index_ls[:,3]>=nrth1)&(index_ls[:,2]>=elev1)&(index_ls[:,2]<=elev2)]=22

    pointFile_df = pd.DataFrame(index_ls,columns=['x','y','z','north','slope','vegClass','properties'])
    return pointFile_df.values

#%% dem snow off (veg) 
elevationMissNoS0f = 1700.
#G:/hhs_DST/lidar_analysis/sagehen/fusion/
filenameIn = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/dem2014snow0ff.tif"
pathNameOut = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/dem2014snow0ff_G.png"
elevationVeg, extent_dem = readPlotDEM(filenameIn,elevationMissNoS0f,pathNameOut)#

def readPlotDEM(filename,elevationMissNo,pathName):#,deletRows
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
    extent=[x0, x1, y1, y0]
    
    plt.figure(figsize=(30,20))
    plt.imshow(elevation, cmap='gist_earth', extent=extent)
    plt.savefig(pathName)
    
    return elevation, extent

dem_groundPoints, dem_groundPointsVeg_df = creatingCentroidGroundpointsFromDem(filenameIn,elevationMissNoS0f)#,pathNameS)
lat, longt, nrows, ncols = creatingLatituteLongitudeFromDem(filenameIn,elevationMissNoS0f)

#calculating slope, aspect, northnes
elevNrth_gp, index_df = calculatingSlopeAspectNorthness(filenameIn,
                                                        'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/slopeSc.tif',
                                                        'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/aspectSc.tif',
                                                         elevationVeg,nrows,ncols,lat,longt)

index_df[:,4][(index_df[:,4]==-9999.)]=31
index_colNm_df = pd.DataFrame(index_df, columns = ['x','y','North','z','slope','z'])
#index_noG30 = index_df[(index_df[:,4]<30.)] #slope less than 30

dim_dem = [np.min(dem_groundPointsVeg_df['x']),np.max(dem_groundPointsVeg_df['x']),np.min(dem_groundPointsVeg_df['y']),np.max(dem_groundPointsVeg_df['y'])]
#%% ascii file snow 0ff
#allcount ==== allcovercount; totalall ==== totalcount
#loading veg ascii file (count)
veg_count_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_count.asc"
veg_count_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_count.png"
elevationVeg_count, extent_countVeg = readPlotDEM(veg_count_path,-4,veg_count_path_out)#,2000,7000,1000,11000
#veg_count_lst, veg_count_df = creatingCentroidGroundpointsFromASCII(veg_count_path,-4,2000,7000,1000,11000)

#reading veg ascii file (totalCount)
veg_totalCount_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_totalCount.asc"
veg_totalCount_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_totalCount.png"
elevationVeg_totalCount, extent_totalCountVeg = readPlotDEM(veg_totalCount_path,0,veg_totalCount_path_out)#,2000,7000,1000,11000)#
#veg_totalCount_lst, veg_totalCount_df = creatingCentroidGroundpointsFromASCII(veg_totalCount_path,0,2000,7000,1000,11000)

#reading veg ascii file (allCoverCount)
veg_allCoverCount_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_allCoverCount.asc"
veg_allCoverCount_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_allCoverCount.png"
elevationVeg_allCoverCount, extent_allCoverCountVeg = readPlotDEM(veg_allCoverCount_path,0,veg_allCoverCount_path_out)#,2000,7000,1000,11000)#
#veg_allCoverCount_lst, veg_allCoverCount_df = creatingCentroidGroundpointsFromASCII(veg_allCoverCount_path,0,2000,7000,1000,11000)

#%% ascii files snow 0n
#loading snow ascii file (count)
snow_count_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_snow_count3.asc"
snow_count_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_snow_count3.png"
elevationsnow_count3, extent_countSnw3 = readPlotDEM(snow_count_path,0,snow_count_path_out)#,1534,5500,1593,6508
#snow_count_lst, snow_count_df = creatingCentroidGroundpointsFromASCII(snow_count_path,0,1534,5500,1593,6508)

#reading snow ascii file (totalCount)
snow_totalCount_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_snow_totalCount3.asc"
snow_totalCount_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_snow_totalCount3.png"
elevationsnow_totalCount3, extent_totalCountSnw3 = readPlotDEM(snow_totalCount_path,0,snow_totalCount_path_out)#,1534,5500,1593,6508
#snow_totalCount_lst, snow_totalCount_df = creatingCentroidGroundpointsFromASCII(snow_totalCount_path,0,1534,5500,1593,6508)

#reading snow ascii file (allCoverCount)
snow_allCoverCount_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_snow_allCoverCount3.asc"
snow_allCoverCount_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_snow_allCoverCount3.png"
elevationsnow_allCoverCount3, extent_allCoverCountSnw3 = readPlotDEM(snow_allCoverCount_path,0,snow_allCoverCount_path_out)#,1534,5500,1593,6508
#snow_allCoverCount_lst, snow_allCoverCount_df = creatingCentroidGroundpointsFromASCII(snow_allCoverCount_path,0,1534,5500,1593,6508)
extent = [730236,738827,4364802,4372274]
extent_dem = [729999.5, 741999.5, 4364800.5, 4374000.5]
extent_countVeg = [728000.0, 740000.0, 4364000.0, 4374000.0]
extent_countSnw = [730235.0, 738828.0, 4364741.0, 4372275.0]
extent_countSnw3 = [728845.0, 740251.0, 4364115.0, 4373021.0]
#%%cutting some part of maps
minX = np.max([np.min(dem_groundPointsVeg_df['x']),np.min(index_colNm_df['x']),
               np.min(veg_allCoverCount_df['x']),np.min(veg_totalCount_df['x']),np.min(veg_count_df['x']),
               np.min(snow_allCoverCount_df['x']),np.min(snow_totalCount_df['x']),np.min(snow_count_df['x'])])#+2000
    
maxX = np.min([np.max(dem_groundPointsVeg_df['x']),np.max(index_colNm_df['x']),
               np.max(veg_allCoverCount_df['x']),np.max(veg_totalCount_df['x']),np.max(veg_count_df['x']),
               np.max(snow_allCoverCount_df['x']),np.max(snow_totalCount_df['x']),np.max(snow_count_df['x'])])#minXVcm+500
    
minY = np.max([np.min(dem_groundPointsVeg_df['y']),np.min(index_colNm_df['y']),
               np.min(veg_allCoverCount_df['y']),np.min(veg_totalCount_df['y']),np.min(veg_count_df['y']),
               np.min(snow_allCoverCount_df['y']),np.min(snow_totalCount_df['y']),np.min(snow_count_df['y'])])#+1500
    
maxY = np.min([np.max(dem_groundPointsVeg_df['y']),np.max(index_colNm_df['y']),
               np.max(veg_allCoverCount_df['y']),np.max(veg_totalCount_df['y']),np.max(veg_count_df['y']),
               np.max(snow_allCoverCount_df['y']),np.max(snow_totalCount_df['y']),np.max(snow_count_df['y'])])#minYVcm+510

#cutting dem file
dem_Veg_intg = cutPart0fMapSort(dem_groundPointsVeg_df,minX,maxX,minY,maxY)
index_intg = cutPart0fASCIISort(index_colNm_df,minX,maxX,minY,maxY)

#cutting veg maps
veg_count_intg = cutPart0fMapSort(veg_count_df,minX,maxX,minY,maxY)
veg_totalCount_intg = cutPart0fMapSort(veg_totalCount_df,minX,maxX,minY,maxY)
veg_allCoverCount_intg = cutPart0fMapSort(veg_allCoverCount_df,minX,maxX,minY,maxY)

#cutting snow maps
snow_count_intg = cutPart0fMapSort(snow_count_df,minX,maxX,minY,maxY)
snow_totalCount_intg = cutPart0fMapSort(snow_totalCount_df,minX,maxX,minY,maxY)
snow_allCoverCount_intg = cutPart0fMapSort(snow_allCoverCount_df,minX,maxX,minY,maxY)

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

aaaa_test = veg_ascii_index[14780:14870,:]
#%%Classification snow 0n #snow classification
#snow_totalCount_intg=all returns %%% snow_count_intg = above 0.15m %%% snow_allCoverCount_intg = above 2m 
snow_ascii_index = np.vstack([dem_Veg_intg['x'],dem_Veg_intg['y'],dem_Veg_intg['z'],
                             index_intg['North'],index_intg['slope'],
                             snow_totalCount_intg['z'],snow_allCoverCount_intg['z'],snow_count_intg['z'],
                             dem_Veg_intg['x']]).T
# N0 retrurn
snow_ascii_index[:,8][(snow_ascii_index[:,5]==0)]=-99
# N0 snow on the ground
snow_ascii_index[:,8][(snow_ascii_index[:,5]>snow_ascii_index[:,6]) & (snow_ascii_index[:,6]==snow_ascii_index[:,7])]=0 #?????
# snow on the Tree
snow_ascii_index[:,8][(snow_ascii_index[:,6]>snow_ascii_index[:,5]) & (snow_ascii_index[:,7]==snow_ascii_index[:,6])]=-9
# snow on the ground
snow_ascii_index[:,8][(snow_ascii_index[:,5]>0) & (snow_ascii_index[:,6]<snow_ascii_index[:,7])]=10

aaaa_test = snow_ascii_index[14780:14870,:]

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

aaaa_test = ascii_index[14780:14870,:]

#%% terrain properties
ascii_index_pr = np.vstack([ascii_index[:,0],ascii_index[:,1],ascii_index[:,2],ascii_index[:,3],ascii_index[:,4],ascii_index[:,7],dem_Veg_intg['x']]).T
newNcols = maxY - minY + 1

#calculating topo_dimension from cut northness
properties = [-0.1,0.1,2100,2400]
ascii_topo_veg_index = classificationNorthnessElevASCII(ascii_index_pr,properties[0],properties[1],properties[2],properties[3])
# slope less than 30
ascii_topo_veg_index_ls = ascii_topo_veg_index[(ascii_topo_veg_index[:,5]>-90)&(ascii_topo_veg_index[:,4]<30)]

aaaa_test = ascii_topo_veg_index_ls[14780:14870,:]
#%%fSCA no topography
ascii_ut_snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==77])
ascii_ut_n0snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==-7])
ascii_fSCA_ut = 100.*ascii_ut_snow/(ascii_ut_snow+ascii_ut_n0snow)

ascii_0p_snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==44])
ascii_0p_n0snow = len(ascii_topo_veg_index_ls[ascii_topo_veg_index_ls[:,5]==-4])
ascii_fSCA_0p = 100.*ascii_0p_snow/(ascii_0p_snow+ascii_0p_n0snow)
#%% all classification based on topography and northness for under canopy snow
# s:sheltered; n:neutral; e:exposed; l:low; m:mid; h:high
#snow under tree
ascii_ut_snow_le = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==11)])
ascii_ut_snow_lf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==21)])
ascii_ut_snow_ls = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==31)])
ascii_ut_snow_me = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==12)])
ascii_ut_snow_mf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==22)])
ascii_ut_snow_ms = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==32)])
ascii_ut_snow_he = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==13)])
ascii_ut_snow_hf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==23)])
ascii_ut_snow_hs = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==77)&(ascii_topo_veg_index_ls[:,6]==33)])
#n0 snow under tree
ascii_ut_n0snow_le = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==11)])
ascii_ut_n0snow_lf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==21)])
ascii_ut_n0snow_ls = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==31)])
ascii_ut_n0snow_me = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==12)])
ascii_ut_n0snow_mf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==22)])
ascii_ut_n0snow_ms = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==32)])
ascii_ut_n0snow_he = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==13)])
ascii_ut_n0snow_hf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==23)])
ascii_ut_n0snow_hs = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-7)&(ascii_topo_veg_index_ls[:,6]==33)])
#fsca under tree
ascii_fSCA_ut_le = 100.*ascii_ut_snow_le/(ascii_ut_snow_le+ascii_ut_n0snow_le)
ascii_fSCA_ut_lf = 100.*ascii_ut_snow_lf/(ascii_ut_snow_lf+ascii_ut_n0snow_lf)
ascii_fSCA_ut_ls = 100.*ascii_ut_snow_ls/(ascii_ut_snow_ls+ascii_ut_n0snow_ls)
ascii_fSCA_ut_me = 100.*ascii_ut_snow_me/(ascii_ut_snow_me+ascii_ut_n0snow_me)
ascii_fSCA_ut_mf = 100.*ascii_ut_snow_mf/(ascii_ut_snow_mf+ascii_ut_n0snow_mf)
ascii_fSCA_ut_ms = 100.*ascii_ut_snow_ms/(ascii_ut_snow_ms+ascii_ut_n0snow_ms)
ascii_fSCA_ut_he = 100.*ascii_ut_snow_he/(ascii_ut_snow_he+ascii_ut_n0snow_he)
ascii_fSCA_ut_hf = 100.*ascii_ut_snow_hf/(ascii_ut_snow_hf+ascii_ut_n0snow_hf)
ascii_fSCA_ut_hs = 100.*ascii_ut_snow_hs/(ascii_ut_snow_hs+ascii_ut_n0snow_hs)

#snow in 0pen
ascii_0p_snow_le = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==11)])
ascii_0p_snow_lf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==21)])
ascii_0p_snow_ls = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==31)])
ascii_0p_snow_me = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==12)])
ascii_0p_snow_mf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==22)])
ascii_0p_snow_ms = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==32)])
ascii_0p_snow_he = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==13)])
ascii_0p_snow_hf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==23)])
ascii_0p_snow_hs = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==44)&(ascii_topo_veg_index_ls[:,6]==33)])
#no snow in 0pen
ascii_0p_n0snow_le = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==11)])
ascii_0p_n0snow_lf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==21)])
ascii_0p_n0snow_ls = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==31)])
ascii_0p_n0snow_me = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==12)])
ascii_0p_n0snow_mf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==22)])
ascii_0p_n0snow_ms = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==32)])
ascii_0p_n0snow_he = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==13)])
ascii_0p_n0snow_hf = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==23)])
ascii_0p_n0snow_hs = len(ascii_topo_veg_index_ls[(ascii_topo_veg_index_ls[:,5]==-4)&(ascii_topo_veg_index_ls[:,6]==33)])
#fsca under tree
ascii_fSCA_0p_le = 100.*ascii_0p_snow_le/(ascii_0p_snow_le+ascii_0p_n0snow_le)
ascii_fSCA_0p_lf = 100.*ascii_0p_snow_lf/(ascii_0p_snow_lf+ascii_0p_n0snow_lf)
ascii_fSCA_0p_ls = 100.*ascii_0p_snow_ls/(ascii_0p_snow_ls+ascii_0p_n0snow_ls)
ascii_fSCA_0p_me = 100.*ascii_0p_snow_me/(ascii_0p_snow_me+ascii_0p_n0snow_me)
ascii_fSCA_0p_mf = 100.*ascii_0p_snow_mf/(ascii_0p_snow_mf+ascii_0p_n0snow_mf)
ascii_fSCA_0p_ms = 100.*ascii_0p_snow_ms/(ascii_0p_snow_ms+ascii_0p_n0snow_ms)
ascii_fSCA_0p_he = 100.*ascii_0p_snow_he/(ascii_0p_snow_he+ascii_0p_n0snow_he)
ascii_fSCA_0p_hf = 100.*ascii_0p_snow_hf/(ascii_0p_snow_hf+ascii_0p_n0snow_hf)
ascii_fSCA_0p_hs = 100.*ascii_0p_snow_hs/(ascii_0p_snow_hs+ascii_0p_n0snow_hs)
#%% ploting
n_groups = 9
x11 = ['exp_low_0p','flt_low_0p','shl_low_0p','exp_mid_0p',
      'flt_mid_0p','shl_mid_0p','exp_hi_0p','flt_hi_0p','shl_hi_0p']
d11 = [ascii_fSCA_0p_le,ascii_fSCA_0p_lf,ascii_fSCA_0p_ls,ascii_fSCA_0p_me,ascii_fSCA_0p_mf,
       ascii_fSCA_0p_ms,ascii_fSCA_0p_he,ascii_fSCA_0p_hf,ascii_fSCA_0p_hs]

x22 = ['exp_low_UT','flt_low_UT','shl_low_UT','exp_mid_UT',
      'flt_mid_UT','shl_mid_UT','exp_hi_UT','flt_hi_UT','shl_hi_UT']
d22 = [ascii_fSCA_ut_le,ascii_fSCA_ut_lf,ascii_fSCA_ut_ls,ascii_fSCA_ut_me,ascii_fSCA_ut_mf,ascii_fSCA_ut_ms,ascii_fSCA_ut_he,ascii_fSCA_ut_hf,ascii_fSCA_ut_hs]
# create plot
fig, ax = plt.subplots(figsize=(20,15))
index = np.arange(n_groups)
bar_width = 0.4
opacity = 0.8

rects1 = plt.bar(index, d11, bar_width, color='tan', label='0pen') #edgecolor='black',

rects2 = plt.bar(index + bar_width, d22, bar_width, color='green', label='underTree')

plt.ylabel('fSCA(%)', fontsize=30)
plt.yticks(fontsize=20)

plt.xticks(index + bar_width/2.,('exp_low','flt_low','shl_low','exp_mid','flt_mid','shl_mid','exp_hi','flt_hi','shl_hi'), fontsize=20)

plt.title('fSCA% in Sagehen Creek 17April2016; under canopy vs. open, exposed vs. sheltered', fontsize=20)
plt.legend(fontsize=20, loc = 'upper left')
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/fSCA_pixel_sagehen_elev_topo10_noG30A18_2.png')



