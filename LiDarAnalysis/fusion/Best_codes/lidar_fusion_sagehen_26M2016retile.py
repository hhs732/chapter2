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

def readPlotDEM(filename,elevationMissNo,pathName):#extent2 = [indx0=4364802,indx1,col0=730002,col1]
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
    plt.savefig(pathName)
    
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
extent2 = [731000,738000,4365000,4372000]

filenameIn = "G:/hhs_DST/lidar_analysis/sagehen/fusion/dem2014snow0ff.tif"
pathNameOut = "G:/hhs_DST/lidar_analysis/sagehen/fusion/dem2014snow0ff_G.png"
elevationVeg = readPlotDEM(filenameIn,elevationMissNoS0f,pathNameOut)#

dem_groundPointsVeg_df = creatingCentroidGroundpointsFromDem(filenameIn,elevationMissNoS0f)#,pathNameS)
dem_groundPointsVeg_df_cut = cutPart0fMapSort(dem_groundPointsVeg_df,extent2[0],extent2[1],extent2[2],extent2[3])

#calculating slope, aspect, northnes
elevNrth_gp, index_df = calculatingSlopeAspectNorthness(filenameIn,elevationMissNoS0f,
                                                        'G:/hhs_DST/lidar_analysis/sagehen/fusion/slopeSc.tif',
                                                        'G:/hhs_DST/lidar_analysis/sagehen/fusion/aspectSc.tif')
index_df[:,4][(index_df[:,4]==-9999.)]=31
index_colNm_df = pd.DataFrame(index_df, columns = ['x','y','North','z','slope','z'])
index_df_cut = cutPart0fASCIISort(index_colNm_df,extent2[0],extent2[1],extent2[2],extent2[3])

#%% ascii file snow 0ff
#allcount ==== allcovercount; totalall ==== totalcount
#loading veg ascii file (count)
veg_count_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_veg_count.asc"
veg_count_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_veg_count.png"
elevationVeg_count = readPlotDEM(veg_count_path,0,veg_count_path_out)#,2000,7000,1000,11000
veg_count_df = creatingCentroidGroundpointsFromDem(veg_count_path,0)
veg_count_df_cut = cutPart0fMapSort(veg_count_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading veg ascii file (totalCount)
veg_totalCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_veg_totalCount.asc"
veg_totalCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_veg_totalCount.png"
elevationVeg_totalCount = readPlotDEM(veg_totalCount_path,0,veg_totalCount_path_out)#,2000,7000,1000,11000)#
veg_totalCount_df = creatingCentroidGroundpointsFromDem(veg_totalCount_path,0)
veg_totalCount_df_cut = cutPart0fMapSort(veg_totalCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading veg ascii file (allCoverCount)
veg_allCoverCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_veg_allCoverCount.asc"
veg_allCoverCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_veg_allCoverCount.png"
elevationVeg_allCoverCount = readPlotDEM(veg_allCoverCount_path,0,veg_allCoverCount_path_out)#,2000,7000,1000,11000)#
veg_allCoverCount_df = creatingCentroidGroundpointsFromDem(veg_allCoverCount_path,0)
veg_allCoverCount_df_cut = cutPart0fMapSort(veg_allCoverCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

veg_density_sc = veg_allCoverCount_df_cut['z']/veg_totalCount_df_cut['z']
veg_density_sc_ls = np.vstack([veg_allCoverCount_df_cut['x'],veg_allCoverCount_df_cut['y'],veg_density_sc]).T
np.save('G:/hhs_DST/lidar_analysis/fusion/veg_density_sc', veg_density_sc_ls)
#%% ascii files snow 0n 23May2016
#loading snow ascii file (count)
#snow_count_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count26Mar.asc"
#snow_count_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count26Mar.png"
#snow_count_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count17Apr.asc"
#snow_count_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count17Apr.png"
#snow_count_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count18May.asc"
#snow_count_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count18May.png"
snow_count_path = "G:/hhs_DST/lidar_analysis/sagehen/fusionSnow/sagehen_snow_retile_count3.asc"
snow_count_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_count.png"

elevationsnow_count = readPlotDEM(snow_count_path,0,snow_count_path_out)#,1534,5500,1593,6508
snow_count_df = creatingCentroidGroundpointsFromDem(snow_count_path,0)
snow_count_df_cut = cutPart0fMapSort(snow_count_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading snow ascii file (totalCount)
#snow_totalCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount26Mar.asc"
#snow_totalCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount26Mar.png"
#snow_totalCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount17Apr.asc"
#snow_totalCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount17Apr.png"
#snow_totalCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount18May.asc"
#snow_totalCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount18May.png"
snow_totalCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusionSnow/sagehen_snow_retile_totalCount3.asc"
snow_totalCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_totalCount.png"

elevationsnow_totalCount = readPlotDEM(snow_totalCount_path,0,snow_totalCount_path_out)#,1534,5500,1593,6508
snow_totalCount_df = creatingCentroidGroundpointsFromDem(snow_totalCount_path,0)
snow_totalCount_df_cut = cutPart0fMapSort(snow_totalCount_df,extent2[0],extent2[1],extent2[2],extent2[3])

#reading snow ascii file (allCoverCount)
#snow_allCoverCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount26Mar.asc"
#snow_allCoverCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount26Mar.png"
#snow_allCoverCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount17Apr.asc"
#snow_allCoverCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount17Apr.png"
#snow_allCoverCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount18May.asc"
#snow_allCoverCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount18May.png"
snow_allCoverCount_path = "G:/hhs_DST/lidar_analysis/sagehen/fusionSnow/sagehen_snow_retile_allCoverCount3.asc"
snow_allCoverCount_path_out = "G:/hhs_DST/lidar_analysis/sagehen/fusion/fusion_allcount/sagehen_snow_allCoverCount.png"

elevationsnow_allCoverCount = readPlotDEM(snow_allCoverCount_path,0,snow_allCoverCount_path_out)#,1534,5500,1593,6508
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

aaaa_test = ascii_index[24780:25870,:]

#%% terrain properties
ascii_index_pr = np.vstack([ascii_index[:,0],ascii_index[:,1],ascii_index[:,2],ascii_index[:,3],ascii_index[:,4],ascii_index[:,7],dem_Veg_intg['x']]).T
#newNcols = maxY - minY + 1
np.save('G:/hhs_DST/lidar_analysis/fusion/ascii_index_pr_sc26mar_retile', ascii_index_pr)

#calculating topo_dimension from cut northness
properties = [-0.1,0.1,2100,2400]
ascii_topo_veg_index = classificationNorthnessElevASCII(ascii_index_pr,properties[0],properties[1],properties[2],properties[3])
# slope less than 30
ascii_topo_veg_index_ls = ascii_topo_veg_index[(ascii_topo_veg_index[:,5]>-90)&(ascii_topo_veg_index[:,4]<30)]

aaaa_test = ascii_topo_veg_index[14780:15780,:]
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

rects1 = plt.bar(index, d11, bar_width, color='ghostwhite', edgecolor='brown', label='0pen') #  beige  snow

rects2 = plt.bar(index + bar_width, d22, bar_width, color='green',label='underTree')#  olivedrab  mediumseagreen  edgecolor='darkgreen', 

plt.ylabel('fSCA(%)', fontsize=30)
plt.yticks(fontsize=30)

plt.xlabel('Terrain properties classification', fontsize=30)
plt.xticks(index + bar_width/2.,('exp_low','flt_low','shl_low','exp_mid','flt_mid','shl_mid','exp_hi','flt_hi','shl_hi'), fontsize=27, rotation = 15)

plt.title('fSCA% in SCWC in 26 March 2016', fontsize=40)
plt.legend(fontsize=30, loc = 'upper left')
plt.savefig('G:/hhs_DST/lidar_analysis/sagehen/fusion/fSCA_pixel_sagehen_elev_topo10_noG30March26.png')



