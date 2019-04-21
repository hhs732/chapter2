import laspy as ls
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import csv
import os
import rasterio
from decimal import Decimal, ROUND_DOWN
from astropy.io import ascii

#test_flag=True

#%% functions
elevationMissNoS0f = 1500.

def readPlotDEM(filename,elevationMissNo,pathName):#,deletRows
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    
    #removing first 1000 rows
    #elevation_df0 = pd.DataFrame(elevation)
    #elevation_df = elevation_df0.drop(elevation_df0.index[:2000], inplace=False)
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows
    extent=[x0, x1, y1, y0]
    
    plt.figure(figsize=(30,20))
    plt.imshow(elevation, cmap='gist_earth', extent=extent)
    plt.savefig(pathName)
    
    return elevation

def readCropPlotASCII(filename,elevationMissNo,pathName):#,deletRows
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    
    #removing first n rows
    elevation_df = pd.DataFrame(elevation)
    elevation_df.drop(elevation_df.index[:2000], inplace=True)
    elevation_df.drop(elevation_df.index[7000:], inplace=True)
    elevation_df.drop(elevation_df.columns[:1000], axis = 1, inplace = True)
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation_df.values.shape
    x1 = x0+1000 + dx * ncols
    y1 = y0-2000 + dy * nrows
    extent=[x0+1000, x1, y1, y0-2000]
    
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


def creatingCentroidGroundpointsFromASCII(tiffFilename,elevationMissNo):#,pathNameforDemImage):
    demset = gdal.Open(tiffFilename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation < 0] = elevationMissNo
    elevation[elevation > 10000.] = elevationMissNo
    
    elevation_df = pd.DataFrame(elevation)
    elevation_df.drop(elevation_df.index[:2000], inplace=True)
    elevation_df.drop(elevation_df.index[7000:], inplace=True)
    elevation_df.drop(elevation_df.columns[:1000], axis = 1, inplace = True)
    
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    X0 = x0+1000
    Y0 = y0-2000
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

def cutPart0fMap(arry,minX,maxX,minY,maxY):
    arry_df = pd.DataFrame(arry)
    points_df_int = arry_df[(arry_df[0] >= minX) & (arry_df[0] <= maxX)]
    points_df_int_sp0 = points_df_int[(points_df_int[1] >= minY) & (points_df_int[1] <= maxY)]
    points_df_sp_intg = pd.concat([points_df_int_sp0], axis=1)
    points_df_sp_intg.index=np.arange(0,len(points_df_sp_intg))
    
    return points_df_sp_intg

#    return coords_df

def createClassIndexBasedDF(pointFile,ncol):
    pointFile_df0 = pd.DataFrame(pointFile,columns=['x','y','z'])
    pointFile_df = pd.concat([pointFile_df0[['x','y']].astype(int),pointFile_df0['z']], axis=1)
    minX = np.min(pointFile_df['x'])
    minY = np.min(pointFile_df['y'])
    pointFile_indx_DF = np.vstack([pointFile_df['x']-minX,pointFile_df['y']-minY,pointFile_df['z']]).T
    pointFile_indx = pointFile_indx_DF[:,1]+ncol*pointFile_indx_DF[:,0]
    pointFile_df.index = pointFile_indx.astype(int)
    return pointFile_indx, pointFile_df
    
def differenceBetwee2classes (primaryClass,secondNumClass): #to define nolowVegClass and openClass
    primaryNumClass = list(np.arange(len(primaryClass)))
    cleanNumClass = list(set(primaryNumClass)-set(secondNumClass))
    cleanClass = []
    for indx in cleanNumClass:
        cleanClass.append(primaryClass[indx])
    return cleanClass, cleanNumClass

def calculatingSlopeAspectNorthness(demFilename,slopeFileName,aspectFileName,elevation,nrows,ncols,lat,Long):
    gdal.DEMProcessing(slopeFileName, demFilename, 'slope')
    with rasterio.open(slopeFileName) as datasets1:
        slope=datasets1.read(1) #slope in degree
    slope_rad = slope * 0.0174533 #slope in radian
    slope_sin = np.sin(slope_rad)

    gdal.DEMProcessing(aspectFileName, demFilename, 'aspect')
    with rasterio.open(aspectFileName) as dataseta:
        aspect=dataseta.read(1)
    aspect_rad = aspect * 0.0174533
    aspect_cos = np.cos(aspect_rad)

    northness=aspect_cos*slope_sin
    northness_col = np.reshape(northness,(nrows*ncols)).T
    elevation_col = np.reshape(elevation,(nrows*ncols)).T
    slope_col = np.reshape(slope,(nrows*ncols)).T
    elevNrth_gp = np.vstack([lat,Long,northness_col,elevation_col]).T
    index_ls = np.vstack([lat,Long,northness_col,elevation_col,slope_col,elevation_col]).T
#    index_ls_no30 = index_ls[index_ls[:,4]<30]
    return elevNrth_gp, index_ls

#classification with northeness and elevation index
def classificationNorthnessElev(index_ls,nrth1,nrth2,elev1,elev2,ncol,minX,minY):#0.25
    index_ls[:,4][(index_ls[:,2]<nrth1)&(index_ls[:,3]<elev1)]=11
    index_ls[:,4][(index_ls[:,2]>nrth2)&(index_ls[:,3]<elev1)]=31
    index_ls[:,4][(index_ls[:,2]<=nrth2)&(index_ls[:,2]>=nrth1)&(index_ls[:,3]<elev1)]=21
    index_ls[:,4][(index_ls[:,2]<nrth1)&(index_ls[:,3]>elev2)]=13
    index_ls[:,4][(index_ls[:,2]>nrth2)&(index_ls[:,3]>elev2)]=33
    index_ls[:,4][(index_ls[:,2]<=nrth2)&(index_ls[:,2]>=nrth1)&(index_ls[:,3]>elev2)]=23
    index_ls[:,4][(index_ls[:,2]<nrth1)&(index_ls[:,3]>=elev1)&(index_ls[:,3]<=elev2)]=12
    index_ls[:,4][(index_ls[:,2]>nrth2)&(index_ls[:,3]>=elev1)&(index_ls[:,3]<=elev2)]=32
    index_ls[:,4][(index_ls[:,2]<=nrth2)&(index_ls[:,2]>=nrth1)&(index_ls[:,3]>=elev1)&(index_ls[:,3]<=elev2)]=22

    pointFile_df = pd.DataFrame(index_ls,columns=['x','y','z','idx','Index','extra'])
    indx_DF = np.vstack([pointFile_df['x']-minX,pointFile_df['y']-minY]).T
    pointFile_indx = (indx_DF[:,1]+ncol*indx_DF[:,0]).astype(int)
    index_ls[:,5] = pointFile_indx
    
    exposed_lowElev = index_ls[index_ls[:,4]==11]
    exposed_lowElev_df = pd.DataFrame(exposed_lowElev,columns = ['x','y','n','e','idx','Index'])
    nuetral_lowElev = index_ls[index_ls[:,4]==21]
    nuetral_lowElev_df = pd.DataFrame(nuetral_lowElev,columns = ['x','y','n','e','idx','Index'])
    sheltered_lowElev = index_ls[index_ls[:,4]==31]
    sheltered_lowElev_df = pd.DataFrame(sheltered_lowElev,columns = ['x','y','n','e','idx','Index'])
    exposed_midElev = index_ls[index_ls[:,4]==12]
    exposed_midElev_df = pd.DataFrame(exposed_midElev,columns = ['x','y','n','e','idx','Index'])
    nuetral_midElev = index_ls[index_ls[:,4]==22]
    nuetral_midElev_df = pd.DataFrame(nuetral_midElev,columns = ['x','y','n','e','idx','Index'])
    sheltered_midElev = index_ls[index_ls[:,4]==32]
    sheltered_midElev_df = pd.DataFrame(sheltered_midElev,columns = ['x','y','n','e','idx','Index'])
    exposed_hiElev = index_ls[index_ls[:,4]==13]
    exposed_hiElev_df = pd.DataFrame(exposed_hiElev,columns = ['x','y','n','e','idx','Index'])
    nuetral_hiElev = index_ls[index_ls[:,4]==23]
    nuetral_hiElev_df = pd.DataFrame(nuetral_hiElev,columns = ['x','y','n','e','idx','Index'])
    sheltered_hiElev = index_ls[index_ls[:,4]==33]
    sheltered_hiElev_df = pd.DataFrame(sheltered_hiElev,columns = ['x','y','n','e','idx','Index'])
    
    topo_dimension = {'sl':sheltered_lowElev_df,'nl':nuetral_lowElev_df,'el':exposed_lowElev_df,
                      'sm':sheltered_midElev_df,'nm':nuetral_midElev_df,'em':exposed_midElev_df,
                      'sh':sheltered_hiElev_df,'nh':nuetral_hiElev_df,'eh':exposed_hiElev_df}
    return topo_dimension

#%% dem snow off (veg) for vcm
#if test_flag is False:    
filenameIn = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/dem2014snow0ff.tif"
pathNameOut = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/dem2014snow0ff_G.png"
elevationVeg = readPlotDEM(filenameIn,elevationMissNoS0f,pathNameOut)#

dem_groundPoints, dem_groundPointsVeg_df = creatingCentroidGroundpointsFromDem(filenameIn,elevationMissNoS0f)#,pathNameS)
lat, longt, nrows, ncols = creatingLatituteLongitudeFromDem(filenameIn,elevationMissNoS0f)

#calculating slope, aspect, northnes
elevNrth_gp, index_df = calculatingSlopeAspectNorthness(filenameIn,
                                                             'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/slopeSc.tif',
                                                             'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/aspectSc.tif',
                                                              elevationVeg,nrows,ncols,lat,longt)
index_df[:,4][(index_df[:,4]==-9999.)]=31
index_noG30 = index_df[(index_df[:,4]<30.)] #slope less than 30
index_noG30_df = pd.DataFrame(index_noG30, columns = ['x','y','Nor','z','slp','z'])

#allcount ==== allcovercount, totalall ==== totalcount
#loading veg ascii file (count)
veg_count_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_count.asc"
veg_count_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_count.png"
#data = ascii.read(veg_count_path)
#df_veg_count = pd.read_table(veg_count_path)
#lst_veg_count = np.loadtxt(veg_count_path, skiprows=6)
elevationVeg_count = readCropPlotASCII(veg_count_path,-4,veg_count_path_out)#
veg_count_lst, veg_count_df = creatingCentroidGroundpointsFromASCII(veg_count_path,-4)

#reading veg ascii file (totalCount)
veg_totalCount_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_totalCount.asc"
veg_totalCount_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_totalCount.png"
elevationVeg_totalCount = readCropPlotASCII(veg_totalCount_path,0,veg_totalCount_path_out)#
veg_totalCount_lst, veg_totalCount_df = creatingCentroidGroundpointsFromASCII(veg_totalCount_path,0)

#reading veg ascii file (allCoverCount)
veg_allCoverCount_path = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_allCoverCount.asc"
veg_allCoverCount_path_out = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/fusion/sagehen_veg_allCoverCount.png"
elevationVeg_allCoverCount = readCropPlotASCII(veg_allCoverCount_path,0,veg_allCoverCount_path_out)#
veg_allCoverCount_lst, veg_allCoverCount_df = creatingCentroidGroundpointsFromASCII(veg_allCoverCount_path,0)

#cutting some part of maps
minX = np.max([np.min(dem_groundPointsVeg_df['x']),np.min(veg_allCoverCount_df['x']),np.min(veg_totalCount_df['x']),np.min(veg_count_df['x']),np.min(index_noG30_df['x'])])#+2000
maxX = np.min([np.max(dem_groundPointsVeg_df['x']),np.max(veg_allCoverCount_df['x']),np.max(veg_totalCount_df['x']),np.max(veg_count_df['x']),np.max(index_noG30_df['x'])])#minXVcm+500
minY = np.max([np.min(dem_groundPointsVeg_df['y']),np.min(veg_allCoverCount_df['y']),np.min(veg_totalCount_df['y']),np.min(veg_count_df['y']),np.min(index_noG30_df['y'])])#+1500
maxY = np.min([np.max(dem_groundPointsVeg_df['y']),np.max(veg_allCoverCount_df['y']),np.max(veg_totalCount_df['y']),np.max(veg_count_df['y']),np.max(index_noG30_df['y'])])#minYVcm+510

#extent=[minX, maxX, minY, maxY]
#plt.figure(figsize=(30,20))
#plt.imshow(elevationVeg_count, cmap='gist_earth', extent=extent)
#plt.imshow(elevationVeg_totalCount, cmap='gist_earth', extent=extent)
#plt.imshow(elevationVeg_allCoverCount, cmap='gist_earth', extent=extent)

dem_Veg_intg = cutPart0fMapSort(dem_groundPointsVeg_df,minX,maxX,minY,maxY)
veg_count_intg = cutPart0fMapSort(veg_count_df,minX,maxX,minY,maxY)
veg_totalCount_intg = cutPart0fMapSort(veg_totalCount_df,minX,maxX,minY,maxY)
veg_allCoverCount_intg = cutPart0fMapSort(veg_allCoverCount_df,minX,maxX,minY,maxY)
index_cut = cutPart0fMap(index_df,minX,maxX,minY,maxY)

veg_returnsAll = veg_totalCount_intg.values
veg_returnsG2m = veg_allCoverCount_intg.values
veg_returnsG15 = veg_count_intg.values
veg_ascii_index = np.vstack([dem_Veg_intg['x'],dem_Veg_intg['y'],dem_Veg_intg['z'],dem_Veg_intg['x']]).T

for ii in range(len(veg_ascii_index)):
    veg_ascii_index[ii,3][(veg_returnsAll['z'][ii]==veg_returnsG2m['z'][ii])]=11








newNcols = maxY - minY + 1
#calculating topo_dimension from cut northness
params = [-0.1,0.1,2000,2300]
topo_dimension_vcm = classificationNorthnessElev(index_cut.values,params[0],params[1],params[2],params[3],newNcols,minX,minY)

#%% classification with grountpoints from veg dem file :VCM
upGroundPointsVegVcm = calculatingUpGroundPoints (coordsVegVcm_intg,dem_groundPointsVegVcm_intg,newNcols)
#classification of snow las file with veg_grountpoints
upGroundPointsSnwVegVcm = calculatingUpGroundPoints (coordsSnwVcm_intg,dem_groundPointsVegVcm_intg,newNcols)
#path2npoutVeg="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/upGroundPointsVegSc18MY.npy"
#np.save(path2npoutVeg,upGroundPointsVegVcm)   
#upLoadUpGroundPointsVeg=np.load(path2npoutVeg)
#
#path2npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/upGroundPointsSnwSc18MY.npy"
#np.save(path2npoutSnw,upGroundPointsSnwVegVcm)   
#upLoadUpGroundPointsSnw=np.load(path2npoutSnw)

#%% vegtation classification from DEM2010 and las2010 snow off
vegClassVcm =  upGroundPointsVegVcm[:]#upLoadUpGroundPointsVeg.copy()

vegClassVcm_indx, vegClassVcm_df = createClassIndexBasedDF(vegClassVcm,newNcols)

#fig3 = plt.figure(figsize=(20,15))
#ax3 = Axes3D(fig3)
#ax3.scatter(dem_groundPointsVegVcm_intg['x'][0:2000], dem_groundPointsVegVcm_intg['y'], dem_groundPointsVegVcm_intg['z'])
#ax3.scatter(coordsVegVcm_intg['x'], coordsVegVcm_intg['y'], coordsVegVcm_intg['z'])
#ax3.scatter(coordsSnwVcm_intg['x'], coordsSnwVcm_intg['y'], coordsSnwVcm_intg['z'])
#ax3.legend()
#vegClassVcm_df['z'] = vegClassVcm_df['z']
#vegMatrix = allTreeReturnVcm_df['Index'].astype(int)).drop_duplicates(keep="first", inplace=False))
#all tree classification based on veg dem
allTreeReturnVcm_df = pd.DataFrame(defineSpecificClassGreater (vegClassVcm_df, 2))
allTreeReturnVcm_indx = pd.DataFrame((allTreeReturnVcm_df['Index'].astype(int)).drop_duplicates(keep="first", inplace=False))

#finding pixels that have trees
allTreeClassVcm_df = vegClassVcm_df.loc[allTreeReturnVcm_indx['Index']]
   
# trees with low branches
lowVegTreeClassVcm = pd.DataFrame(defineLowVegClass2(allTreeClassVcm_df))

lowVegTreeClassVcm_indx = pd.DataFrame((lowVegTreeClassVcm['Index'].astype(int)).drop_duplicates(keep="first", inplace=False))

# trees with no low branches (excluding all pixels that have low blanches, pixels that have a no low branches)
allTreeClassNoLowVegVcm_indx = (pd.concat([allTreeReturnVcm_indx,lowVegTreeClassVcm_indx])).drop_duplicates(keep=False, inplace=False)
allTreeClassNoLowVegVcm_df = allTreeClassVcm_df.loc[allTreeClassNoLowVegVcm_indx['Index']]
#path4npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/allTreeClassNoLowVegSc18MY.npy"
#np.save(path4npoutSnw,allTreeClassNoLowVegVcm_df) 

#testing tree with no low branch class
allTreeClassNoLowVeg_test1 = allTreeClassNoLowVegVcm_df[((allTreeClassNoLowVegVcm_df['z']<0.15))]
allTreeClassNoLowVeg_test1_pixle = allTreeClassNoLowVeg_test1[['x','y']].drop_duplicates(keep='first', inplace=False)
allTreeClassNoLowVeg_pixle_indx = np.vstack([allTreeClassNoLowVeg_test1_pixle['x']-minXVcm,allTreeClassNoLowVeg_test1_pixle['y']-minYVcm]).T
pointFile_indx1 = allTreeClassNoLowVeg_pixle_indx[:,1]+newNcols*allTreeClassNoLowVeg_pixle_indx[:,0]

allTreeClassNoLowVeg_test2 = allTreeClassNoLowVegVcm_df[((allTreeClassNoLowVegVcm_df['z']>2))]
allTreeClassNoLowVeg_test2_pixle = allTreeClassNoLowVeg_test2[['x','y']].drop_duplicates(keep='first', inplace=False)
allTreeClassNoLowVeg_pixle_indx2 = np.vstack([allTreeClassNoLowVeg_test2_pixle['x']-minXVcm,allTreeClassNoLowVeg_test2_pixle['y']-minYVcm]).T
pointFile_indx2 = allTreeClassNoLowVeg_pixle_indx2[:,1]+newNcols*allTreeClassNoLowVeg_pixle_indx2[:,0]

set(pointFile_indx1) < set(pointFile_indx2)

#all UTen classification based on veg dem
all0penClassVcm = (pd.concat([vegClassVcm_df,allTreeClassVcm_df])).drop_duplicates(keep=False, inplace=False)
all0penClassVcm_indx = ((pd.DataFrame(all0penClassVcm.index)).astype(int)).drop_duplicates(keep='first', inplace=False)
# 0pen with low branches
lowVeg0penClassVcm = pd.DataFrame(defineLowVegClass2(all0penClassVcm))
lowVeg0penClassVcm.index = lowVeg0penClassVcm['Index']

# 0pen with no grass (excluding all pixels that have low blanches from open)
all0penClassNoGrass_indx = (pd.concat([lowVeg0penClassVcm['Index'],all0penClassVcm_indx[0]])).drop_duplicates(keep=False, inplace=False)
all0penClassNoGrass_df = all0penClassVcm.loc[all0penClassNoGrass_indx]
#path5npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/all0penClassNoGrassSc18MY.npy"
#np.save(path5npoutSnw,all0penClassNoGrass_df) 
plt.hist(all0penClassVcm['z'])
all0penClassNoGrass_df_fail = all0penClassNoGrass_df[(all0penClassNoGrass_df['z']>0.15)]
#%% snow classification from DEM2010 (snow off) and las2010 snow on 
snowClassVcm = upGroundPointsSnwVegVcm[:] #upLoadUpGroundPointsSnw.copy()
snowClassVcm_indx, snowClassVcm_df = createClassIndexBasedDF(snowClassVcm,newNcols)

# tree snow pixcels that do not have a low branch
allSnowNoLowVegClassVcm_df = (snowClassVcm_df.loc[allTreeClassNoLowVegVcm_indx['Index'].values]).dropna()
#path6npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/allSnowNoLowVegClassSc18MY.npy"
#np.save(path6npoutSnw,allSnowNoLowVegClassVcm_df)

# snow under canopy
#groupby_mean = allSnowNoLowVegClassVcm_df.groupby(['x','y']).mean()

allSnowUnderTreeClassVcm_df0 = pd.DataFrame(defineSpecificClassLess (allSnowNoLowVegClassVcm_df, 2))
allSnowUnderTreeClassVcm_df0.index = allSnowUnderTreeClassVcm_df0['Index']
allSnowUnderTreeClassVcm_df = pd.DataFrame(defineSpecificClassGreater (allSnowUnderTreeClassVcm_df0[['x','y','z']], 0.15))
allSnowUnderTreeClassVcm_df.index = allSnowUnderTreeClassVcm_df['Index']

allNoSnowUnderTreeClassVcm_df0 = pd.DataFrame(defineSpecificClassLess (allSnowUnderTreeClassVcm_df0[['x','y','z']], 0.15))
allNoSnowUnderTreeClassVcm_df0.index = allNoSnowUnderTreeClassVcm_df0['Index']
allNoSnowUnderTreeClassVcm_df = pd.DataFrame(defineSpecificClassGreater (allNoSnowUnderTreeClassVcm_df0[['x','y','z']], -0.3))
allNoSnowUnderTreeClassVcm_df.index = allNoSnowUnderTreeClassVcm_df['Index']

#test snow under the trees (finding duplicates)
allSnowUnderTreeClassVcm_df_pix = allSnowUnderTreeClassVcm_df[['x','y']].drop_duplicates(keep='first', inplace=False)
allNoSnowUnderTreeClassVcm_df_pix = allNoSnowUnderTreeClassVcm_df[['x','y']].drop_duplicates(keep='first', inplace=False)
concat2DF = pd.concat([allSnowUnderTreeClassVcm_df_pix,allNoSnowUnderTreeClassVcm_df_pix])
concat2DF["is_duplicate"]= concat2DF.duplicated()
#indexing
SnowUnderTreeDuplicate_indx = (pd.DataFrame(((concat2DF[(concat2DF["is_duplicate"] == True)]).index).values)).drop_duplicates(keep='first', inplace=False)
allSnowUnderTreeClassVcm_indx = (pd.DataFrame((allSnowUnderTreeClassVcm_df.index).values)).drop_duplicates(keep='first', inplace=False)
allNoSnowUnderTreeClassVcm_indx = (pd.DataFrame(allNoSnowUnderTreeClassVcm_df.index.values)).drop_duplicates(keep='first', inplace=False)

#final snow under the tree and no snow under tree pixels
SnowUnderTreePix_indx = (pd.concat([SnowUnderTreeDuplicate_indx,allSnowUnderTreeClassVcm_indx])).drop_duplicates(keep=False, inplace=False) #no duplicate
SnowUnderTreeClassVcm_df = allSnowUnderTreeClassVcm_df.loc[SnowUnderTreePix_indx[0].values] #no duplicate
#path7npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/allsnowUnderTreeClassSc18MY.npy"
#np.save(path7npoutSnw,SnowUnderTreeClassVcm_df)

NoSnowUnderTreePix_indx = (pd.concat([SnowUnderTreeDuplicate_indx,allNoSnowUnderTreeClassVcm_indx])).drop_duplicates(keep=False, inplace=False)
NoSnowUnderTreeClassVcm_df = allNoSnowUnderTreeClassVcm_df.loc[NoSnowUnderTreePix_indx[0].values]
#path8npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/allNoSnowUnderTreeClassSc18MY.npy"
#np.save(path8npoutSnw,NoSnowUnderTreeClassVcm_df)
#%% #snow on the ground with no Grass
allSnow0penGroundClassVcm_df0 = (snowClassVcm_df.loc[all0penClassNoGrass_df.index.drop_duplicates(keep='first')]).dropna()
#excluding returns with z more than 2m (conflicting pixels)
allSnow0penGroundClassVcm_df1 = pd.DataFrame(defineSpecificClassLess(allSnow0penGroundClassVcm_df0,5)) 
allSnow0penGroundClassVcm_df1.index = allSnow0penGroundClassVcm_df1['Index']
allSnow0penGroundClassVcm_df = pd.DataFrame(defineSpecificClassGreater (allSnow0penGroundClassVcm_df1[['x','y','z']], 0.15))
allSnow0penGroundClassVcm_df.index = allSnow0penGroundClassVcm_df['Index']

#no snow on the ground
allNoSnow0penClassVcm0 = pd.DataFrame(defineSpecificClassLess (allSnow0penGroundClassVcm_df1[['x','y','z']], 0.15))
allNoSnow0penClassVcm0.index = allNoSnow0penClassVcm0['Index']
allNoSnow0penClassVcm = pd.DataFrame(defineSpecificClassGreater (allNoSnow0penClassVcm0[['x','y','z']], -0.3))
allNoSnow0penClassVcm.index = allNoSnow0penClassVcm['Index']

#finding duplicates for snow in 0pen sites
allSnow0penClassVcm_df_pix = allSnow0penGroundClassVcm_df[['x','y']].drop_duplicates(keep='first', inplace=False)
allNoSnow0penClassVcm_df_pix = allNoSnow0penClassVcm[['x','y']].drop_duplicates(keep='first', inplace=False)
concat2DF2 = pd.concat([allSnow0penClassVcm_df_pix,allNoSnow0penClassVcm_df_pix])
concat2DF2["is_duplicate"]= concat2DF2.duplicated()
#indexing
Snow0penDuplicate_indx = (pd.DataFrame(((concat2DF2[(concat2DF2["is_duplicate"] == True)]).index).values)).drop_duplicates(keep='first', inplace=False)
allSnow0penClassVcm_indx = (pd.DataFrame((allSnow0penGroundClassVcm_df.index).values)).drop_duplicates(keep='first', inplace=False)
allNoSnow0penClassVcm_indx = (pd.DataFrame(allNoSnow0penClassVcm.index.values)).drop_duplicates(keep='first', inplace=False)

#final snow in 0pen sites and no snow in 0pen sites
Snow0penPix_indx = (pd.concat([Snow0penDuplicate_indx,allSnow0penClassVcm_indx])).drop_duplicates(keep=False, inplace=False)
Snow0penClassVcm_df = allSnow0penGroundClassVcm_df.loc[Snow0penPix_indx[0].values]
#path9npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/allSnow0penGroundClassSc18MY.npy"
#np.save(path9npoutSnw,Snow0penClassVcm_df)

NoSnow0penPix_indx = (pd.concat([Snow0penDuplicate_indx,allNoSnow0penClassVcm_indx])).drop_duplicates(keep=False, inplace=False)
NoSnow0penClassVcm_df = allNoSnow0penClassVcm.loc[NoSnow0penPix_indx[0].values]
#path10npoutSnw="C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/allNoSnow0penClassSc18MY.npy"
#np.save(path10npoutSnw,NoSnow0penClassVcm_df)

#%%fSCA no topography
UT_snowPixels = len(SnowUnderTreeClassVcm_df[['x','y']].drop_duplicates(keep='first'))
UT_NoSnowPixels = len(NoSnowUnderTreeClassVcm_df[['x','y']].drop_duplicates(keep='first'))
UT_fSCA = 100.*UT_snowPixels/(UT_snowPixels+UT_NoSnowPixels)

OP_snowPixels = len(Snow0penClassVcm_df[['x','y']].drop_duplicates(keep='first'))
OP_NoSnowPixels = len(NoSnow0penClassVcm_df[['x','y']].drop_duplicates(keep='first'))
OP_fSCA = 100.*OP_snowPixels/(OP_NoSnowPixels+OP_snowPixels)
#%% all classification based on topography and northness for under canopy snow
# s:sheltered; n:neutral; e:exposed; l:low; m:mid; h:high
# calculation of fSCA for each of 9 catagories and under canopy and in 0pen
UnTree_el_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
UT_el = len(UnTree_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_el_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
NS_UT_el = len(NS_UnTree_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_el = 0#100.*(float(UT_el)/(UT_el+NS_UT_el))

UnTree_em_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
UT_em = len(UnTree_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_em_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
NS_UT_em = len(NS_UnTree_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_em = 100.*(float(UT_em)/(UT_em+NS_UT_em))

UnTree_eh_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
UT_eh = len(UnTree_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_eh_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
NS_UT_eh = len(NS_UnTree_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_eh = 100.*(float(UT_eh)/(UT_eh+NS_UT_eh))

UnTree_nl_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
UT_nl = len(UnTree_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_nl_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
NS_UT_nl = len(NS_UnTree_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_nl = 0#100.*(float(UT_nl)/(UT_nl+NS_UT_nl))

UnTree_nm_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
UT_nm = len(UnTree_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_nm_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
NS_UT_nm = len(NS_UnTree_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_nm = 100.*(float(UT_nm)/(UT_nm+NS_UT_nm))

UnTree_nh_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
UT_nh = len(UnTree_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_nh_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
NS_UT_nh = len(NS_UnTree_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_nh = 100.*(float(UT_nh)/(UT_nh+NS_UT_nh))

UnTree_sl_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
UT_sl = len(UnTree_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_sl_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
NS_UT_sl = len(NS_UnTree_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_sl = 0#100.*(float(UT_sl)/(UT_sl+NS_UT_sl))

UnTree_sm_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
UT_sm = len(UnTree_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_sm_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
NS_UT_sm = len(NS_UnTree_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_sm = 100.*(float(UT_sm)/(UT_sm+NS_UT_sm))

UnTree_sh_vcm_df = SnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
UT_sh = len(UnTree_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_sh_vcm_df = NoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
NS_UT_sh = len(NS_UnTree_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_sh = 100.*(float(UT_sh)/(UT_sh+NS_UT_sh))

#%% all classification based on topography and northness for 0pen space snow
#calculation of fSCA for each of 9 catagories for 0pen sites
Open_el_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
Op_el = len(Open_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_el_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
NS_Op_el = len(NS_Open_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_el = 0#100.*(Op_el)/(Op_el+NS_Op_el)

Open_em_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
Op_em = len(Open_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_em_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
NS_Op_em = len(NS_Open_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_em = 100.*(Op_em)/(Op_em+NS_Op_em)

Open_eh_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
Op_eh = len(Open_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_eh_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
NS_Op_eh = len(NS_Open_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_eh = 100.*(Op_eh)/(Op_eh+NS_Op_eh)

Open_nl_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
Op_nl = len(Open_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_nl_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
NS_Op_nl = len(NS_Open_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_nl = 0#100.*(Op_nl)/(Op_nl+NS_Op_nl)

Open_nm_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
Op_nm = len(Open_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_nm_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
NS_Op_nm = len(NS_Open_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_nm = 100.*(Op_nm)/(Op_nm+NS_Op_nm)

Open_nh_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
Op_nh = len(Open_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_nh_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
NS_Op_nh = len(NS_Open_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_nh = 100.*(Op_nh)/(Op_nh+NS_Op_nh)

Open_sl_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
Op_sl = len(Open_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_sl_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
NS_Op_sl = len(NS_Open_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_sl = 0#100.*(Op_sl)/(Op_sl+NS_Op_sl)

Open_sm_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
Op_sm = len(Open_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_sm_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
NS_Op_sm = len(NS_Open_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_sm = 100.*(Op_sm)/(Op_sm+NS_Op_sm)

Open_sh_vcm_df = Snow0penClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
Op_sh = len(Open_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_sh_vcm_df = NoSnow0penClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
NS_Op_sh = len(NS_Open_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_sh = 100.*(Op_sh)/(Op_sh+NS_Op_sh)

#%% ploting
n_groups = 9
x11 = ['exp_low_0p','flt_low_0p','shl_low_0p','exp_mid_0p',
      'flt_mid_0p','shl_mid_0p','exp_hi_0p','flt_hi_0p','shl_hi_0p']
d11 = [FSCA_0p_el,FSCA_0p_nl,FSCA_0p_sl,FSCA_0p_em,FSCA_0p_nm,FSCA_0p_sm,FSCA_0p_eh,FSCA_0p_nh,FSCA_0p_sh]

x22 = ['exp_low_UT','flt_low_UT','shl_low_UT','exp_mid_UT',
      'flt_mid_UT','shl_mid_UT','exp_hi_UT','flt_hi_UT','shl_hi_UT']
d22 = [FSCA_ut_el,FSCA_ut_nl,FSCA_ut_sl,FSCA_ut_em,FSCA_ut_nm,FSCA_ut_sm,FSCA_ut_eh,FSCA_ut_nh,FSCA_ut_sh]
# create plot
fig, ax = plt.subplots(figsize=(15,10))
index = np.arange(n_groups)
bar_width = 0.4
opacity = 0.8

rects1 = plt.bar(index, d11, bar_width,color='snow',edgecolor='tan',label='0pen')

rects2 = plt.bar(index + bar_width, d22, bar_width, color='green',label='underTree')

plt.ylabel('fSCA(%)', fontsize=30)
plt.yticks(fontsize=20)

plt.xticks(index + bar_width/2.,('exp_low','flt_low','shl_low','exp_mid','flt_mid','shl_mid','exp_hi','flt_hi','shl_hi'), fontsize=20)

plt.title('fSCA% in Sagehen Creek 18MaY2016; under canopy vs. open, exposed vs. sheltered',fontsize=20)
plt.legend(fontsize=20)
plt.savefig('G:/hhs_DST/lidar_analysis/sagehen/fSCA_pixel_sagehen_elev_topo15_noG30MY18_test.png')



