import laspy as ls
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import csv
import os
import rasterio

test_flag=True

#%% functions
elevationMissNoS0f = 1900.

def readPlotDEM(filename,elevationMissNo):#,pathName
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
 
#    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
#    nrows, ncols = elevation.shape
#    x1 = x0 + dx * ncols
#    y1 = y0 + dy * nrows
#    extent=[x0, x1, y1, y0]
    
#    plt.figure(figsize=(30,20))
#    plt.imshow(elevation, cmap='gist_earth', extent=extent)
#    plt.savefig(pathName)
    
    return elevation

def readDEMt0findBoundray(filename,elevationMissNo,pathName):
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
 
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows
    extent=[x0, x1, y1, y0]
   
    return extent

def creatingCentroidGroundpointsFromDem(tiffFilename,elevationMissNo):#,pathNameforDemImage):
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
    
    latitude_rp = np.tile(latitude, nrows)
    longitude_rp = np.repeat(longitude, ncols)
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
    
    latitude_rp = np.tile(latitude, nrows)
    longitude_rp = np.repeat(longitude, ncols)
    
    return latitude_rp, longitude_rp, nrows, ncols

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

def readLasFile (lasFilePath):
    infile = ls.file.File(lasFilePath, mode="r")
    coords = np.vstack((infile.x, infile.y, infile.z)).T
    coords_df0 = pd.DataFrame(coords,columns=['x','y','z'])
    coords_df = pd.concat([coords_df0[['x','y']].astype(int),coords_df0['z']], axis=1)
    coords_df.sort_values(by=['x','y'],inplace=True)
    coords_df.index=np.arange(0,len(coords_df)) 
    
    return coords_df

def calculatingUpGroundPoints (coords_df,dem_groundPoints_df,ncols):
    print("please wait ...")
    coords_bnry = np.vstack([coords_df['x']-np.min(coords_df['x']),coords_df['y']-np.min(coords_df['y']),coords_df['z']]).T
    coords_bnry_df = pd.DataFrame(coords_bnry,columns=['x','y','z'])
    N=len(coords_bnry_df)
    M=len(dem_groundPoints_df)
    print("total length: %s" %N)
    upGroundPoints = [0]*N
    nn=0
    min_x=np.min(coords_df['x'])
    min_y=np.min(coords_df['y'])
    z=dem_groundPoints_df['z']
    for row in coords_bnry_df.itertuples():
        if nn%100000==0:
            print ("processed %s" %nn)
        
        indx = row.y + row.x*(ncols)
        
        if indx<M:
            upG = row.z-z[indx] # dem_groundPoints_df['z'][indx]
        else:
            upG = row.z-z[M-1]
        
        upGroundPoints[nn]=[row.x+min_x,row.y+min_y,upG]
        nn+=1
    return upGroundPoints

def createClassIndexBasedDF(pointFile,ncol):
    pointFile_df0 = pd.DataFrame(pointFile,columns=['x','y','z'])
    pointFile_df = pd.concat([pointFile_df0[['x','y']].astype(int),pointFile_df0['z']], axis=1)
    minX = np.min(pointFile_df['x'])
    minY = np.min(pointFile_df['y'])
    pointFile_indx_DF = np.vstack([pointFile_df['x']-minX,pointFile_df['y']-minY,pointFile_df['z']]).T
    pointFile_indx = pointFile_indx_DF[:,1]+ncol*pointFile_indx_DF[:,0]
    pointFile_df.index = pointFile_indx.astype(int)
    return pointFile_indx, pointFile_df
    
def defineSpecificClassGreater(classesG, specific0bjectHeightG):
    specificClassG = []
    for row in classesG.itertuples():
        if row.z>specific0bjectHeightG:
            specificClassG.append(row)
    return specificClassG

def defineSpecificClassLess (classesL, specific0bjectHeightL):
    specificClassL = []
    for row in classesL.itertuples():
        if row.z<specific0bjectHeightL:
            specificClassL.append(row)
    return specificClassL

def defineLowVegClass2(classes):
    lowVegClass = []
    for row in classes.itertuples():
        if row.z<2 and row.z>0.15:
            lowVegClass.append(row)
    return lowVegClass

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
    
    sheltered_lowElev = index_ls[index_ls[:,4]==11]
    sheltered_lowElev_df = pd.DataFrame(sheltered_lowElev,columns = ['x','y','n','e','idx','Index'])
    nuetral_lowElev = index_ls[index_ls[:,4]==21]
    nuetral_lowElev_df = pd.DataFrame(nuetral_lowElev,columns = ['x','y','n','e','idx','Index'])
    exposed_lowElev = index_ls[index_ls[:,4]==31]
    exposed_lowElev_df = pd.DataFrame(exposed_lowElev,columns = ['x','y','n','e','idx','Index'])
    sheltered_midElev = index_ls[index_ls[:,4]==12]
    sheltered_midElev_df = pd.DataFrame(sheltered_midElev,columns = ['x','y','n','e','idx','Index'])
    nuetral_midElev = index_ls[index_ls[:,4]==22]
    nuetral_midElev_df = pd.DataFrame(nuetral_midElev,columns = ['x','y','n','e','idx','Index'])
    exposed_midElev = index_ls[index_ls[:,4]==32]
    exposed_midElev_df = pd.DataFrame(exposed_midElev,columns = ['x','y','n','e','idx','Index'])
    sheltered_hiElev = index_ls[index_ls[:,4]==13]
    sheltered_hiElev_df = pd.DataFrame(sheltered_hiElev,columns = ['x','y','n','e','idx','Index'])
    nuetral_hiElev = index_ls[index_ls[:,4]==23]
    nuetral_hiElev_df = pd.DataFrame(nuetral_hiElev,columns = ['x','y','n','e','idx','Index'])
    exposed_hiElev = index_ls[index_ls[:,4]==33]
    exposed_hiElev_df = pd.DataFrame(exposed_hiElev,columns = ['x','y','n','e','idx','Index'])
    
    topo_dimension = {'sl':sheltered_lowElev_df,'nl':nuetral_lowElev_df,'el':exposed_lowElev_df,
                      'sm':sheltered_midElev_df,'nm':nuetral_midElev_df,'em':exposed_midElev_df,
                      'sh':sheltered_hiElev_df,'nh':nuetral_hiElev_df,'eh':exposed_hiElev_df}
    return topo_dimension

#%% dem snow off (veg) for vcm
#if test_flag is False:    
filenameS0fvcm = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/gis/dem2014snow0ff.tif"
elevationVegvcm = readPlotDEM(filenameS0fvcm,elevationMissNoS0f)#,pathNameS0fvcm

dem_groundPointsVegVcm, dem_groundPointsVegVcm_df = creatingCentroidGroundpointsFromDem(filenameS0fvcm,elevationMissNoS0f)#,pathNameS)
lat_vcm, long_vcm, nrows_vcm, ncols_vcm = creatingLatituteLongitudeFromDem(filenameS0fvcm,elevationMissNoS0f)

#calculating slope, aspect, northnes
elevNrth_gp_vcm, index_vcm = calculatingSlopeAspectNorthness(filenameS0fvcm,
                                                             'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/slopeSc.tif',
                                                             'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/aspectSc.tif',
                                                              elevationVegvcm,nrows_vcm,ncols_vcm,lat_vcm,long_vcm)
index_vcm[:,4][(index_vcm[:,4]==-9999.)]=31
index_noG30 = index_vcm[(index_vcm[:,4]<30.01)] #slope less than 30

filename_slope = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/slopeSc.tif"
slope_degree = readPlotDEM(filename_slope,0)

#las file snow off (snw) for vcm
lasFilePath_vegVcm1 = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/sagehen/snw0ff_sc_points.las"
coordsVegVcm_df = readLasFile (lasFilePath_vegVcm1)

#las file snow on (snw) for vcm
#from LiDar_DemBasedClassification_scSnow0n import coordsnow0nSc18Y
from LiDar_DemBasedClassification_scSnow0n import coordsnow0nSc18Y
#from LiDar_DemBasedClassification_scSnow0n import coordsnow0nSc17A
coordsSnwVcm_df = coordsnow0nSc18Y.copy()
coordsSnwVcm_df['z'] = coordsSnwVcm_df['z']-0.38

#cutting some part of dem map
minXVcm = np.max([np.min(dem_groundPointsVegVcm_df['x']),np.min(coordsSnwVcm_df['x']),np.min(coordsVegVcm_df['x'])]) #361000-1500
maxXVcm = np.min([np.max(dem_groundPointsVegVcm_df['x']),np.max(coordsSnwVcm_df['x']),np.max(coordsVegVcm_df['x'])])#minXVcm+2000
#maxXVcm = minXVcm + 500.
minYVcm = np.max([np.min(dem_groundPointsVegVcm_df['y']),np.min(coordsSnwVcm_df['y']),np.min(coordsVegVcm_df['y'])])#3971000-1500
maxYVcm = np.min([np.max(dem_groundPointsVegVcm_df['y']),np.max(coordsSnwVcm_df['y']),np.max(coordsVegVcm_df['y'])])#minYVcm+2000
#maxYVcm = minYVcm + 600.

dem_groundPointsVegVcm_intg = cutPart0fMapSort(dem_groundPointsVegVcm_df,minXVcm,maxXVcm,minYVcm,maxYVcm)
coordsVegVcm_intg = cutPart0fMapSort(coordsVegVcm_df,minXVcm,maxXVcm,minYVcm,maxYVcm)
coordsSnwVcm_intg = cutPart0fMapSort(coordsSnwVcm_df,minXVcm,maxXVcm,minYVcm,maxYVcm)
indexVcm_cut = cutPart0fMap(index_noG30,minXVcm,maxXVcm,minYVcm,maxYVcm)

newNcols = maxYVcm - minYVcm + 1
#calculating topo_dimension from cut northness
params = [-0.15,0.15,2020,2120]
topo_dimension_vcm = classificationNorthnessElev(indexVcm_cut.values,params[0],params[1],params[2],params[3],newNcols,minXVcm,minYVcm)

#elevLim = 2650
#higher_elevLim_index = dem_groundPointsVegVcm_intg[dem_groundPointsVegVcm_intg['z']>elevLim].index

#%% classification with grountpoints from veg dem file :VCM
upGroundPointsVegVcm = calculatingUpGroundPoints (coordsVegVcm_intg,dem_groundPointsVegVcm_intg,newNcols)
#classification of snow las file with veg_grountpoints
upGroundPointsSnwVegVcm = calculatingUpGroundPoints (coordsSnwVcm_intg,dem_groundPointsVegVcm_intg,newNcols)

path2npoutVeg="G:/hhs/lidar_analysis/sagehen/upGroundPointsVegSc32.npy"
np.save(path2npoutVeg,upGroundPointsVegVcm)   
upLoadUpGroundPointsVeg=np.load(path2npoutVeg)

path2npoutSnw="G:/hhs/lidar_analysis/sagehen/upGroundPointsSnwSc32.npy"
np.save(path2npoutSnw,upGroundPointsSnwVegVcm)   
upLoadUpGroundPointsSnw=np.load(path2npoutSnw)

#%% vegtation classification from DEM2010 and las2010 snow off
vegClassVcm =  upGroundPointsVegVcm[:]#upLoadUpGroundPointsVeg.copy()
vegClassVcm_indx, vegClassVcm_df = createClassIndexBasedDF(vegClassVcm,newNcols)

#all tree classification based on veg dem
allTreeReturnVcm_df = pd.DataFrame(defineSpecificClassGreater (vegClassVcm_df, 2))
allTreeReturnVcm_indx = pd.DataFrame((allTreeReturnVcm_df['Index'].astype(int)).drop_duplicates(keep="first", inplace=False))

#finding pixels that have trees
allTreeClassVcm_df = vegClassVcm_df.loc[allTreeReturnVcm_indx['Index']]
path3npoutSnw="G:/hhs/lidar_analysis/sagehen/allTreeClassSc3.npy"
np.save(path3npoutSnw,allTreeClassVcm_df)   
   
# trees with low branches
lowVegTreeClassVcm = pd.DataFrame(defineLowVegClass2(allTreeClassVcm_df))
lowVegTreeClassVcm_indx = pd.DataFrame((lowVegTreeClassVcm['Index'].astype(int)).drop_duplicates(keep="first", inplace=False))

# trees with no low branches (excluding all pixels that have low blanches, pixels that have a no low branches)
allTreeClassNoLowVegVcm_indx = (pd.concat([allTreeReturnVcm_indx,lowVegTreeClassVcm_indx])).drop_duplicates(keep=False, inplace=False)
allTreeClassNoLowVegVcm_df = allTreeClassVcm_df.loc[allTreeClassNoLowVegVcm_indx['Index']]
path4npoutSnw="G:/hhs/lidar_analysis/sagehen/allTreeClassNoLowVegSc3.npy"
np.save(path4npoutSnw,allTreeClassNoLowVegVcm_df) 

#all UTen classification based on veg dem
all0penClassVcm = (pd.concat([vegClassVcm_df,allTreeClassVcm_df])).drop_duplicates(keep=False, inplace=False)
#all0penClassVcm.set_index(np.arange(0,len(all0penClassVcm)))

# 0pen with low branches
lowVeg0penClassVcm = pd.DataFrame(defineLowVegClass2(all0penClassVcm))
lowVeg0penClassVcm3c = lowVeg0penClassVcm.drop(['Index'], axis=1)#)

# 0pen with no grass (excluding all pixels that have low blanches from open)
all0penClassNoGrass_df = (pd.concat([lowVeg0penClassVcm3c,all0penClassVcm])).drop_duplicates(keep=False, inplace=False)
path5npoutSnw="G:/hhs/lidar_analysis/sagehen/all0penClassNoGrassSc3.npy"
np.save(path5npoutSnw,all0penClassNoGrass_df) 
#all0penClassNoGrass_df.set_index(np.arange(0,len(all0penClassNoGrass_df)))

#%% snow classification from DEM2010 (snow off) and las2010 snow on 
snowClassVcm = upGroundPointsSnwVegVcm[:] #upLoadUpGroundPointsSnw.copy()
snowClassVcm_indx, snowClassVcm_df = createClassIndexBasedDF(snowClassVcm,newNcols)
#snowClassVcm_df1 = snowClassVcm_df.loc[higher_elevLim_index]

# tree snow pixcels that do not have a low branch
allSnowNoLowVegClassVcm_df = snowClassVcm_df.loc[allTreeClassNoLowVegVcm_indx['Index'].values]
path6npoutSnw="G:/hhs/lidar_analysis/sagehen/allSnowNoLowVegClassSc3.npy"
np.save(path6npoutSnw,allSnowNoLowVegClassVcm_df)

# snow under canopy
allsnowUnderTreeClassVcm_df0 = pd.DataFrame(defineSpecificClassLess (allSnowNoLowVegClassVcm_df, 2))
allsnowUnderTreeClassVcm_df0.index = allsnowUnderTreeClassVcm_df0['Index']
allsnowUnderTreeClassVcm_df = pd.DataFrame(defineSpecificClassGreater (allsnowUnderTreeClassVcm_df0[['x','y','z']], 0.15))
allsnowUnderTreeClassVcm_df.index = allsnowUnderTreeClassVcm_df['Index']
path7npoutSnw="G:/hhs/lidar_analysis/sagehen/allsnowUnderTreeClassSc3.npy"
np.save(path7npoutSnw,allsnowUnderTreeClassVcm_df)

allNoSnowUnderTreeClassVcm_df0 = pd.DataFrame(defineSpecificClassLess (allsnowUnderTreeClassVcm_df0[['x','y','z']], 0.15))
allNoSnowUnderTreeClassVcm_df0.index = allNoSnowUnderTreeClassVcm_df0['Index']
allNoSnowUnderTreeClassVcm_df = pd.DataFrame(defineSpecificClassGreater (allNoSnowUnderTreeClassVcm_df0[['x','y','z']], -0.3))
allNoSnowUnderTreeClassVcm_df.index = allNoSnowUnderTreeClassVcm_df['Index']
path8npoutSnw="G:/hhs/lidar_analysis/sagehen/allNoSnowUnderTreeClassSc3.npy"
np.save(path8npoutSnw,allNoSnowUnderTreeClassVcm_df)

#snow on the ground with no Grass
allSnow0penGroundClassVcm_df0 = snowClassVcm_df.loc[all0penClassNoGrass_df.index.drop_duplicates(keep='first')]
#excluding returns with z more than 2m (conflicting pixels)
allSnow0penGroundClassVcm_df1 = pd.DataFrame(defineSpecificClassLess(allSnow0penGroundClassVcm_df0,5)) 
allSnow0penGroundClassVcm_df1.index = allSnow0penGroundClassVcm_df1['Index']
allSnow0penGroundClassVcm_df = pd.DataFrame(defineSpecificClassGreater (allSnow0penGroundClassVcm_df1[['x','y','z']], 0.15))
allSnow0penGroundClassVcm_df.index = allSnow0penGroundClassVcm_df['Index']
path9npoutSnw="G:/hhs/lidar_analysis/sagehen/allSnow0penGroundClassSc3.npy"
np.save(path9npoutSnw,allSnow0penGroundClassVcm_df)

#no snow on the ground
allNoSnow0penClassVcm0 = pd.DataFrame(defineSpecificClassLess (allSnow0penGroundClassVcm_df1[['x','y','z']], 0.15))
allNoSnow0penClassVcm0.index = allNoSnow0penClassVcm0['Index']
allNoSnow0penClassVcm = pd.DataFrame(defineSpecificClassGreater (allNoSnow0penClassVcm0[['x','y','z']], -0.3))
allNoSnow0penClassVcm.index = allNoSnow0penClassVcm['Index']
path10npoutSnw="G:/hhs/lidar_analysis/sagehen/allNoSnow0penClassSc3.npy"
np.save(path10npoutSnw,allNoSnow0penClassVcm)
#%%fSCA no topography
UT_snowPixels = len(allsnowUnderTreeClassVcm_df[['x','y']].drop_duplicates(keep='first'))
UT_NoSnowPixels = len(allNoSnowUnderTreeClassVcm_df[['x','y']].drop_duplicates(keep='first'))
UT_fSCA = 100.*UT_snowPixels/(UT_snowPixels+UT_NoSnowPixels)

OP_snowPixels = len(allSnow0penGroundClassVcm_df[['x','y']].drop_duplicates(keep='first'))
OP_NoSnowPixels = len(allNoSnow0penClassVcm[['x','y']].drop_duplicates(keep='first'))
OP_fSCA = 100.*OP_snowPixels/(OP_NoSnowPixels+OP_snowPixels)
#%% all classification based on topography and northness for under canopy snow
# s:sheltered; n:neutral; e:exposed; l:low; m:mid; h:high
# calculation of fSCA for each of 9 catagories and under canopy and in 0pen
UnTree_el_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
UT_el = len(UnTree_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_el_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
NS_UT_el = len(NS_UnTree_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_el = 100.*(float(UT_el)/(UT_el+NS_UT_el))

UnTree_em_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
UT_em = len(UnTree_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_em_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
NS_UT_em = len(NS_UnTree_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_em = 100.*(float(UT_em)/(UT_em+NS_UT_em))

UnTree_eh_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
UT_eh = len(UnTree_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_eh_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
NS_UT_eh = len(NS_UnTree_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_eh = 100.*(float(UT_eh)/(UT_eh+NS_UT_eh))

UnTree_nl_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
UT_nl = len(UnTree_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_nl_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
NS_UT_nl = len(NS_UnTree_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_nl = 100.*(float(UT_nl)/(UT_nl+NS_UT_nl))

UnTree_nm_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
UT_nm = len(UnTree_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_nm_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
NS_UT_nm = len(NS_UnTree_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_nm = 100.*(float(UT_nm)/(UT_nm+NS_UT_nm))

UnTree_nh_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
UT_nh = len(UnTree_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_nh_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
NS_UT_nh = len(NS_UnTree_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_nh = 100.*(float(UT_nh)/(UT_nh+NS_UT_nh))

UnTree_sl_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
UT_sl = len(UnTree_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_sl_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
NS_UT_sl = len(NS_UnTree_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_sl = 100.*(float(UT_sl)/(UT_sl+NS_UT_sl))

UnTree_sm_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
UT_sm = len(UnTree_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_sm_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
NS_UT_sm = len(NS_UnTree_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_sm = 100.*(float(UT_sm)/(UT_sm+NS_UT_sm))

UnTree_sh_vcm_df = allsnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
UT_sh = len(UnTree_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_UnTree_sh_vcm_df = allNoSnowUnderTreeClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
NS_UT_sh = len(NS_UnTree_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_ut_sh = 100.*(float(UT_sh)/(UT_sh+NS_UT_sh))

#%% all classification based on topography and northness for 0pen space snow
#calculation of fSCA for each of 9 catagories for 0pen sites
Open_el_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
Op_el = len(Open_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_el_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['el']['Index'].astype(int)].dropna()
NS_Op_el = len(NS_Open_el_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_el = 100.*(Op_el)/(Op_el+NS_Op_el)

Open_em_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
Op_em = len(Open_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_em_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['em']['Index'].astype(int)].dropna()
NS_Op_em = len(NS_Open_em_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_em = 100.*(Op_em)/(Op_em+NS_Op_em)

Open_eh_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
Op_eh = len(Open_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_eh_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['eh']['Index'].astype(int)].dropna()
NS_Op_eh = len(NS_Open_eh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_eh = 100.*(Op_eh)/(Op_eh+NS_Op_eh)

Open_nl_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
Op_nl = len(Open_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_nl_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['nl']['Index'].astype(int)].dropna()
NS_Op_nl = len(NS_Open_nl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_nl = 100.*(Op_nl)/(Op_nl+NS_Op_nl)

Open_nm_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
Op_nm = len(Open_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_nm_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['nm']['Index'].astype(int)].dropna()
NS_Op_nm = len(NS_Open_nm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_nm = 100.*(Op_nm)/(Op_nm+NS_Op_nm)

Open_nh_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
Op_nh = len(Open_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_nh_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['nh']['Index'].astype(int)].dropna()
NS_Op_nh = len(NS_Open_nh_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_nh = 100.*(Op_nh)/(Op_nh+NS_Op_nh)

Open_sl_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
Op_sl = len(Open_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_sl_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['sl']['Index'].astype(int)].dropna()
NS_Op_sl = len(NS_Open_sl_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_sl = 100.*(Op_sl)/(Op_sl+NS_Op_sl)

Open_sm_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
Op_sm = len(Open_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_sm_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['sm']['Index'].astype(int)].dropna()
NS_Op_sm = len(NS_Open_sm_vcm_df[['x','y']].drop_duplicates(keep='first'))
FSCA_0p_sm = 100.*(Op_sm)/(Op_sm+NS_Op_sm)

Open_sh_vcm_df = allSnow0penGroundClassVcm_df.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
Op_sh = len(Open_sh_vcm_df[['x','y']].drop_duplicates(keep='first'))
NS_Open_sh_vcm_df = allNoSnow0penClassVcm.loc[topo_dimension_vcm['sh']['Index'].astype(int)].dropna()
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

plt.title('fSCA% in Sagehen Creek; under canopy vs. open, exposed vs. sheltered',fontsize=20)
plt.legend(fontsize=20)
plt.savefig('G:/hhs/lidar_analysis/sagehen/fSCA_pixel_sagehen_elev_topo15_noG30Y18.png')



