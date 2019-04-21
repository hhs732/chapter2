import laspy as ls
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import csv
import os
import rasterio
#%%example
A= pd.DataFrame(np.array([[2,22,1],[3,23,1],[1,21,0],[4,24,9]]),columns=['x','y','z'])

B= pd.DataFrame(np.array([[1,21,2],[2,22,3],[3,23,3],[10,200,3],[4,24,2]]),columns=['x','y','z'])

C=B.copy()
R,c2=B.shape
indx = []
for i in range(R):
    row=B[['x','y']].iloc[i]
    
    R=np.where(A[['x','y']]==row)[0]
    print (R)
    if len(R)>0:
        if R[0]==R[1]:
            print ("row %i of B is the same as row %i of A" %(i,R[0]))
            C['z'][i]=A['z'][R[0]]
            indx.append(i)
    else:
        C['z'][i]=-1

#print (C)
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

# Grab just the X dimension from the file, and scale it.
def scaled_x_dimension(las_file):
    x_dimension = las_file.X
    scale = las_file.header.scale[0]
    offset = las_file.header.offset[0]
    return(x_dimension*scale + offset)
    
def defineSpecificClassGreater(classesG, specific0bjectHeightG):
    specificClassG = classesG[classesG[:,2]>=specific0bjectHeightG]
    return specificClassG

def defineSpecificClassLess (classesL, specific0bjectHeightL):
    specificClassL = classesL[classesL[:,2]<=specific0bjectHeightL]
    return specificClassL

def defineLowVegClass(classes):
    lowVegClass = classes[(classes[:,2]<2) & (classes[:,2]>0.15)]
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
    elevNrth_gp = np.vstack([lat,Long,northness_col,elevation_col]).T
    index_ls = np.vstack([lat,Long,northness_col,elevation_col,northness_col]).T
    
    return elevNrth_gp, index_ls

#classification with northeness and elevation index
def classificationNorthnessElev(index_ls,elev1,elev2):
    index_ls[:,4][(index_ls[:,2]<-0.33)&(index_ls[:,3]<elev1)]=11
    index_ls[:,4][(index_ls[:,2]>0.33)&(index_ls[:,3]<elev1)]=31
    index_ls[:,4][(index_ls[:,2]<=0.33)&(index_ls[:,2]>=-0.33)&(index_ls[:,3]<elev1)]=21
    index_ls[:,4][(index_ls[:,2]<-0.33)&(index_ls[:,3]>elev2)]=13
    index_ls[:,4][(index_ls[:,2]>0.33)&(index_ls[:,3]>elev2)]=33
    index_ls[:,4][(index_ls[:,2]<=0.33)&(index_ls[:,2]>=-0.33)&(index_ls[:,3]>elev2)]=23
    index_ls[:,4][(index_ls[:,2]<-0.33)&(index_ls[:,3]>=elev1)&(index_ls[:,3]<=elev2)]=12
    index_ls[:,4][(index_ls[:,2]>0.33)&(index_ls[:,3]>=elev1)&(index_ls[:,3]<=elev2)]=32
    index_ls[:,4][(index_ls[:,2]<=0.33)&(index_ls[:,2]>=-0.33)&(index_ls[:,3]>=elev1)&(index_ls[:,3]<=elev2)]=22

    sheltered_lowElev = index_ls[index_ls[:,4]==11]
    #sheltered_lowElev_df = pd.DataFrame(sheltered_lowElev,columns = ['x','y','n','e','idx'])
    nuetral_lowElev = index_ls[index_ls[:,4]==21]
    #nuetral_lowElev_df = pd.DataFrame(nuetral_lowElev,columns = ['x','y','n','e','idx'])
    exposed_lowElev = index_ls[index_ls[:,4]==31]
    #exposed_lowElev_df = pd.DataFrame(exposed_lowElev,columns = ['x','y','n','e','idx'])
    sheltered_midElev = index_ls[index_ls[:,4]==12]
    #sheltered_midElev_df = pd.DataFrame(sheltered_midElev,columns = ['x','y','n','e','idx'])
    nuetral_midElev = index_ls[index_ls[:,4]==22]
    #nuetral_midElev_df = pd.DataFrame(nuetral_midElev,columns = ['x','y','n','e','idx'])
    exposed_midElev = index_ls[index_ls[:,4]==32]
    #exposed_midElev_df = pd.DataFrame(exposed_midElev,columns = ['x','y','n','e','idx'])
    sheltered_hiElev = index_ls[index_ls[:,4]==13]
    #sheltered_hiElev_df = pd.DataFrame(sheltered_hiElev,columns = ['x','y','n','e','idx'])
    nuetral_hiElev = index_ls[index_ls[:,4]==23]
    #nuetral_hiElev_df = pd.DataFrame(nuetral_hiElev,columns = ['x','y','n','e','idx'])
    exposed_hiElev = index_ls[index_ls[:,4]==33]
    #exposed_hiElev_df = pd.DataFrame(exposed_hiElev,columns = ['x','y','n','e','idx'])
    
    topo_dimension = {'sl':sheltered_lowElev,'nl':nuetral_lowElev,'el':exposed_lowElev,
                      'sm':sheltered_midElev,'nm':nuetral_midElev,'em':exposed_midElev,
                      'sh':sheltered_hiElev,'nh':nuetral_hiElev,'eh':exposed_hiElev}
    return topo_dimension

def cutPart0fMap(points_df,minX,maxX,minY,maxY):
    points_df_int = points_df[(points_df['x'] >= minX) & (points_df['x'] <= maxX)]
    points_df_int_sp0 = points_df_int[(points_df_int['y'] >= minY) & (points_df_int['y'] <= maxY)]
    #points_df_int_sp = points_df_int_sp0.values
    points_df_sp_intg = pd.concat([points_df_int_sp0[['x','y']].astype(int),points_df_int_sp0['z']], axis=1)
    points_df_sp_intg.sort_values(by=['x','y'],inplace=True)
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

#def calculatingUpGroundPoints (coords_df,dem_groundPoints_df):
#    upGroundPoints = []
#    for row in coords_df.itertuples():
#        indexFinl = np.where(((dem_groundPoints_df['x']==row.x) & (dem_groundPoints_df['y']==row.y)))[0]#
#        if (np.size(indexFinl))>0:
#            upG = row.z-dem_groundPoints_df['z'][indexFinl[0]]
#            upGroundPoints.append([row.x,row.y,upG])#
#    return upGroundPoints
   

#%% dem snow off (veg) for vcm
filenameS0fvcm = 'W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/vcm_snw0ff_dem.tif' #path to raster
#pathNameS0fvcm = 'W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/dem_snow0ff_vcm.png'
elevationVegvcm = readPlotDEM(filenameS0fvcm,elevationMissNoS0f)#,pathNameS0fvcm
dem_groundPointsVegVcm, dem_groundPointsVegVcm_df = creatingCentroidGroundpointsFromDem(filenameS0fvcm,elevationMissNoS0f)#,pathNameS)
lat_vcm, long_vcm, nrows_vcm, ncols_vcm = creatingLatituteLongitudeFromDem(filenameS0fvcm,elevationMissNoS0f)

#calculating slope, aspect, northness, and topo_dimension
elevNrth_gp_vcm, index_vcm = calculatingSlopeAspectNorthness(filenameS0fvcm,'slopeVcm.tif','aspectVcm.tif',elevationVegvcm,nrows_vcm,ncols_vcm,lat_vcm,long_vcm)
topo_dimension_vcm = classificationNorthnessElev(index_vcm,2830,3140)

#cutting some part of dem map
minXVcm = np.min(lat_vcm)+2000
maxXVcm = minXVcm+150
minYVcm = np.min(long_vcm)+2500
maxYVcm = minYVcm+150
dem_groundPointsVegVcm_intg = cutPart0fMap(dem_groundPointsVegVcm_df,minXVcm,maxXVcm,minYVcm,maxYVcm)

#las file snow off (snw) for vcm
lasFilePath_vegVcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/jemezSnow0ffpoints.las"
coordsVegVcm_df = readLasFile (lasFilePath_vegVcm)
#validation las file for snow 0n and snow off (road map)
coordsVegVcm_df['z'] = coordsVegVcm_df['z']-0.081
coordsVegVcm_intg = cutPart0fMap(coordsVegVcm_df,minXVcm,maxXVcm,minYVcm,maxYVcm)

dem_groundPointsVegVcm_intg_ls = np.vstack([dem_groundPointsVegVcm_intg['x'],dem_groundPointsVegVcm_intg['y'],dem_groundPointsVegVcm_intg['z']]).T

upGroundPointsVegVcm0 = np.vstack([coordsVegVcm_intg['x'],coordsVegVcm_intg['y'],coordsVegVcm_intg['z'],coordsVegVcm_intg['z']]).T
for xyvalue in dem_groundPointsVegVcm_intg_ls:
    upGroundPointsVegVcm0[:,3][(upGroundPointsVegVcm0[:,0]==xyvalue[0])&(upGroundPointsVegVcm0[:,1]==xyvalue[1])]=xyvalue[2] #
upGroundPointsVegVcm = np.vstack([upGroundPointsVegVcm0[:,0],upGroundPointsVegVcm0[:,1],upGroundPointsVegVcm0[:,2]-upGroundPointsVegVcm0[:,3]]).T

#las file snow on (snw) for vcm
lasFilePath_snwVcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/jemezSnow0npoints.las"
coordsSnwVcm_df = readLasFile (lasFilePath_snwVcm)
coordsSnwVcm_intg = cutPart0fMap(coordsSnwVcm_df,minXVcm,maxXVcm,minYVcm,maxYVcm)

upGroundPointsSnwVcm0 = np.vstack([coordsSnwVcm_intg['x'],coordsSnwVcm_intg['y'],coordsSnwVcm_intg['z'],coordsSnwVcm_intg['z']]).T
for xyvalue in dem_groundPointsVegVcm_intg_ls:
    upGroundPointsSnwVcm0[:,3][(upGroundPointsSnwVcm0[:,0]==xyvalue[0])&(upGroundPointsSnwVcm0[:,1]==xyvalue[1])]=xyvalue[2] #
upGroundPointsSwnVcm = np.vstack([upGroundPointsSnwVcm0[:,0],upGroundPointsSnwVcm0[:,1],upGroundPointsSnwVcm0[:,2]-upGroundPointsSnwVcm0[:,3]]).T
#%% classification with grountpoints from veg dem file :VCM
#upGroundPointsVegVcm = calculatingUpGroundPoints (coordsVegVcm_intg,dem_groundPointsVegVcm_intg)
# classification of snow las file with veg_grountpoints
#upGroundPointsSnwVegVcm = calculatingUpGroundPoints (coordsSnwVcm_intg,dem_groundPointsVegVcm_intg)
   
#%% vegtation classification from DEM2010 and las2010 snow off
vegClassVcm = upGroundPointsVegVcm.copy()
#vegClassVcm_df = pd.DataFrame(vegClassVcm,columns=['x','y','z'])

#all tree classification based on veg dem
#negVegClassVcm = defineSpecificClassLess (vegClassVcm_df, 0)
allTreeReturnVcm = defineSpecificClassGreater (vegClassVcm, 2)

#finding pixels that have trees
allTreeClassVcm0 = np.vstack([vegClassVcm[:,0],vegClassVcm[:,1],vegClassVcm[:,2],vegClassVcm[:,2]]).T
for xyvalue in allTreeReturnVcm:
    allTreeClassVcm0[:,3][(allTreeClassVcm0[:,0]==xyvalue[0])&(allTreeClassVcm0[:,1]==xyvalue[1])]=1 #
allTreeClassVcm = allTreeClassVcm0[allTreeClassVcm0[:,3]==1]

# trees with low branches
lowVegTreeReturnVcm = defineLowVegClass(allTreeClassVcm)

#finding pixels that have a no low branches
for xyvalue in lowVegTreeReturnVcm:
    allTreeClassVcm[:,3][(allTreeClassVcm[:,0]==xyvalue[0])&(allTreeClassVcm[:,1]==xyvalue[1])]=2 #
allTreeClassNoLowVegVcm = allTreeClassVcm[allTreeClassVcm[:,3]==1]

# trees with no low branches
#concTreeClsLowVegCls = pd.concat([allTreeClassVcm_df,lowVegTreeClassVcm[['x','y','z']]])
#allTreeClassNoLowVeg = concTreeClsLowVegCls.drop_duplicates(keep=False, inplace=False)

#all 0pen classification based on veg dem
all0penClassVcm = allTreeClassVcm0[allTreeClassVcm0[:,3]!=1]

# 0pen with with no grass
all0penClassNoLowVeg = all0penClassVcm[all0penClassVcm[:,3]<=0.15]

#%% snow classification from DEM2010 (snow off) and las2010 snow on 
snowClassVcm = upGroundPointsSwnVcm.copy()

# tree pixcels that do not have a low branch
allTreeSnowNoLowVegClassVcm0 = np.vstack([snowClassVcm[:,0],snowClassVcm[:,1],snowClassVcm[:,2],snowClassVcm[:,2]]).T
for xyvalue in allTreeClassNoLowVegVcm:
    allTreeSnowNoLowVegClassVcm0[:,3][(allTreeSnowNoLowVegClassVcm0[:,0]==xyvalue[0])&(allTreeSnowNoLowVegClassVcm0[:,1]==xyvalue[1])]=1 #
allTreeSnowNoLowVegClassVcm = allTreeSnowNoLowVegClassVcm0[allTreeSnowNoLowVegClassVcm0[:,3]==1]

# snow under canopy
allsnowUnderTreeClassVcm = defineSpecificClassLess (allTreeSnowNoLowVegClassVcm, 2)
#allsnowUnderTreeClassVcm2 = defineSpecificClassGreater (allsnowUnderTreeClassVcm, 0.15))

#snow on the ground with no low branch
for xyvalue in all0penClassNoLowVeg:
    allTreeSnowNoLowVegClassVcm0[:,3][(allTreeSnowNoLowVegClassVcm0[:,0]==xyvalue[0])&(allTreeSnowNoLowVegClassVcm0[:,1]==xyvalue[1])]=2 #
allSnow0nGroundClassVcm = allTreeSnowNoLowVegClassVcm0[allTreeSnowNoLowVegClassVcm0[:,3]==2]

#no snow on the ground
allNosnowClassVcm = defineSpecificClassLess (allSnow0nGroundClassVcm, 0.15)
#%% all classification based on topography and northness
for xyvalue in topo_dimension_vcm['eh']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=1 #
Snow0nGround_eh_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==1]

for xyvalue in topo_dimension_vcm['el']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=2 #
Snow0nGround_el_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==2]

for xyvalue in topo_dimension_vcm['em']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=3 #
Snow0nGround_em_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==3]

for xyvalue in topo_dimension_vcm['nh']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=4 #
Snow0nGround_nh_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==4]

for xyvalue in topo_dimension_vcm['nl']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=5 #
Snow0nGround_nl_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==5]

for xyvalue in topo_dimension_vcm['nm']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=6 #
Snow0nGround_nm_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==6]
    
for xyvalue in topo_dimension_vcm['sh']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=7 #
Snow0nGround_sh_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==7]

for xyvalue in topo_dimension_vcm['sl']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=8 #
Snow0nGround_sl_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==8]

for xyvalue in topo_dimension_vcm['sm']:
    allSnow0nGroundClassVcm[:,3][(allSnow0nGroundClassVcm[:,0]==xyvalue[0])&(allSnow0nGroundClassVcm[:,1]==xyvalue[1])]=9 #
Snow0nGround_sm_vcm = allSnow0nGroundClassVcm[allSnow0nGroundClassVcm[:,3]==9]

#%%
for xyvalue in topo_dimension_vcm['eh']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=1 #
SnowUnderTree_eh_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==1]

for xyvalue in topo_dimension_vcm['el']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=2 #
SnowUnderTree_el_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==2]

for xyvalue in topo_dimension_vcm['em']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=3 #
SnowUnderTree_em_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==3]

for xyvalue in topo_dimension_vcm['nh']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=4 #
SnowUnderTree_nh_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==4]

for xyvalue in topo_dimension_vcm['nl']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=5 #
SnowUnderTree_nl_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==5]

for xyvalue in topo_dimension_vcm['nm']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=6 #
SnowUnderTree_nm_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==6]
    
for xyvalue in topo_dimension_vcm['sh']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=7 #
SnowUnderTree_sh_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==7]

for xyvalue in topo_dimension_vcm['sl']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=8 #
SnowUnderTree_sl_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==8]

for xyvalue in topo_dimension_vcm['sm']:
    allsnowUnderTreeClassVcm[:,3][(allsnowUnderTreeClassVcm[:,0]==xyvalue[0])&(allsnowUnderTreeClassVcm[:,1]==xyvalue[1])]=9 #
SnowUnderTree_sm_vcm = allsnowUnderTreeClassVcm[allsnowUnderTreeClassVcm[:,3]==9]
#%% box plot
d1 = [UnTree_sheltered_midElev_vcm_df['z'],UnTree_neutral_midElev_vcm_df['z'],UnTree_exposed_midElev_vcm_df['z'],
      open_sheltered_midElev_vcm_df['z'],open_neutral_midElev_vcm_df['z'],open_exposed_midElev_vcm_df['z']]

fig, ax = plt.subplots(1,1, figsize=(20,15))
bp1 = ax.boxplot(d1, patch_artist=True)
bp1['boxes'][0].set(color='darkseagreen', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][1].set(color='green', linewidth=2, facecolor = 'olive', hatch = '/')
bp1['boxes'][2].set(color='darkgreen', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][3].set(color='orange', linewidth=2, facecolor = 'olive', hatch = '/')
bp1['boxes'][4].set(color='red', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][5].set(color='darkred', linewidth=2, facecolor = 'olive', hatch = '/')

plt.title('snow depth in Jemez_VCM; under canopy vs. open, exposed vs. sheltered',fontsize=30)

plt.xticks([1,2,3,4,5,6], ['canopy_sheltered','canopy_neutral','canopy_exposed','0pen_sheltered','0pen_neutral','0pen_exposed'], fontsize=20, rotation=30)#'open2016','open2017',
plt.yticks(fontsize=30)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('snow depth (m)', fontsize=30)
plt.savefig('snowDepth_JemezVCM_elev_topo.png')
#%% ploting
#figVcm = plt.figure(figsize=(20,15))
#axVcm = Axes3D(figVcm)
#axVcm.scatter(dem_groundPointsVegVcm_sp[:, 0], dem_groundPointsVegVcm_sp[:, 1], dem_groundPointsVegVcm_sp[:, 2])
#axVcm.scatter(coordsSnwVcm_sp[:, 0], coordsSnwVcm_sp[:, 1], coordsSnwVcm_sp[:, 2])
##axVcm.scatter(coordsVegVcm_sp[:, 0], coordsVegVcm_sp[:, 1], coordsVegVcm_sp[:, 2])
#axVcm.legend()
#plt.title('snow Lidar data JemezVcm',fontsize=30)
#
##for flcl in failclass2:
##    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
#plt.savefig('H:\Chapter2_snow_forest\LIDAR_analysis\JemezVcm\dem_groundPointsVeg&coordsVegSnwVcm.png')

