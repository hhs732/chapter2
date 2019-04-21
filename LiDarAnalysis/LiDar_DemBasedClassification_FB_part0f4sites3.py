import laspy as ls
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import csv
import os

class K_Means:
    def __init__(self, numOfClusters=2, init_centroids=None):
        self.numOfClusters = numOfClusters
        self.centroids={}        
        for i in range(self.numOfClusters):
            self.centroids[i] = init_centroids[i]

    def fit(self,data,cols,cole):
        self.classifications = {}

        for i in range(self.numOfClusters):
            self.classifications[i] = []

        for featureset in data:
            distances = [np.linalg.norm(featureset[cols:cole]-self.centroids[centroid]) for centroid in self.centroids]
            classification = distances.index(min(distances))
            self.classifications[classification].append(featureset)

    def predict(self,data):
        distances = [np.linalg.norm(data-self.centroids[centroid]) for centroid in self.centroids]
        classification = distances.index(min(distances))
        return classification

def readPlotDEM(filename,elevationMissNo,pathName):
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
    
    plt.imshow(elevation, cmap='gist_earth', extent=extent)
    plt.savefig(pathName)
    
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
    
    return dem_groundPoints

# Grab just the X dimension from the file, and scale it.
def scaled_x_dimension(las_file):
    x_dimension = las_file.X
    scale = las_file.header.scale[0]
    offset = las_file.header.offset[0]
    return(x_dimension*scale + offset)
    
def lidarDiffGrndPoints(classes,dem_groundPoints):
    upGroundPoints = []
    classes_rplc = [] #returns filled class with [0,0,0]
    pureDiff = []
    for clss in range (len (classes)):
        if len(classes[clss])==0:
            nokhale = [np.array([0,0,0])]
            classes_rplc.append(nokhale)
        else: classes_rplc.append(classes[clss])
        
        pureDiff.append(classes_rplc[clss]-dem_groundPoints[clss])  
     
        eachpoint = []
        for ep in range(len(classes_rplc[clss])):
            height = classes_rplc[clss][ep][2]-dem_groundPoints[clss][2]
            eachpoint.append(np.vstack([classes_rplc[clss][ep][0],classes_rplc[clss][ep][1],height]).T)

        upGroundPoints.append(eachpoint)
    return upGroundPoints, classes_rplc, pureDiff

def classificationTest(pureDiff):
    failclass=[]
    for xyidx in range (len(pureDiff)):
        for xycl in range (len(pureDiff[xyidx])):
            if ((abs(pureDiff[xyidx][xycl][0])>0.5) or (abs(pureDiff[xyidx][xycl][1])>0.5)):
                failclass.append(xyidx)
                break
    return failclass

def defineSpecificClassGreater(classesG, specific0bjectHeightG):
    specificClassG = []
    numSpecificClassG = []
    for vgcl in range (len (classesG)):
        for pnt in range (len (classesG[vgcl])):
            if classesG[vgcl][pnt][0][2]>specific0bjectHeightG:
                numSpecificClassG.append(vgcl)
                specificClassG.append(classesG[vgcl])
                break
    return specificClassG, numSpecificClassG

def defineSpecificClassLess (classesL, specific0bjectHeightL):
    specificClassL = []
    numSpecificClassL = []
    for vgcl in range (len (classesL)):
        for pnt in range (len (classesL[vgcl])):
            if classesL[vgcl][pnt][0][2]<specific0bjectHeightL:
                numSpecificClassL.append(vgcl)
                specificClassL.append(classesL[vgcl])
                break
    return specificClassL, numSpecificClassL

def defineLowVegClass(classes):
    lowVegClass = []
    lowVegNumClass = []
    for lvgcl in range (len (classes)):
        for crdnt in range (len (classes[lvgcl])):
            if classes[lvgcl][crdnt][0][2]<2 and classes[lvgcl][crdnt][0][2]>0.15:
                lowVegNumClass.append(lvgcl)
                lowVegClass.append(classes[lvgcl])
                break
    return lowVegClass,lowVegNumClass

def differenceBetwee2classes (primaryClass,secondNumClass): #to define nolowVegClass and openClass
    primaryNumClass = list(np.arange(len(primaryClass)))
    cleanNumClass = list(set(primaryNumClass)-set(secondNumClass))
    cleanClass = []
    for indx in cleanNumClass:
        cleanClass.append(primaryClass[indx])
    return cleanClass, cleanNumClass

#%% DEM file (.tif) reading and creating centroid (ground points) from Dem files for sagehen
elevationMissNoS0f = 1900.

#lidar data for whole sagehen, snow on 2016
from LiDar_DemBasedClassification_scSnow0n import coordsnow0nSc

filenameS0fT1 = "W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT1\T1_snow0ff_dem.tif" #path to raster
pathNameS0fT1 = "W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT1\dem_snow0ff_sagehenT1.png"
elevationVegT1 = readPlotDEM(filenameS0fT1,elevationMissNoS0f,pathNameS0fT1)
dem_groundPointsVegT1 = creatingCentroidGroundpointsFromDem(filenameS0fT1,elevationMissNoS0f)#,pathNameS)

# LiDar Data reading and Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array
infileVegT1 = ls.file.File("W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT1\T1_snow0ff_points.las", mode="r")
coordsVegT1 = np.vstack((infileVegT1.x, infileVegT1.y, infileVegT1.z)).T
coordsVegT1_df = pd.DataFrame(np.vstack((infileVegT1.x, infileVegT1.y, infileVegT1.z)).T,columns=['x','y','z'])
#finding boarder of las file veg T1
maxX = coordsVegT1_df.iloc[coordsVegT1_df['x'].idxmax()][0:1]
minX = coordsVegT1_df.iloc[coordsVegT1_df['x'].idxmin()][0:1]
maxY = coordsVegT1_df.iloc[coordsVegT1_df['y'].idxmax()][1:2]
minY = coordsVegT1_df.iloc[coordsVegT1_df['y'].idxmin()][1:2]
#borderT1_init = pd.DataFrame([maxX,maxY,minX,minY])
#borderT1 = pd.DataFrame([[borderT1_init['x'].max(),borderT1_init['y'].min()],[borderT1_init['x'].max(),borderT1_init['y'].max()],
#                        [borderT1_init['x'].min(),borderT1_init['y'].min()],[borderT1_init['x'].min(),borderT1_init['y'].max()]])
#lidar data for sagehenT1, snow on 2016
coordsnow0nScT1_int = coordsnow0nSc[(coordsnow0nSc['x'] >= minX[0]) & (coordsnow0nSc['x'] <= maxX[0])]
coordSnwT10 = coordsnow0nScT1_int[(coordsnow0nScT1_int['y'] >= minY[0]) & (coordsnow0nScT1_int['y'] <= maxY[0])]
coordsSnwT1 = coordSnwT10.values
#%% classification with veggrountpoints from veg dem file :T1
centroids_newT1=dem_groundPointsVegT1[:,0:2]  
kT1 = np.size(dem_groundPointsVegT1[:,0])
# instantiate a class
clf1T1 = K_Means(numOfClusters=kT1,init_centroids=centroids_newT1)
# fit kmean class to data
clf1T1.fit(coordsVegT1,0,2)
# get classification 
classesVegT1 = clf1T1.classifications

upGroundPointsVegT1, classes_rplcVegT1, pureVegClassT1 = lidarDiffGrndPoints(classesVegT1,dem_groundPointsVegT1)

#%% classification of snow las file with veg_grountpoints
clfsT1 = K_Means(numOfClusters=kT1,init_centroids=centroids_newT1)
# fit kmean class to data
clfsT1.fit(coordsSnwT1,0,2)
# get classification 
classeSnwT1 = clfsT1.classifications

upGroundPointsSnwVegT1, classes_rplcVST1, pureVegsnowClassT1 = lidarDiffGrndPoints(classeSnwT1,dem_groundPointsVegT1)
#%% vegtation classification from DEM2014 and las2014
vegClassT1 = upGroundPointsVegT1[:]
#all tree classification
allTreeClassT1, treeNumClassT1 = defineSpecificClassGreater (vegClassT1, 2)
        
negVegClassT1, negVegNumClassT1 = defineSpecificClassLess (vegClassT1, 0)

# trees with low branches
lowVegTreeClassT1, lowVegNumClassT1 = defineLowVegClass(allTreeClassT1)

# trees with no low blanches
# "*************tall canopy no snow class*****************"
nolowVegTreeClassT1, nolowVegTreeNumClassT1 = differenceBetwee2classes (allTreeClassT1,lowVegNumClassT1)

# open space (no trees, no return between 0.15 to 2)
# all low veg
allLowVegClassT1, allLowVegNumClassT1 = defineLowVegClass(vegClassT1)

#open places
notOpenNumClassT1 = list(set(allLowVegNumClassT1).union(set(treeNumClassT1)))
# "******************open no snow class*******************"
allOpenClassT1, allOpenNumClassT1 = differenceBetwee2classes (vegClassT1,notOpenNumClassT1)

#%% snow classification
vegSnowClassT1, vegSnowNumClassT1 = defineSpecificClassGreater (upGroundPointsSnwVegT1, -1)

#no snow on the ground
nosnowClassT1, noSnowNumClassT1 = defineSpecificClassLess (vegSnowClassT1, 0.15)
#snow on the ground or on the trees
allSnowClassT1, allSnowNumClassT1 = differenceBetwee2classes (vegSnowClassT1,noSnowNumClassT1)
# "******************tall canopy snow class*******************"
#snow on the tall canopy >2m
treeSnowClassT1, treeSnowNumClassT1 = defineSpecificClassGreater (allSnowClassT1, 2)
# "******************open snow class**************************"
#snow on the ground
groundSnowClassT1, groundSnowNumClassT1 = differenceBetwee2classes (allSnowClassT1,treeSnowNumClassT1)

#%% ploting
fig3 = plt.figure(figsize=(20,15))
ax3 = Axes3D(fig3)
ax3.scatter(dem_groundPointsVegT1[:, 0], dem_groundPointsVegT1[:, 1], dem_groundPointsVegT1[:, 2])
#ax3.scatter(coordsVegT1[:, 0], coordsVegT1[:, 1], coordsVegT1[:, 2])
ax3.scatter(coordsSnwT1[:, 0], coordsSnwT1[:, 1], coordsSnwT1[:, 2])
ax3.legend()
#for flcl in failclass2:
#    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
plt.savefig('H:\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT1\dem_groundPointsVegT1&coordsSnwT1.png')

#%%
#%%
#%%
#%%
#%% DEM file (.tif) reading and creating centroid (ground points) from Dem files for sagehen
filenameS0fT4 = "W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT4\T4_snow0ff_dem.tif" #path to raster
pathNameS0fT4 = "W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT4\dem_snow0ff_sagehenT4.png"
elevationVegT4 = readPlotDEM(filenameS0fT4,elevationMissNoS0f,pathNameS0fT4)
dem_groundPointsVegT4 = creatingCentroidGroundpointsFromDem(filenameS0fT4,elevationMissNoS0f)#,pathNameS)

# LiDar Data reading and Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array
infileVegT4 = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/sagehenT4/T4_snow0ff_points.las", mode="r")
coordsVegT4 = np.vstack((infileVegT4.x, infileVegT4.y, infileVegT4.z)).T
coordsVegT4_df = pd.DataFrame(np.vstack((infileVegT4.x, infileVegT4.y, infileVegT4.z)).T,columns=['x','y','z'])
#finding boarder of las file for T4
maxX4 = coordsVegT4_df.iloc[coordsVegT4_df['x'].idxmax()][0:1]
minX4 = coordsVegT4_df.iloc[coordsVegT4_df['x'].idxmin()][0:1]
maxY4 = coordsVegT4_df.iloc[coordsVegT4_df['y'].idxmax()][1:2]
minY4 = coordsVegT4_df.iloc[coordsVegT4_df['y'].idxmin()][1:2]
borderT4 = pd.DataFrame([[maxX4[0],minY4[0]],[maxX4[0],maxY4[0]],[minX4[0],minY4[0]],[minX4[0],maxY4[0]]])
#lidar data for sagehenT4, snow on 2016
coordsnow0nScT4_int = coordsnow0nSc[(coordsnow0nSc['x'] >= minX4[0]) & (coordsnow0nSc['x'] <= maxX4[0])]
coordSnwT40 = coordsnow0nScT4_int[(coordsnow0nScT4_int['y'] >= minY4[0]) & (coordsnow0nScT4_int['y'] <= maxY4[0])]
coordsSnwT4 = coordSnwT40.values

#%% classification with grountpoints from veg dem file :T4
centroids_newT4=dem_groundPointsVegT4[:,0:2]  
kT4 = np.size(dem_groundPointsVegT4[:,0])
# instantiate a class
clf1T4 = K_Means(numOfClusters=kT4,init_centroids=centroids_newT4)
# fit kmean class to data
clf1T4.fit(coordsVegT4,0,2)
# get classification 
classesVegT4 = clf1T4.classifications

upGroundPointsVegT4, classes_rplcVegT4, pureVegClassT4 = lidarDiffGrndPoints(classesVegT4,dem_groundPointsVegT4)

#%% classification of snow las file with veg_grountpoints
clfsT4 = K_Means(numOfClusters=kT4,init_centroids=centroids_newT4)
# fit kmean class to data
clfsT4.fit(coordsSnwT4,0,2)
# get classification 
classeSnwT4 = clfsT4.classifications

upGroundPointsSnwVegT4, classes_rplcVST4, pureVegsnowClassT4 = lidarDiffGrndPoints(classeSnwT4,dem_groundPointsVegT4)
#%% vegtation classification from DEM2014 and las2014
vegClassT4 = upGroundPointsVegT4[:]
#all tree classification
allTreeClassT4, treeNumClassT4 = defineSpecificClassGreater (vegClassT4, 2)
        
negVegClassT4, negVegNumClassT4 = defineSpecificClassLess (vegClassT4, 0)

# trees with low branches
lowVegTreeClassT4, lowVegNumClassT4 = defineLowVegClass(allTreeClassT4)

# trees with no low blanches
# "*************tall canopy no snow class*****************"
nolowVegTreeClassT4, nolowVegTreeNumClassT4 = differenceBetwee2classes (allTreeClassT4,lowVegNumClassT4)

# open space (no trees, no return between 0.15 to 2)
# all low veg
allLowVegClassT4, allLowVegNumClassT4 = defineLowVegClass(vegClassT4)

#open places
notOpenNumClassT4 = list(set(allLowVegNumClassT4).union(set(treeNumClassT4)))
# "******************open no snow class*******************"
allOpenClassT4, allOpenNumClassT4 = differenceBetwee2classes (vegClassT4,notOpenNumClassT4)

#%% snow classification
vegSnowClassT4, vegSnowNumClassT4 = defineSpecificClassGreater (upGroundPointsSnwVegT4, -1)

#no snow on the ground
nosnowClassT4, noSnowNumClassT4 = defineSpecificClassLess (vegSnowClassT4, 0.15)
#snow on the ground or on the trees
allSnowClassT4, allSnowNumClassT4 = differenceBetwee2classes (vegSnowClassT4,noSnowNumClassT4)
# "******************tall canopy snow class*******************"
#snow on the tall canopy >2m
treeSnowClassT4, treeSnowNumClassT4 = defineSpecificClassGreater (allSnowClassT4, 2)
# "******************open snow class**************************"
#snow on the ground
groundSnowClassT4, groundSnowNumClassT4 = differenceBetwee2classes (allSnowClassT4,treeSnowNumClassT4)
#%% ploting
fig4 = plt.figure(figsize=(20,15))
ax4 = Axes3D(fig4)
ax4.scatter(dem_groundPointsVegT4[:, 0], dem_groundPointsVegT4[:, 1], dem_groundPointsVegT4[:, 2])
#ax4.scatter(coordsSnwT4[:, 0], coordsSnwT4[:, 1], coordsSnwT4[:, 2])
ax4.scatter(coordsVegT4[:, 0], coordsVegT4[:, 1], coordsVegT4[:, 2])
ax4.legend()
#for flcl in failclass2:
#    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
plt.savefig('H:\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\sagehenT4\dem_groundPointsVegT1&coordsVegT4.png')
#%%
#%%
#%%
#%%
#%% dem snow off (veg) for vcm
filenameS0fvcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/jemez_vcm/vcm_snow0ff_dem.tif" #path to raster
pathNameS0fvcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/jemez_vcm/dem_snow0ff_vcm.png"
elevationVegvcm = readPlotDEM(filenameS0fvcm,elevationMissNoS0f,pathNameS0fvcm)
dem_groundPointsVegVcm = creatingCentroidGroundpointsFromDem(filenameS0fvcm,elevationMissNoS0f)#,pathNameS)
#dem snow on (snw) for vcm
filenameS0nVcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/jemez_vcm/vcm_snow0n_dem.tif" #path to raster
pathNameS0nVcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/jemez_vcm/dem_snow0n_vcm.png"
elevationSnowVcm = readPlotDEM(filenameS0nVcm,elevationMissNoS0f,pathNameS0nVcm)
dem_groundPointsSnowVcm = creatingCentroidGroundpointsFromDem(filenameS0nVcm,elevationMissNoS0f)#,pathNameS)
#las file snow off (snw) for vcm
infileVegVcm = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/jemez_vcm/vcm_snow0ff_points.las", mode="r")
coordsVegVcm = np.vstack((infileVegVcm.x, infileVegVcm.y, infileVegVcm.z)).T
#las file snow on (snw) for vcm
infileSnwVcm = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/jemez_vcm/vcm_snow0n_points.las", mode="r")
coordsSnwVcm = np.vstack((infileSnwVcm.x, infileSnwVcm.y, infileSnwVcm.z)).T

#%% classification with grountpoints from veg dem file :VCM
centroids_newVcm=dem_groundPointsVegVcm[:,0:2]  
kVcm = np.size(dem_groundPointsVegVcm[:,0])
# instantiate a class
clf1Vcm = K_Means(numOfClusters=kVcm,init_centroids=centroids_newVcm)
# fit kmean class to data
clf1Vcm.fit(coordsVegVcm,0,2)
# get classification 
classesVegVcm = clf1Vcm.classifications

upGroundPointsVegVcm, classes_rplcVegVcm, pureVegClassVcm = lidarDiffGrndPoints(classesVegVcm,dem_groundPointsVegVcm)

#%% classification of snow las file with veg_grountpoints
clfsVcm = K_Means(numOfClusters=kVcm,init_centroids=centroids_newVcm)
# fit kmean class to data
clfsVcm.fit(coordsSnwVcm,0,2)
# get classification 
classeSnwVcm = clfsVcm.classifications

upGroundPointsSnwVegVcm, classes_rplcVSVcm, pureVegsnowClassVcm = lidarDiffGrndPoints(classeSnwVcm,dem_groundPointsVegVcm)
#%% vegtation classification from DEM2014 and las2014
vegClassVcm = upGroundPointsVegVcm[:]
#all tree classification
allTreeClassVcm, treeNumClassVcm = defineSpecificClassGreater (vegClassVcm, 2)
        
negVegClassVcm, negVegNumClassVcm = defineSpecificClassLess (vegClassVcm, 0)

# trees with low branches
lowVegTreeClassVcm, lowVegNumClassVcm = defineLowVegClass(allTreeClassVcm)

# trees with no low blanches
# "*************tall canopy no snow class*****************"
nolowVegTreeClassVcm, nolowVegTreeNumClassVcm = differenceBetwee2classes (allTreeClassVcm,lowVegNumClassVcm)

# open space (no trees, no return between 0.15 to 2)
# all low veg
allLowVegClassVcm, allLowVegNumClassVcm = defineLowVegClass(vegClassVcm)

#open places
notOpenNumClassVcm = list(set(allLowVegNumClassVcm).union(set(treeNumClassVcm)))
# "******************open no snow class*******************"
allOpenClassVcm, allOpenNumClassVcm = differenceBetwee2classes (vegClassVcm,notOpenNumClassVcm)

#%% snow classification
vegSnowClassVcm, vegSnowNumClassVcm = defineSpecificClassGreater (upGroundPointsSnwVegVcm, -1)

#no snow on the ground
nosnowClassVcm, noSnowNumClassVcm = defineSpecificClassLess (vegSnowClassVcm, 0.15)
#snow on the ground or on the trees
allSnowClassVcm, allSnowNumClassVcm = differenceBetwee2classes (vegSnowClassVcm,noSnowNumClassVcm)
# "******************tall canopy snow class*******************"
#snow on the tall canopy >2m
treeSnowClassVcm, treeSnowNumClassVcm = defineSpecificClassGreater (allSnowClassVcm, 2)
# "******************open snow class**************************"
#snow on the ground
groundSnowClassVcm, groundSnowNumClassVcm = differenceBetwee2classes (allSnowClassVcm,treeSnowNumClassVcm)
#%% ploting
figVcm = plt.figure(figsize=(20,15))
axVcm = Axes3D(figVcm)
axVcm.scatter(dem_groundPointsVegVcm[:, 0], dem_groundPointsVegVcm[:, 1], dem_groundPointsVegVcm[:, 2])
#axVcm.scatter(coordsSnwVcm[:, 0], coordsSnwVcm[:, 1], coordsSnwVcm[:, 2])
axVcm.scatter(coordsVegVcm[:, 0], coordsVegVcm[:, 1], coordsVegVcm[:, 2])
axVcm.legend()
#for flcl in failclass2:
#    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
plt.savefig('H:\Chapter2_snow_forest\LIDAR_analysis\chunk0fdata\jemez_vcm\dem_groundPointsVegT1&coordsVegVcm.png')
#%%
#%%
#%%
#%%
#%%
#%% dem snow off (veg) for NR1
filenameS0fNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/nr1_snow0ff_dem.tif" #path to raster
pathNameS0fNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/dem_snow0ff_Nr.png"
elevationVegNr1 = readPlotDEM(filenameS0fNr1,elevationMissNoS0f,pathNameS0fNr1)
dem_groundPointsVegNr1 = creatingCentroidGroundpointsFromDem(filenameS0fNr1,elevationMissNoS0f)#,pathNameS)
# dem snow on (snw) for NR1
filenameS0nNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/nr1_snow0n_dem.tif" #path to raster
pathNameS0nNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/dem_snow0n_Nr1.png"
elevationSnowNr1 = readPlotDEM(filenameS0nNr1,elevationMissNoS0f,pathNameS0nNr1)
dem_groundPointsSnowNr1 = creatingCentroidGroundpointsFromDem(filenameS0nNr1,elevationMissNoS0f)#,pathNameS)
#las file snow off (snw) for NR1
infileVegNr1 = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/nr1_snow0ff_points.las", mode="r")
coordsVegNr1 = np.vstack((infileVegNr1.x, infileVegNr1.y, infileVegNr1.z)).T
#las file snow on (snw) for NR1
infileSnwNr1 = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/nr1_snow0n_points.las", mode="r")
coordsSnwNr1 = np.vstack((infileSnwNr1.x, infileSnwNr1.y, infileSnwNr1.z)).T
#%% classification with grountpoints from veg dem file :NR1
centroids_newNr1=dem_groundPointsVegNr1[:,0:2]  
kNr1 = np.size(dem_groundPointsVegNr1[:,0])
# instantiate a class
clf1Nr1 = K_Means(numOfClusters=kNr1,init_centroids=centroids_newNr1)
# fit kmean class to data
clf1Nr1.fit(coordsVegNr1,0,2)
# get classification 
classesVegNr1 = clf1Nr1.classifications

upGroundPointsVegNr1, classes_rplcVegNr1, pureVegClassNr1 = lidarDiffGrndPoints(classesVegNr1,dem_groundPointsVegNr1)

#%% classification of snow las file with veg_grountpoints
clfsNr1 = K_Means(numOfClusters=kNr1,init_centroids=centroids_newNr1)
# fit kmean class to data
clfsNr1.fit(coordsSnwNr1,0,2)
# get classification 
classeSnwNr1 = clfsNr1.classifications

upGroundPointsSnwVegNr1, classes_rplcVSNr1, pureVegsnowClassNr1 = lidarDiffGrndPoints(classeSnwNr1,dem_groundPointsVegNr1)

#%% vegtation classification from DEM2014 and las2014
vegClassNr1 = upGroundPointsVegNr1[:]
#all tree classification
allTreeClassNr1, treeNumClassNr1 = defineSpecificClassGreater (vegClassNr1, 2)
        
negVegClassNr1, negVegNumClassNr1 = defineSpecificClassLess (vegClassNr1, 0)

# trees with low branches
lowVegTreeClassNr1, lowVegNumClassNr1 = defineLowVegClass(allTreeClassNr1)

# trees with no low blanches
# "*************tall canopy no snow class*****************"
nolowVegTreeClassNr1, nolowVegTreeNumClassNr1 = differenceBetwee2classes (allTreeClassNr1,lowVegNumClassNr1)

# open space (no trees, no return between 0.15 to 2)
# all low veg
allLowVegClassNr1, allLowVegNumClassNr1 = defineLowVegClass(vegClassNr1)

#open places
notOpenNumClassNr1 = list(set(allLowVegNumClassNr1).union(set(treeNumClassNr1)))
# "******************open no snow class*******************"
allOpenClassNr1, allOpenNumClassNr1 = differenceBetwee2classes (vegClassNr1,notOpenNumClassNr1)

#%% snow classification
vegSnowClassNr1, vegSnowNumClassNr1 = defineSpecificClassGreater (upGroundPointsSnwVegNr1, -1)

#no snow on the ground
nosnowClassNr1, noSnowNumClassNr1 = defineSpecificClassLess (vegSnowClassNr1, 0.15)
#snow on the ground or on the trees
allSnowClassNr1, allSnowNumClassNr1 = differenceBetwee2classes (vegSnowClassNr1,noSnowNumClassNr1)
# "******************tall canopy snow class*******************"
#snow on the tall canopy >2m
treeSnowClassNr1, treeSnowNumClassNr1 = defineSpecificClassGreater (allSnowClassNr1, 2)
# "******************open snow class**************************"
#snow on the ground
groundSnowClassNr1, groundSnowNumClassNr1 = differenceBetwee2classes (allSnowClassNr1,treeSnowNumClassNr1)
#%%
figNr1 = plt.figure(figsize=(20,15))
axNr1 = Axes3D(figNr1)
axNr1.scatter(dem_groundPointsVegNr1[:, 0], dem_groundPointsVegNr1[:, 1], dem_groundPointsVegNr1[:, 2])
#axNr1.scatter(coordsSnwNr1[:, 0], coordsSnwNr1[:, 1], coordsSnwNr1[:, 2])
axNr1.scatter(coordsVegNr1[:, 0], coordsVegNr1[:, 1], coordsVegNr1[:, 2])
axNr1.legend()
#for flcl in failclass2:
#    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
plt.savefig("H:/Chapter2_snow_forest/LIDAR_analysis/chunk0fdata/niwot_nr1/dem_groundPointsVegT1&coordsVegNr1.png")
#%% ploting
#fig3 = plt.figure(figsize=(20,15))
#ax3 = Axes3D(fig3)
##ax3.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
#ax3.scatter(coordsnow[:, 0], coordsnow[:, 1], coordsnow[:, 2])
#ax3.scatter(coordsVeg[:, 0], coordsVeg[:, 1], coordsVeg[:, 2])
#ax3.legend()
##for flcl in failclass2:
##    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
#plt.savefig('lidardata\dem_lidar_snow&veg.png')





























