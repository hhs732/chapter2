import laspy as ls
import numpy as np
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

#def extract2equalPartof2dems_clipDem(inputfile,muskfile):
#    for indx1 in range (len(inputfile)):#(len(dem_groundPointsVeg)):
#        for indx2 in range (len(muskfile)):#(len(sagehen_demGroundPoint_snow)):
#            if (inputfile[indx1,0]==muskfile[indx2,0]) and (inputfile[indx1,1]==muskfile[indx2,1]):
#                sagehen_demGroundPoint_Veg == dem_groundPointsSnow[np.where((inputfile[indx1,0]==muskfile[indx2,0]) and (inputfile[indx1,1]==muskfile[indx2,1]))]


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
    for cls in range (len (classes)):
        if len(classes[cls])==0:
            nokhale = [np.array([0,0,0])]
            classes_rplc.append(nokhale)
        else: classes_rplc.append(classes[cls])
        
        pureDiff.append(classes_rplc[cls]-dem_groundPoints[cls])  
     
        eachpoint = []
        for ep in range(len(classes_rplc[cls])):
            height = classes_rplc[cls][ep][2]-dem_groundPoints[cls][2]
            eachpoint.append(np.vstack([classes_rplc[cls][ep][0],classes_rplc[cls][ep][1],height]).T)

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

def defineSpecificClassGreater (classesG, specific0bjectHeightG):
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
filenameS0 = "lidardata\SagehenLasDEM\demSnow0n2008.tif" #path to raster
elevationMissNoS0 = 1900.
pathNameS0 = 'lidardata\SagehenLasDEM\dem_snowOn_sagehen.png'
elevationS0 = readPlotDEM(filenameS0,elevationMissNoS0,pathNameS0)
dem_groundPointsSnow = creatingCentroidGroundpointsFromDem(filenameS0,elevationMissNoS0)#,pathNameS)
sagehen_demGroundPoint_snow = dem_groundPointsSnow[np.where(dem_groundPointsSnow[:,2]>1900.)] #SagehenGP = [z for z in dem_groundPointsSnow if z[2]>1900]

filenameS0f = "lidardata\SagehenLasDEM\demSnow0ff2014.tif" #path to raster
elevationMissNoS0f = 1900.
pathNameS0f = 'lidardata\SagehenLasDEM\dem_snow0ff_sagehen.png'
elevationVeg = readPlotDEM(filenameS0f,elevationMissNoS0f,pathNameS0f)
dem_groundPointsVeg = creatingCentroidGroundpointsFromDem(filenameS0f,elevationMissNoS0f)#,pathNameS)
#sagehen_demGroundPoint_veg = dem_groundPointsVeg[np.where((dem_groundPointsVeg[:,0]==sagehen_demGroundPoint_snow[:,0]))] #and (dem_groundPointsVeg[:,1]==sagehen_demGroundPoint_snow[:,1]))] #SagehenGP = [z for z in dem_groundPointsSnow if z[2]>1900]
print 'now'
#%%
#sagehen_demGroundPoint_VegIndex = []
#for z in sagehen_demGroundPoint_snow:
#      sagehen_demGroundPoint_VegIndex.append(np.where(dem_groundPointsVeg[:,0:2]==z[0:2]))
#
#print 'now*'
#%%# LiDar Data reading and Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array; snow on 2016
infileSnow1 = ls.file.File("lidardata\SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_222743_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow1 = np.vstack((infileSnow1.x, infileSnow1.y, infileSnow1.z)).T
infileSnow2 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_223059_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow2 = np.vstack((infileSnow2.x, infileSnow2.y, infileSnow2.z)).T
infileSnow3 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_223357_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow3 = np.vstack((infileSnow3.x, infileSnow3.y, infileSnow3.z)).T
infileSnow4 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_223635_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow4 = np.vstack((infileSnow4.x, infileSnow4.y, infileSnow4.z)).T
infileSnow5 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_223932_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow5 = np.vstack((infileSnow5.x, infileSnow5.y, infileSnow5.z)).T
infileSnow6 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_221955_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow6 = np.vstack((infileSnow6.x, infileSnow6.y, infileSnow6.z)).T
infileSnow7 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_222224_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow7 = np.vstack((infileSnow7.x, infileSnow7.y, infileSnow7.z)).T
infileSnow8 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_222504_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow8 = np.vstack((infileSnow8.x, infileSnow8.y, infileSnow8.z)).T
infileSnow9 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_222743_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow9 = np.vstack((infileSnow9.x, infileSnow9.y, infileSnow9.z)).T
infileSnow10 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_223059_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow10 = np.vstack((infileSnow10.x, infileSnow10.y, infileSnow10.z)).T
infileSnow11 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_223358_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow11 = np.vstack((infileSnow11.x, infileSnow11.y, infileSnow11.z)).T
infileSnow12 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_223635_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow12 = np.vstack((infileSnow12.x, infileSnow12.y, infileSnow12.z)).T #***
infileSnow13 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 2 - 160326_223932_2 - originalpoints_dem_filter.las", mode="r"); coordsSnow13 = np.vstack((infileSnow13.x, infileSnow13.y, infileSnow13.z)).T
infileSnow14 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_221954_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow14 = np.vstack((infileSnow14.x, infileSnow14.y, infileSnow14.z)).T
infileSnow15 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_222224_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow15 = np.vstack((infileSnow15.x, infileSnow15.y, infileSnow15.z)).T
infileSnow16 = ls.file.File("lidardata/SagehenLasDEM/las2016snow0n/USCASH20160326f2a1 - Channel 1 - 160326_222504_1 - originalpoints_dem_filter.las", mode="r"); coordsSnow16 = np.vstack((infileSnow16.x, infileSnow16.y, infileSnow16.z)).T

coordsSnow = np.vstack((coordsSnow1,coordsSnow2,coordsSnow3,coordsSnow4,coordsSnow5,coordsSnow6,coordsSnow7,coordsSnow8,coordsSnow9,
                        coordsSnow10,coordsSnow11,coordsSnow12,coordsSnow13,coordsSnow14,coordsSnow15,coordsSnow16))
#%%# LiDar Data reading and Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array
infileVeg1 = ls.file.File("lidardata\SagehenLasDEM\lasVeg20141.las", mode="r")
infileVeg2 = ls.file.File("lidardata\SagehenLasDEM\lasVeg20142.las", mode="r")
infileVeg3 = ls.file.File("lidardata\SagehenLasDEM\lasVeg20143.las", mode="r")

coordsVeg1 = np.vstack((infileVeg1.x, infileVeg1.y, infileVeg1.z)).T
coordsVeg2 = np.vstack((infileVeg2.x, infileVeg2.y, infileVeg2.z)).T
coordsVeg3 = np.vstack((infileVeg3.x, infileVeg3.y, infileVeg3.z)).T
coordsVeg = np.vstack((coordsVeg1,coordsVeg2,coordsVeg3))

# calculating the nearest neighbors of a set of points, you might want to use a highly optimized package like FLANN 
#datasetVeg1 = np.vstack([infileVeg1.X, infileVeg1.Y, infileVeg1.Z]).T
#minLat,maxLat = np.min(datasetVeg[:,0]),np.max(datasetVeg[:,0])
#minLon, maxLon = np.min(datasetVeg[:,1]),np.max(datasetVeg[:,1])

print 'now2'
#%% classification with new grountpoints from dem file
centroids_new=dem_groundPointsVeg[:,0:2]  
k1 = np.size(dem_groundPointsVeg[:,0])
# instantiate a class
clf1 = K_Means(numOfClusters=k1,init_centroids=centroids_new)
# fit kmean class to data
clf1.fit(coordsVeg,0,2)
# get classification 
classesVeg = clf1.classifications

print 'now3'
#%% test classification
upGroundPointsVeg, classes_rplcVeg, pureVegClass = lidarDiffGrndPoints(classesVeg,coordsVeg)
failclassVeg = classificationTest(pureVegClass)

print 'now3'
#%% ploting
#fig = plt.figure(figsize=(20,15))
#ax = Axes3D(fig)
##ax.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
##ax.scatter(coordsVeg[:, 0], coordsVeg[:, 1], coordsVeg[:, 2])
#for flcl in failclassVeg:
#    #ax.scatter([x[0] for x in classes1[flcl]], [x[1] for x in classes1[flcl]])
#    ax.scatter([x[0] for x in classesVeg[flcl]], [x[1] for x in classesVeg[flcl]], [x[2] for x in classesVeg[flcl]])
##ax.scatter([x[0] for x in classes1[25]], [x[1] for x in classes1[25]])
#plt.savefig('lidardata\SagehenLasDEM\dem_lidar_classificationtest.png')
#%% vegtation classification from DEM2014 and las2014
vegClass = upGroundPointsVeg[:]
#all tree classification
allTreeClass, treeNumClass = defineSpecificClassGreater (vegClass, 2)
        
negVegClass, negVegNumClass = defineSpecificClassLess (vegClass, 0)

# trees with low branches
lowVegTreeClass, lowVegNumClass = defineLowVegClass(allTreeClass)

# trees with no low blanches
# "*************tall canopy no snow class*****************"
nolowVegTreeClass, nolowVegTreeNumClass = differenceBetwee2classes (allTreeClass,lowVegNumClass)

# open space (no trees, no return between 0.15 to 2)
# all low veg
allLowVegClass, allLowVegNumClass = defineLowVegClass(vegClass)

#open places
notOpenNumClass = list(set(allLowVegNumClass).union(set(treeNumClass)))
# "******************open no snow class*******************"
allOpenClass, allOpenNumClass = differenceBetwee2classes (vegClass,notOpenNumClass)
#test=[]
#for i in range (len(nolowVegTreeClass)):
#    for j in range(len(nolowVegTreeClass[i])):
#        test.append([i,nolowVegTreeClass[i][j][0][2]])

print 'now4'
#%%#%% snow on lidar files
infileSnow = ls.file.File("lidardata\SagehenLasDEM\las2016snow0n\USCASH20160326f2a1 - Channel 2 - 160326_223635_2 - originalpoints_dem_filter.las", mode="r")
coordsnow = np.vstack((infileSnow.x, infileSnow.y, infileSnow.z)).T
datasetSnow = np.vstack([infileSnow.X, infileSnow.Y, infileSnow.Z]).T

print 'now5'
#%% classification with new grountpoints from dem file
clfs = K_Means(numOfClusters=k1,init_centroids=centroids_new)
# fit kmean class to data
clfs.fit(coordsnow,0,2)
# get classification 
classesnow = clfs.classifications

print 'now6'
#%% test classification
vegSnowClass, classes_rplcVS, pureVegsnowClass = lidarDiffGrndPoints(classesnow,dem_groundPointsVeg)
failclassSnow = classificationTest(pureVegsnowClass)
vegSnowClass2, vegSnowNumClass = defineSpecificClassGreater (vegSnowClass, -1)

print 'now7'
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
#%% raw classification
rawClass = infileSnow.raw_classification
return_arr = infileSnow.num_returns

print 'now8'
#%% snow classification

#no snow on the ground
nosnowClass, noSnowNumClass = defineSpecificClassLess (vegSnowClass2, 0.15)
#snow on the ground or on the trees
allSnowClass, allSnowNumClass = differenceBetwee2classes (vegSnowClass2,noSnowNumClass)
# "******************tall canopy snow class*******************"
#snow on the tall canopy >2m
treeSnowClass, treeSnowNumClass = defineSpecificClassGreater (allSnowClass, 2)
# "******************open snow class**************************"
#snow on the ground
groundSnowClass, groundSnowNumClass = differenceBetwee2classes (allSnowClass,treeSnowNumClass)

#for gsclss in range(len(groundSnow0lowBranchClass)):
#    for gsdim in range (len(groundSnow0lowBranchClass[gsclss])):
#        if groundSnow0lowBranchClass[gsclss][gsdim][0][2]>=2:
#            print gsclss
#            #groundSnow0lowBranchClass2 = np.delete(groundSnow0lowBranchClass,gsclss[gsdim])

print 'now9'
#%% ploting
#fig4 = plt.figure(figsize=(20,15))
#ax4 = Axes3D(fig4)
#ax4.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1])#, dem_groundPoints[:, 2])

#for clss1 in allOpenClass:
#    for dim1 in range (len(clss1)):
#        ax4.scatter(clss1[dim1][0][0], clss1[dim1][0][1], clss1[dim1][0][2], c='orange')
#for clss in groundSnowClass:
#    for dim in range (len(clss)):
#        ax4.scatter(clss[dim][0][0], clss[dim][0][1], clss[dim][0][2], c='green')




































