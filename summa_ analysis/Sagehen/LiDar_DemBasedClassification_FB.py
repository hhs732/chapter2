import laspy as ls
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from mpl_toolkits.mplot3d import Axes3D

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
 
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows
    extent=[x0, x1, y1, y0]
    
    plt.imshow(elevation, cmap='gist_earth', extent=extent)
    plt.savefig(pathName)
    
    return x0, y0, elevation, nrows, ncols

def creatingCentroidGroundpointsFromDem(tiffFilename,elevationMissNo,pathNameforDemImage):
    demset = gdal.Open(tiffFilename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
 
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = elevation.shape
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows
    extent=[x0, x1, y1, y0]
    
    plt.imshow(elevation, cmap='gist_earth', extent=extent)
    plt.savefig(pathNameforDemImage)
    
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
#%% DEM file (.tif) reading (test)
#filename = "lidardata\demvegTest.tif" #path to raster
#demset = gdal.Open(filename)
#lyr = gdal.GDALDEMProcessingOptions_swigregister(demset)
#
#dem_hs_ds = gdal.DEMProcessing('', demset, 'hillshade', format='MEM')
#
#band = demset.GetRasterBand(1)
#BandType=gdal.GetDataTypeName(band.DataType)
#minBand = band.GetMinimum()
#maxBand = band.GetMaximum()
#minMax = band.ComputeRasterMinMax(True)
#
#cols = demset.RasterXSize
#rows = demset.RasterYSize
#print("Size is {} x {} x {}".format(demset.RasterXSize,demset.RasterYSize,demset.RasterCount))

#elevationMissNo = 1960.
#pathName = 'lidardata\dem_bareearth_test.png'
#x0, y0, elevation, nrows, ncols = readPlotDEM(filename,elevationMissNo,pathName)
#%% LiDar Data reading
#infileVeg = ls.file.File("lidardata\lasvegTest.las", mode="r")
## Grab all of the points from the file.
#point_records = infileVeg.points
#
#scaled_x = scaled_x_dimension(infileVeg)
#
## Find out what the point format looks like.
#pointformat = infileVeg.point_format
#for spec in infileVeg.point_format:
#    print(spec.name)
#
##Lets take a look at the header also.
#headerformat = infileVeg.header.header_format
#for spec in headerformat:
#    print(spec.name)
#%%we’re interested only in the last return from each pulse in order to do ground detection. 
#num_returns = infileVeg.num_returns
#return_num = infileVeg.return_num
#ground_points = infileVeg.points[num_returns == return_num]
## list and array of groundpoints from lidar file
#groundPoints_ls = ground_points.tolist()
#groundPoints_arr = []
#for i in range (len(groundPoints_ls)):
#    GPlist = np.array(groundPoints_ls[i])
#    groundPoints_arr.append(GPlist[0,0:3])
#groundPoints_arr = np.array(groundPoints_arr)
#groundPoints_lidar = groundPoints_arr.copy()
#%% DEM file (.tif) reading and creating centroid (ground points) from Dem files (test)
filename = "lidardata\demvegTest.tif" #path to raster
elevationMissNo = 1960.
pathName = 'lidardata\dem_bareearth_test.png'
dem_groundPoints = creatingCentroidGroundpointsFromDem(filename,elevationMissNo,pathName)

#%% difference between 2 years dems (snow on and snow off)
filenameSnow = "lidardata\demsnowTest.tif" #path to raster
elevationMissNoSnow = 1960
outImagePathName = 'lidardata\dem_bareearth_test3.png'
dem_groundPointsSnow = creatingCentroidGroundpointsFromDem(filenameSnow,elevationMissNoSnow,outImagePathName)

fig1 = plt.figure(figsize=(20,15))
ax1 = Axes3D(fig1)
ax1.scatter(dem_groundPointsSnow[:, 0], dem_groundPointsSnow[:, 1], dem_groundPointsSnow[:, 2])
ax1.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
plt.legend()
plt.savefig('lidardata\demscomaprison.png')

#%%# LiDar Data reading and Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array
infileVeg = ls.file.File("lidardata\lasvegTest.las", mode="r")

coordsVeg = np.vstack((infileVeg.x, infileVeg.y, infileVeg.z)).T
# calculating the nearest neighbors of a set of points, you might want to use a highly optimized package like FLANN 
datasetVeg = np.vstack([infileVeg.X, infileVeg.Y, infileVeg.Z]).T
#minLat,maxLat = np.min(datasetVeg[:,0]),np.max(datasetVeg[:,0])
#minLon, maxLon = np.min(datasetVeg[:,1]),np.max(datasetVeg[:,1])
#%% classification with new grountpoints from dem file
centroids_new=dem_groundPoints[:,0:2]  
k1 = np.size(dem_groundPoints[:,0])
# instantiate a class
clf1 = K_Means(numOfClusters=k1,init_centroids=centroids_new)
# fit kmean class to data
clf1.fit(coordsVeg,0,2)
# get classification 
classesVeg = clf1.classifications
#%% test classification
upGroundPointsVeg, classes_rplcVeg, pureVegClass = lidarDiffGrndPoints(classesVeg,dem_groundPoints)
failclassVeg = classificationTest(pureVegClass)

#%% ploting
fig = plt.figure(figsize=(20,15))
ax = Axes3D(fig)
#ax.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
#ax.scatter(coordsVeg[:, 0], coordsVeg[:, 1], coordsVeg[:, 2])
for flcl in failclassVeg:
    #ax.scatter([x[0] for x in classes1[flcl]], [x[1] for x in classes1[flcl]])
    ax.scatter([x[0] for x in classesVeg[flcl]], [x[1] for x in classesVeg[flcl]], [x[2] for x in classesVeg[flcl]])
#ax.scatter([x[0] for x in classes1[25]], [x[1] for x in classes1[25]])
plt.savefig('lidardata\dem_lidar_classificationtest.png')
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
#%%#%% snow on lidar files
infileSnow = ls.file.File("lidardata\lassnowTest.las", mode="r")
coordsnow = np.vstack((infileSnow.x, infileSnow.y, infileSnow.z)).T
datasetSnow = np.vstack([infileSnow.X, infileSnow.Y, infileSnow.Z]).T

#%% classification with new grountpoints from dem file
clfs = K_Means(numOfClusters=k1,init_centroids=centroids_new)
# fit kmean class to data
clfs.fit(coordsnow,0,2)
# get classification 
classesnow = clfs.classifications
#%% test classification
vegSnowClass, classes_rplcVS, pureVegsnowClass = lidarDiffGrndPoints(classesnow,dem_groundPoints)
failclassSnow = classificationTest(pureVegsnowClass)
vegSnowClass2, vegSnowNumClass = defineSpecificClassGreater (vegSnowClass, -1)
#%% ploting
fig3 = plt.figure(figsize=(20,15))
ax3 = Axes3D(fig3)
ax3.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
ax3.scatter(coordsnow[:, 0], coordsnow[:, 1], coordsnow[:, 2])
ax3.scatter(coordsVeg[:, 0], coordsVeg[:, 1], coordsVeg[:, 2])
ax3.legend()
#for flcl in failclass2:
#    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
   plt.savefig('lidardata\dem_lidar_snow&veg.png')
#%% raw classification
rawClass = infileSnow.raw_classification
return_arr = infileSnow.num_returns
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






















#%%
#Y
#Z
#intensity
#flag_byte
#raw_classification
#scan_angle_rank
#user_data
#pt_src_id
#gps_time
#file_sig ???????????????????
#file_source_id
#global_encoding
#proj_id_1 ??????????????
#proj_id_2   ????????????????
#proj_id_3   ?????????
#proj_id_4    ???????????/
#version_major
#version_minor
#system_id
#software_id
#created_day
#created_year
#header_size
#data_offset
#num_variable_len_recs
#data_format_id
#data_record_length
#point_records_count
#point_return_count
#x_scale
#y_scale
#z_scale
#x_offset
#y_offset
#z_offset
#x_max
#x_min
#y_max
#y_min
#z_max
#z_min




































