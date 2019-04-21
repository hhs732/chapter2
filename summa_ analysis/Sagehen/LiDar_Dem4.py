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
#%% DEM file (.tif) reading (test)
#import gzip
#with gzip.open("lidardata\sagehen_demveg.tin.gz", 'rb') as f:
#     for line in f:        
#         print(line)
driver = gdal.GetDriverByName('GTiff')
filename = "lidardata\demvegTest.tif" #path to raster
demset = gdal.Open(filename)
lyr = gdal.GDALDEMProcessingOptions_swigregister(demset)

dem_hs_ds = gdal.DEMProcessing('', demset, 'hillshade', format='MEM')

band = demset.GetRasterBand(1)
BandType=gdal.GetDataTypeName(band.DataType)
minBand = band.GetMinimum()
maxBand = band.GetMaximum()
minMax = band.ComputeRasterMinMax(True)

elevation = band.ReadAsArray()
elevation[elevation == -9999.] = 1960

minelev = elevation[elevation > 1960] 
#plt.imshow(elevation, cmap='gist_earth')

cols = demset.RasterXSize
rows = demset.RasterYSize
print("Size is {} x {} x {}".format(demset.RasterXSize,demset.RasterYSize,demset.RasterCount))

geotransform = demset.GetGeoTransform()
x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
nrows, ncols = elevation.shape
x1 = x0 + dx * ncols
y1 = y0 + dy * nrows
extent=[x0, x1, y1, y0]
plt.imshow(elevation, cmap='gist_earth', extent=extent)
plt.savefig('lidardata\dem_bareearth_test.png')
#%% whole sagehen DEM
#filename = "lidardata\sagehenDem.tif" #path to raster
#demset1 = gdal.Open(filename)
#
#band1 = demset1.GetRasterBand(1)
#elevation1 = band1.ReadAsArray()
#elevation1[elevation1 >10000] = 1900
#
#x01, dx1, dxdy1, y01, dydx1, dy1 = demset1.GetGeoTransform()
#nrows1, ncols1 = elevation1.shape
#x11 = x01 + dx1 * ncols1
#y11 = y01 + dy1 * nrows1
#extent1=[x01, x11, y11, y01]
#plt.imshow(elevation1, cmap='gist_earth', extent=extent1)
#plt.savefig('lidardata\sagehenDem.png')
#%% creating centroid (ground points) from Dem files (test)
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

#%% difference between 2 years dems (snow on and snow off)
filename1 = "lidardata\demsnowTest.tif" #path to raster
demset1 = gdal.Open(filename1)
band1 = demset1.GetRasterBand(1)
elevation1 = band1.ReadAsArray()
elevation1[elevation1 == -9999.] = 1960

cols1 = demset1.RasterXSize
rows1 = demset1.RasterYSize

geotransform1 = demset1.GetGeoTransform()
x01, dx1, dxdy1, y01, dydx1, dy1 = demset1.GetGeoTransform()
nrows1, ncols1 = elevation1.shape
x11 = x01 + dx1 * ncols1
y11 = y01 + dy1 * nrows1
extent1=[x01, x11, y11, y01]
plt.imshow(elevation1, cmap='gist_earth', extent=extent1)
plt.savefig('lidardata\dem_bareearth_test3.png')

latitude1 =[]
for x in range (ncols1):
    latitude1.append(x+x01)
longitude1 = []
for y in range (nrows1):
    longitude1.append(y01-y)

latitude_rp1 = np.tile(latitude1, nrows1)
longitude_rp1 = np.repeat(longitude1, ncols1)
elevation_rp1 = np.reshape(elevation1,(nrows1*ncols1)).T
dem_groundPoints1 = np.vstack([latitude_rp1,longitude_rp1,elevation_rp1]).T

fig1 = plt.figure(figsize=(20,15))
ax1 = Axes3D(fig1)
ax1.scatter(dem_groundPoints1[:, 0], dem_groundPoints1[:, 1], dem_groundPoints1[:, 2])
ax1.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
plt.legend()
plt.savefig('lidardata\demscomaprison.png')
#%%
infileVeg = ls.file.File("lidardata\lasvegTest.las", mode="r")
# Grab all of the points from the file.
point_records = infileVeg.points
# Grab just the X dimension from the file, and scale it.
def scaled_x_dimension(las_file):
    x_dimension = las_file.X
    scale = las_file.header.scale[0]
    offset = las_file.header.offset[0]
    return(x_dimension*scale + offset)

scaled_x = scaled_x_dimension(infileVeg)
#%%# Find out what the point format looks like.
pointformat = infileVeg.point_format
for spec in infileVeg.point_format:
    print(spec.name)

#Lets take a look at the header also.
headerformat = infileVeg.header.header_format
for spec in headerformat:
    print(spec.name)
#%%# Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array
coordsVeg = np.vstack((infileVeg.x, infileVeg.y, infileVeg.z)).T
# calculating the nearest neighbors of a set of points, you might want to use a highly optimized package like FLANN 
datasetVeg = np.vstack([infileVeg.X, infileVeg.Y, infileVeg.Z]).T
minLat,maxLat = np.min(datasetVeg[:,0]),np.max(datasetVeg[:,0])
minLon, maxLon = np.min(datasetVeg[:,1]),np.max(datasetVeg[:,1])
#%%we’re interested only in the last return from each pulse in order to do ground detection. 
num_returns = infileVeg.num_returns
return_num = infileVeg.return_num
ground_points = infileVeg.points[num_returns == return_num]
# list and array of groundpoints from lidar file
groundPoints_ls = ground_points.tolist()
groundPoints_arr = []
for i in range (len(groundPoints_ls)):
    GPlist = np.array(groundPoints_ls[i])
    groundPoints_arr.append(GPlist[0,0:3])
groundPoints_arr = np.array(groundPoints_arr)
groundPoints_lidar = groundPoints_arr.copy()

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
testClass1 = []
classes_rplc = []#[x if np.size(x)!= 0 else [0,0,0] for x in classes1]
for cls in range (len (dem_groundPoints)):
    if len(classesVeg[cls])==0:
        nokhale = [np.array([0,0])]
        classes_rplc.append(nokhale)
    else: classes_rplc.append(classesVeg[cls])

    testClass1.append(abs(classes_rplc[cls]-dem_groundPoints[cls]))

failclass=[]
for xyidx in range (len(testClass1)):
    for xycl in range (len(testClass1[xyidx])):
        if ((testClass1[xyidx][xycl][0]>0.5) or (testClass1[xyidx][xycl][1]>0.5)):
            failclass.append(xyidx)
            break
#%% ploting
coordinatelidar = coordsVeg.copy()
datasetLidar = datasetVeg.copy()

fig = plt.figure(figsize=(20,15))
ax = Axes3D(fig)
ax.scatter(dem_groundPoints[:, 0], dem_groundPoints[:, 1], dem_groundPoints[:, 2])
ax.scatter(coordinatelidar[:, 0], coordinatelidar[:, 1], coordinatelidar[:, 2])
for flcl in failclass:
    #ax.scatter([x[0] for x in classes1[flcl]], [x[1] for x in classes1[flcl]])
    ax.scatter([x[0] for x in classesVeg[flcl]], [x[1] for x in classesVeg[flcl]], [x[2] for x in classesVeg[flcl]])
#ax.scatter([x[0] for x in classes1[25]], [x[1] for x in classes1[25]])
plt.savefig('lidardata\dem_lidar_classificationtest.png')
#%% vegtation classification from DEM2014 and las2014

vegClass2 = []
classes_rplc2 = []#[x if np.size(x)!= 0 else [0,0,0] for x in classes1]
for cls in range (len (dem_groundPoints)):
    if len(classesVeg[cls])==0:
        nokhale = [np.array([0,0,0])]
        classes_rplc2.append(nokhale)
    else: classes_rplc2.append(classesVeg[cls])
    vegClass2.append(classes_rplc2[cls]-dem_groundPoints[cls])
    
vegClass1 = []
classes_rplc = []#[x if np.size(x)!= 0 else [0,0,0] for x in classes1]
for cls in range (len (dem_groundPoints)):
    if len(classesVeg[cls])==0:
        nokhale = [np.array([0,0,0])]
        classes_rplc.append(nokhale)
    else: classes_rplc.append(classesVeg[cls])
    test = []
    for tst in range(len(classes_rplc[cls])):
        height = classes_rplc[cls][tst][2]-dem_groundPoints[cls][2]
        test.append(np.vstack([classes_rplc[cls][tst][0],classes_rplc[cls][tst][1],height]).T)
    vegClass1.append(test)

#all tree classification
allTreeClass = []
numTreeClass = []
for vgcl in range (len (dem_groundPoints)):
    for pnt in range (len (vegClass1[vgcl])):
        if vegClass1[vgcl][pnt][2]>2:
            #print 'hooorrraaa tree'
            numTreeClass.append(vgcl)
            allTreeClass.append(vegClass1[vgcl])
            break
        
negVegTreeClass = []
numNegVegClass = []
for tcl in range (len (allTreeClass)):
    for pn in range (len (allTreeClass[tcl])):
        if allTreeClass[tcl][pn][2]<-0.15:
            numNegVegClass.append(tcl)
            negVegTreeClass.append(allTreeClass[tcl])
            break        
# trees with low branches
lowVegTreeClass = []
numlowVegClass = []
for tcl in range (len (allTreeClass)):
    for pn in range (len (allTreeClass[tcl])):
        if allTreeClass[tcl][pn][2]>0.15 and allTreeClass[tcl][pn][2]<2:
            numlowVegClass.append(tcl)
            lowVegTreeClass.append(allTreeClass[tcl])
            break
# trees with no low blanches
allTreeClassNum = list(np.arange(len(allTreeClass)))
numNolowVegClass = list(set(allTreeClassNum)-set(numlowVegClass))
nolowVegTreeClass = []
for indx in numNolowVegClass:
    nolowVegTreeClass.append(allTreeClass[indx])    
# open space (no trees, no return between 0.15 -2)
# all low veg
allLowVegClass = []
numallLowVegClass = []
for lvgcl in range (len (dem_groundPoints)):
    for crdnt in range (len (vegClass1[lvgcl])):
        if vegClass1[lvgcl][crdnt][2]<2 and vegClass1[lvgcl][crdnt][2]>0.15:
            numallLowVegClass.append(lvgcl)
            allLowVegClass.append(vegClass1[lvgcl])
            break
#open places
numNotopen = list(set(numallLowVegClass).union(set(numTreeClass)))
allClassNum = list(np.arange(len(vegClass1)))
numopen =  list(set(allClassNum)-set(numNotopen))
allOpenClass = []
for idop in numopen:
    allOpenClass.append(vegClass1[idop])

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
testClass2 = []
classes_rplc2 = []#[x if np.size(x)!= 0 else [0,0,0] for x in classes1]
for cls in range (len (dem_groundPoints)):
    if len(classesnow[cls])==0:
        nokhale = [np.array([0,0,0])]
        classes_rplc2.append(nokhale)
    else: classes_rplc2.append(classesnow[cls])

    testClass2.append(abs(classes_rplc2[cls]-dem_groundPoints[cls]))

failclass2=[]
for xyidx in range (len(testClass2)):
    for xycl in range (len(testClass2[xyidx])):
        if ((testClass2[xyidx][xycl][0]>0.5) or (testClass2[xyidx][xycl][1]>0.5)):
            failclass2.append(xyidx)
            break
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
vegsnowClass = []
classes_rplc2 = []#[x if np.size(x)!= 0 else [0,0,0] for x in classes1]
for cls in range (len (dem_groundPoints)):
    if len(classesnow[cls])==0:
        nokhale = [np.array([0,0,0])]
        classes_rplc2.append(nokhale)
    else: classes_rplc.append(classesnow[cls])

    vegsnowClass.append(classes_rplc2[cls]-dem_groundPoints[cls])




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




































