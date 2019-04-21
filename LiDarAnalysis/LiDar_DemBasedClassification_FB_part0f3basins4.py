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

elevationMissNoS0f = 1900.

def readPlotDEM(filename,elevationMissNo):#,pathName
    demset = gdal.Open(filename)
    band = demset.GetRasterBand(1)
    elevation = band.ReadAsArray()
    elevation[elevation == -9999.] = elevationMissNo
    elevation[elevation > 100000] = elevationMissNo
 
    #x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    #nrows, ncols = elevation.shape
    #x1 = x0 + dx * ncols
    #y1 = y0 + dy * nrows
    #extent=[x0, x1, y1, y0]
    
    #plt.figure(figsize=(30,20))
    #plt.imshow(elevation, cmap='gist_earth', extent=extent)
    #plt.savefig(pathName)
    
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
#%%example
A= pd.DataFrame(np.array([[2,22,1],[3,23,1],[1,21,0],[4,24,9]]),columns=['x','y','z'])

B= pd.DataFrame(np.array([[1,21,2],[2,22,3],[3,23,3],[10,200,3],[4,24,2]]),columns=['x','y','z'])
#B= np.array([[1,21,2],
#            [2,22,3],
#            [10,200,3],
#            [3,23,3],
#            [4,24,2]
#            ])
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

print (C)


#%% dem snow off (veg) for vcm
minXVcm = 362000.
maxXVcm = 362010.
minYVcm = 3972985.
maxYVcm = 3973000.

filenameS0fvcm = 'W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/vcm_snw0ff_dem.tif' #path to raster
#pathNameS0fvcm = 'W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/dem_snow0ff_vcm.png'
elevationVegvcm = readPlotDEM(filenameS0fvcm,elevationMissNoS0f)#,pathNameS0fvcm
dem_groundPointsVegVcm = creatingCentroidGroundpointsFromDem(filenameS0fvcm,elevationMissNoS0f)#,pathNameS)
dem_groundPointsVegVcm_df0 = pd.DataFrame(dem_groundPointsVegVcm,columns=['x','y','z'])
dem_groundPointsVegVcm_df = pd.concat([dem_groundPointsVegVcm_df0[['x','y']].astype(int),dem_groundPointsVegVcm_df0['z']], axis=1)
dem_groundPointsVegVcm_df.sort_values(by=['x','y'],inplace=True)
dem_groundPointsVegVcm_df.index=np.arange(0,len(dem_groundPointsVegVcm_df))
   
dem_groundPointsVegVcm_int = dem_groundPointsVegVcm_df[(dem_groundPointsVegVcm_df['x'] >= minXVcm) & (dem_groundPointsVegVcm_df['x'] <= maxXVcm)]
dem_groundPointsVegVcm_sp0 = dem_groundPointsVegVcm_int[(dem_groundPointsVegVcm_int['y'] >= minYVcm) & (dem_groundPointsVegVcm_int['y'] <= maxYVcm)]
dem_groundPointsVegVcm_sp = dem_groundPointsVegVcm_sp0.values
dem_groundPointsVegVcm_intg = pd.concat([dem_groundPointsVegVcm_sp0[['x','y']].astype(int),dem_groundPointsVegVcm_sp0['z']], axis=1)
dem_groundPointsVegVcm_intg.sort_values(by=['x','y'],inplace=True)
dem_groundPointsVegVcm_intg.index=np.arange(0,len(dem_groundPointsVegVcm_intg))

#dem snow on (snw) for vcm
#filenameS0nVcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/vcm_snw0n_dem.tif" #path to raster
#pathNameS0nVcm = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/dem_snow0n_vcm.png"
#elevationSnowVcm = readPlotDEM(filenameS0nVcm,elevationMissNoS0f,pathNameS0nVcm)
#dem_groundPointsSnowVcm = creatingCentroidGroundpointsFromDem(filenameS0nVcm,elevationMissNoS0f)#,pathNameS)

#las file snow off (snw) for vcm
infileVegVcm = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/jemezSnow0ffpoints.las", mode="r")
coordsVegVcm = np.vstack((infileVegVcm.x, infileVegVcm.y, infileVegVcm.z)).T
coordsVegVcm_df0 = pd.DataFrame(coordsVegVcm,columns=['x','y','z'])
coordsVegVcm_df = pd.concat([coordsVegVcm_df0[['x','y']].astype(int),coordsVegVcm_df0['z']], axis=1)
coordsVegVcm_df.sort_values(by=['x','y'],inplace=True)
coordsVegVcm_df.index=np.arange(0,len(coordsVegVcm_df))

minXVcm2 = 362000.
maxXVcm2 = 362012.
minYVcm2 = 3972990.
maxYVcm2 = 3973000.

coordsVegVcm_df_int = coordsVegVcm_df[(coordsVegVcm_df['x'] >= minXVcm2) & (coordsVegVcm_df['x'] <= maxXVcm2)]
coordsVegVcm_sp0 = coordsVegVcm_df_int[(coordsVegVcm_df_int['y'] >= minYVcm2) & (coordsVegVcm_df_int['y'] <= maxYVcm2)]
coordsVegVcm_sp = coordsVegVcm_sp0.values
coordsVegVcm_intg = pd.concat([coordsVegVcm_sp0[['x','y']].astype(int),coordsVegVcm_sp0['z']], axis=1)
coordsVegVcm_intg.sort_values(by=['x','y'],inplace=True)
coordsVegVcm_intg.index=np.arange(0,len(coordsVegVcm_intg))
#
##las file snow on (snw) for vcm
#infileSnwVcm = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/JemezVcm/jemezSnow0npoints.las", mode="r")
#coordsSnwVcm = np.vstack((infileSnwVcm.x, infileSnwVcm.y, infileSnwVcm.z)).T
#coordsSnwVcm_df = pd.DataFrame(coordsSnwVcm,columns=['x','y','z'])
#
#coordsnow0nVcm_int = coordsSnwVcm_df[(coordsSnwVcm_df['x'] >= minXVcm) & (coordsSnwVcm_df['x'] <= maxXVcm)]
#coordSnwVcm_sp0 = coordsnow0nVcm_int[(coordsnow0nVcm_int['y'] >= minYVcm) & (coordsnow0nVcm_int['y'] <= maxYVcm)]
#coordsSnwVcm_sp = coordSnwVcm_sp0.values
#
#coordsSnwVcm_intg = pd.concat([coordSnwVcm_sp0[['x','y']].astype(int),coordSnwVcm_sp0['z']], axis=1)
#coordsSnwVcm_intg.sort_values(by=['x','y'],inplace=True)
#coordsSnwVcm_intg.index=np.arange(0,len(coordsSnwVcm_intg))

#%%
#upGroundPointsVegVcm2=coordsVegVcm_intg.copy()
#R=len(coordsVegVcm_intg)
#indx = []
upGroundPointsVegVcm = []
#print dem_groundPointsVegVcm_df['y'][0]
for row in coordsVegVcm_intg.itertuples():
    indexFinl = np.where(((dem_groundPointsVegVcm_intg['x']==row.x) & (dem_groundPointsVegVcm_intg['y']==row.y)))[0]#
    if (np.size(indexFinl))>0:
        upG = row.z-dem_groundPointsVegVcm_intg['z'][indexFinl[0]]
        upGroundPointsVegVcm.append([row.x,row.y,upG])#
    #else: upGroundPointsVegVcm.append([row.x,row.y,-999])
     


#for indx, row in coordsVegVcm_intg[['x','y']].iterrows():
#    R=np.where((dem_groundPointsVegVcm_intg['x']==row.x) & (dem_groundPointsVegVcm_intg['y']==row.y))[0]
#    #test.append(R)
#    dem_groundPointsVegVcm_intg[['x','y']].loc[R[0]]
#    upGroundPointsVegVcm2['z'][indx] = coordsVegVcm_intg['z'][indx]-dem_groundPointsVegVcm_intg['z'][R]
#    print len(test)


#for i in range(R):
#    row=coordsVegVcm_intg[['x','y']].iloc[i]
#    
#    r=np.where(dem_groundPointsVegVcm_intg[['x','y']]==row)[0]
#    print (r)
#    if len(r)>0:
#        #if R[0]==R[1]:
#            #print ("row %i of B is the same as row %i of A" %(i,R[0]))
#        upGroundPointsVegVcm2['z'][i] = coordsVegVcm_intg['z'][i]-dem_groundPointsVegVcm_intg['z'][r[0]]
#        indx.append(i)
#    else:
#        upGroundPointsVegVcm2['z'][i]=-9999
        
#(coordsVegVcm_intg['x'].isin(dem_groundPointsVegVcm_intg['x']))&(coordsVegVcm_intg['y'].isin(dem_groundPointsVegVcm_intg['y']))

#%% classification with grountpoints from veg dem file :VCM
#centroids_newVcm=dem_groundPointsVegVcm_sp[:,0:2]  
#kVcm = np.size(dem_groundPointsVegVcm_sp[:,0])
## instantiate a class
#clf1Vcm = K_Means(numOfClusters=kVcm,init_centroids=centroids_newVcm)
## fit kmean class to data
#clf1Vcm.fit(coordsVegVcm_sp,0,2)
## get classification 
#classesVegVcm = clf1Vcm.classifications
#
#upGroundPointsVegVcm, classes_rplcVegVcm, pureVegClassVcm = lidarDiffGrndPoints(classesVegVcm,dem_groundPointsVegVcm_sp)
#
##%% classification of snow las file with veg_grountpoints
#clfsVcm = K_Means(numOfClusters=kVcm,init_centroids=centroids_newVcm)
## fit kmean class to data
#clfsVcm.fit(coordsSnwVcm_sp,0,2)
## get classification 
#classeSnwVcm = clfsVcm.classifications
#
#upGroundPointsSnwVegVcm, classes_rplcVSVcm, pureVegsnowClassVcm = lidarDiffGrndPoints(classeSnwVcm,dem_groundPointsVegVcm_sp)
##%% vegtation classification from DEM2014 and las2014
#vegClassVcm = upGroundPointsVegVcm[:]
##all tree classification
#allTreeClassVcm, treeNumClassVcm = defineSpecificClassGreater (vegClassVcm, 2)
#        
#negVegClassVcm, negVegNumClassVcm = defineSpecificClassLess (vegClassVcm, 0)
#
## trees with low branches
#lowVegTreeClassVcm, lowVegNumClassVcm = defineLowVegClass(allTreeClassVcm)
#
## trees with no low blanches
## "*************tall canopy no snow class*****************"
#nolowVegTreeClassVcm, nolowVegTreeNumClassVcm = differenceBetwee2classes (allTreeClassVcm,lowVegNumClassVcm)
#
## open space (no trees, no return between 0.15 to 2)
## all low veg
#allLowVegClassVcm, allLowVegNumClassVcm = defineLowVegClass(vegClassVcm)
#
##open places
#notOpenNumClassVcm = list(set(allLowVegNumClassVcm).union(set(treeNumClassVcm)))
## "******************open no snow class*******************"
#allOpenClassVcm, allOpenNumClassVcm = differenceBetwee2classes (vegClassVcm,notOpenNumClassVcm)
#
##%% snow classification
#vegSnowClassVcm, vegSnowNumClassVcm = defineSpecificClassGreater (upGroundPointsSnwVegVcm, -1)
#
##no snow on the ground
#nosnowClassVcm, noSnowNumClassVcm = defineSpecificClassLess (vegSnowClassVcm, 0.15)
##snow on the ground or on the trees
#allSnowClassVcm, allSnowNumClassVcm = differenceBetwee2classes (vegSnowClassVcm,noSnowNumClassVcm)
## "******************tall canopy snow class*******************"
##snow on the tall canopy >2m
#treeSnowClassVcm, treeSnowNumClassVcm = defineSpecificClassGreater (allSnowClassVcm, 2)
## "******************open snow class**************************"
##snow on the ground
#groundSnowClassVcm, groundSnowNumClassVcm = differenceBetwee2classes (allSnowClassVcm,treeSnowNumClassVcm)
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
#%%
##%%
##%%
##%%
##%%
##%% dem snow off (veg) for NR1 Niwot
#minXNr1 = 454900.
#maxXNr1 = 455000.
#minYNr1 = 4430500.
#maxYNr1 = 4430600.
#
#filenameS0fNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/Niwot/nr1_snw0ff_dem.tif" #path to raster
#pathNameS0fNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/Niwot/dem_snow0ff_Nr1.png"
#elevationVegNr1 = readPlotDEM(filenameS0fNr1,elevationMissNoS0f,pathNameS0fNr1)
#dem_groundPointsVegNr1 = creatingCentroidGroundpointsFromDem(filenameS0fNr1,elevationMissNoS0f)#,pathNameS)
#
#dem_groundPointsVegNr1_df = pd.DataFrame(dem_groundPointsVegNr1,columns=['x','y','z'])
#dem_groundPointsVegNr1_int = dem_groundPointsVegNr1_df[(dem_groundPointsVegNr1_df['x'] >= minXNr1) & (dem_groundPointsVegNr1_df['x'] <= maxXNr1)]
#dem_groundPointsVegNr1_sp0 = dem_groundPointsVegNr1_int[(dem_groundPointsVegNr1_int['y'] >= minYNr1) & (dem_groundPointsVegNr1_int['y'] <= maxYNr1)]
#dem_groundPointsVegNr1_sp = dem_groundPointsVegNr1_sp0.values
#
## dem snow on (snw) for NR1
##filenameS0nNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/Niwot/nr1_snw0n9_dem.tif" #path to raster
##pathNameS0nNr1 = "W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/Niwot/dem_snow0n_Nr1.png"
##elevationSnowNr1 = readPlotDEM(filenameS0nNr1,elevationMissNoS0f,pathNameS0nNr1)
##dem_groundPointsSnowNr1 = creatingCentroidGroundpointsFromDem(filenameS0nNr1,elevationMissNoS0f)#,pathNameS)
#
##las file snow off (snw) for NR1
#infileVegNr1 = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/Niwot/niwotSnow0ffpoints.las", mode="r")
#coordsVegNr1 = np.vstack((infileVegNr1.x, infileVegNr1.y, infileVegNr1.z)).T
#
#coordsVegNr1_df = pd.DataFrame(coordsVegNr1,columns=['x','y','z'])
#coordsVegNr1_df_int = coordsVegNr1_df[(coordsVegNr1_df['x'] >= minXNr1) & (coordsVegNr1_df['x'] <= maxXNr1)]
#coordsVegNr1_sp0 = coordsVegNr1_df_int[(coordsVegNr1_df_int['y'] >= minYNr1) & (coordsVegNr1_df_int['y'] <= maxYNr1)]
#coordsVegNr1_sp = coordsVegNr1_sp0.values
#
##las file snow on (snw) for NR1
#infileSnwNr1 = ls.file.File("W:/_users/hsafa/Chapter2_snow_forest/LIDAR_analysis/Niwot/niwotSnowOn09may2010points.las", mode="r")
#coordsSnwNr1 = np.vstack((infileSnwNr1.x, infileSnwNr1.y, infileSnwNr1.z)).T
#
#coordsSnwNr1_df = pd.DataFrame(coordsSnwNr1,columns=['x','y','z'])
#coordsnow0nNr1_int = coordsSnwNr1_df[(coordsSnwNr1_df['x'] >= minXNr1) & (coordsSnwNr1_df['x'] <= maxXNr1)]
#coordSnwNr10 = coordsnow0nNr1_int[(coordsnow0nNr1_int['y'] >= minYNr1) & (coordsnow0nNr1_int['y'] <= maxYNr1)]
#coordsSnwNr1_sp = coordSnwNr10.values
##%% classification with grountpoints from veg dem file :NR1
#centroids_newNr1=dem_groundPointsVegNr1_sp[:,0:2]  
#kNr1 = np.size(dem_groundPointsVegNr1_sp[:,0])
## instantiate a class
#clf1Nr1 = K_Means(numOfClusters=kNr1,init_centroids=centroids_newNr1)
## fit kmean class to data
#clf1Nr1.fit(coordsVegNr1_sp,0,2)
## get classification 
#classesVegNr1 = clf1Nr1.classifications
#
#upGroundPointsVegNr1, classes_rplcVegNr1, pureVegClassNr1 = lidarDiffGrndPoints(classesVegNr1,dem_groundPointsVegNr1_sp)
#
##%% classification of snow las file with veg_grountpoints
#clfsNr1 = K_Means(numOfClusters=kNr1,init_centroids=centroids_newNr1)
## fit kmean class to data
#clfsNr1.fit(coordsSnwNr1_sp,0,2)
## get classification 
#classeSnwNr1 = clfsNr1.classifications
#
#upGroundPointsSnwVegNr1, classes_rplcVSNr1, pureVegsnowClassNr1 = lidarDiffGrndPoints(classeSnwNr1,dem_groundPointsVegNr1_sp)
#
##%% vegtation classification from DEM2014 and las2014
#vegClassNr1 = upGroundPointsVegNr1[:]
##all tree classification
#allTreeClassNr1, treeNumClassNr1 = defineSpecificClassGreater (vegClassNr1, 2)
#        
#negVegClassNr1, negVegNumClassNr1 = defineSpecificClassLess (vegClassNr1, 0)
#
## trees with low branches
#lowVegTreeClassNr1, lowVegNumClassNr1 = defineLowVegClass(allTreeClassNr1)
#
## trees with no low blanches
## "*************tall canopy no snow class*****************"
#nolowVegTreeClassNr1, nolowVegTreeNumClassNr1 = differenceBetwee2classes (allTreeClassNr1,lowVegNumClassNr1)
#
## open space (no trees, no return between 0.15 to 2)
## all low veg
#allLowVegClassNr1, allLowVegNumClassNr1 = defineLowVegClass(vegClassNr1)
#
##open places
#notOpenNumClassNr1 = list(set(allLowVegNumClassNr1).union(set(treeNumClassNr1)))
## "******************open no snow class*******************"
#allOpenClassNr1, allOpenNumClassNr1 = differenceBetwee2classes (vegClassNr1,notOpenNumClassNr1)
#
##%% snow classification
#vegSnowClassNr1, vegSnowNumClassNr1 = defineSpecificClassGreater (upGroundPointsSnwVegNr1, -1)
#
##no snow on the ground
#nosnowClassNr1, noSnowNumClassNr1 = defineSpecificClassLess (vegSnowClassNr1, 0.15)
##snow on the ground or on the trees
#allSnowClassNr1, allSnowNumClassNr1 = differenceBetwee2classes (vegSnowClassNr1,noSnowNumClassNr1)
## "******************tall canopy snow class*******************"
##snow on the tall canopy >2m
#treeSnowClassNr1, treeSnowNumClassNr1 = defineSpecificClassGreater (allSnowClassNr1, 2)
## "******************open snow class**************************"
##snow on the ground
#groundSnowClassNr1, groundSnowNumClassNr1 = differenceBetwee2classes (allSnowClassNr1,treeSnowNumClassNr1)
##%%
#figNr1 = plt.figure(figsize=(20,15))
#axNr1 = Axes3D(figNr1)
#axNr1.scatter(dem_groundPointsVegNr1_sp[:, 0], dem_groundPointsVegNr1_sp[:, 1], dem_groundPointsVegNr1_sp[:, 2])
#axNr1.scatter(coordsSnwNr1_sp[:, 0], coordsSnwNr1_sp[:, 1], coordsSnwNr1_sp[:, 2])
##axNr1.scatter(coordsVegNr1_sp[:, 0], coordsVegNr1_sp[:, 1], coordsVegNr1_sp[:, 2])
#axNr1.legend()
#plt.title('snow Lidar data sagehen ',fontsize=30)
#
##for flcl in failclass2:
##    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
#plt.savefig("H:/Chapter2_snow_forest/LIDAR_analysis/Niwot/dem_groundPointsVegT1&coordsVegNr1.png")
#
##%%
##%%
##%%
##%%
##%%
##%% DEM file (.tif) reading and creating centroid (ground points) from Dem files for sagehen
#
#filenameS0fSc = "W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\Sagehen\sc_veg_dem.tif" #path to raster
#pathNameS0fSc = "W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\Sagehen\dem_snow0ff_sc.png"
#elevationVegSc = readPlotDEM(filenameS0fSc,elevationMissNoS0f,pathNameS0fSc)
#dem_groundPointsVegSc = creatingCentroidGroundpointsFromDem(filenameS0fSc,elevationMissNoS0f)#,pathNameS)
#dem_groundPointsVegSc_df = pd.DataFrame(dem_groundPointsVegSc,columns=['x','y','z'])
#dem_groundPointsVegSc_intg = pd.concat([dem_groundPointsVegSc_df[['x','y']].astype(int),dem_groundPointsVegSc_df['z']],axis=1)
#dem_groundPointsVegSc_intg.sort_values(by=['x'])
#
## LiDar Data reading and Grab the scaled x, y, and z dimensions and stick them together in an nx3 numpy array
#infileVegSc = ls.file.File("W:\_users\hsafa\Chapter2_snow_forest\LIDAR_analysis\Sagehen\sc_snw0ff_points.las", mode="r")
#coordsVegSc = np.vstack((infileVegSc.x, infileVegSc.y, infileVegSc.z)).T
#coordsVegSc_df = pd.DataFrame(coordsVegSc,columns=['x','y','z'])
#coordsVegSc_intg = pd.concat([coordsVegSc_df[['x','y']].astype(int),coordsVegSc_df['z']],axis=1)
#coordsVegSc_intg.sort_values(by=['x','y'])
#
##lidar data for whole sagehen, snow on 2016
#from LiDar_DemBasedClassification_scSnow0n import coordsnow0nSc26M
#coordsnow0nSc26M_df = pd.DataFrame(coordsnow0nSc26M,columns=['x','y','z'])
#coordsnow0nSc26M_intg = pd.concat([coordsnow0nSc26M_df[['x','y']].astype(int),coordsnow0nSc26M_df['z']],axis=1)
#coordsnow0nSc26M_intg.sort_values(by=['x','y'])
#
#from LiDar_DemBasedClassification_scSnow0n import coordsnow0nSc18A
#coordsnow0nSc18M_df = pd.DataFrame(coordsnow0nSc18A,columns=['x','y','z'])
#coordsnow0nSc18M_intg = pd.concat([coordsnow0nSc18M_df[['x','y']].astype(int),coordsnow0nSc18M_df['z']],axis=1)
#coordsnow0nSc18M_intg.sort_values(by=['x','y'])
##%% some part of data
#minXSc = 736400.
#maxXSc = 736500.
#minYSc = 4368100.
#maxYSc = 4368200.
#
#dem_groundPointsVegSc_int = dem_groundPointsVegSc_df[(dem_groundPointsVegSc_df['x'] >= minXSc) & (dem_groundPointsVegSc_df['x'] <= maxXSc)]
#dem_groundPointsVegSc_sp0 = dem_groundPointsVegSc_int[(dem_groundPointsVegSc_int['y'] >= minYSc) & (dem_groundPointsVegSc_int['y'] <= maxYSc)]
#dem_groundPointsVegSc_sp = dem_groundPointsVegSc_sp0.values
#
#coordsVegSc_df_int = coordsVegSc_df[(coordsVegSc_df['x'] >= minXSc) & (coordsVegSc_df['x'] <= maxXSc)]
#coordsVegSc_sp0 = coordsVegSc_df_int[(coordsVegSc_df_int['y'] >= minYSc) & (coordsVegSc_df_int['y'] <= maxYSc)]
#coordsVegSc_sp = coordsVegSc_sp0.values
#
#lidar data for sagehen, snow on 2016, 26 March
#coordsnow0nSc26M_int = coordsnow0nSc26M[(coordsnow0nSc26M['x'] >= minXSc) & (coordsnow0nSc26M['x'] <= maxXSc)]
#coordSnwSc26M0 = coordsnow0nSc26M_int[(coordsnow0nSc26M_int['y'] >= minYSc) & (coordsnow0nSc26M_int['y'] <= maxYSc)]
#coordsSnwSc26M = coordSnwSc26M0.values
#
#lidar data for sagehen, snow on 2016, 18 May
#coordsnow0nSc18M_int = coordsnow0nSc18A[(coordsnow0nSc18A['x'] >= minXSc) & (coordsnow0nSc18A['x'] <= maxXSc)]
#coordSnwSc18M0 = coordsnow0nSc18M_int[(coordsnow0nSc18M_int['y'] >= minYSc) & (coordsnow0nSc18M_int['y'] <= maxYSc)]
#coordsSnwSc18M = coordSnwSc18M0.values
##%% classification with veggrountpoints from veg dem file :T1
#centroids_newT1=dem_groundPointsVegSc_sp[:,0:2]  
#kT1 = np.size(dem_groundPointsVegSc_sp[:,0])
## instantiate a class
#clf1T1 = K_Means(numOfClusters=kT1,init_centroids=centroids_newT1)
## fit kmean class to data
#clf1T1.fit(coordsVegSc_sp,0,2)
## get classification 
#classesVegSc = clf1T1.classifications
#
#upGroundPointsVegSc, classes_rplcVegSc, pureVegClassSc = lidarDiffGrndPoints(classesVegSc,dem_groundPointsVegSc_sp)
#
##%% classification of snow las file with veg_grountpoints
#clfsT1 = K_Means(numOfClusters=kT1,init_centroids=centroids_newT1)
## classification for sagehen 26 March
#clfsT1.fit(coordsSnwSc26M,0,2)
#classeSnwSc26M = clfsT1.classifications
#upGroundPointsSnwVegSc26M, classes_rplcVSSc26M, pureVegsnowClassSc26M = lidarDiffGrndPoints(classeSnwSc26M,dem_groundPointsVegSc_sp)
## classification for sagehen 18 May
#clfsT1.fit(coordsSnwSc18M,0,2)
#classeSnwSc18M = clfsT1.classifications
#upGroundPointsSnwVegSc18M, classes_rplcVSSc18M, pureVegsnowClassSc18M = lidarDiffGrndPoints(classeSnwSc18M,dem_groundPointsVegSc_sp)
##%% vegtation classification from DEM2014 and las2014
#vegClassSc = upGroundPointsVegSc[:]
##all tree classification
#allTreeClassSc, treeNumClassSc = defineSpecificClassGreater (vegClassSc, 2)
#        
#negVegClassSc, negVegNumClassSc = defineSpecificClassLess (vegClassSc, 0)
#
## trees with low branches
#lowVegTreeClassSc, lowVegNumClassSc = defineLowVegClass(allTreeClassSc)
#
## trees with no low blanches
## "*************tall canopy no snow class*****************"
#nolowVegTreeClassSc, nolowVegTreeNumClassSc = differenceBetwee2classes (allTreeClassSc,lowVegNumClassSc)
#
## open space (no trees, no return between 0.15 to 2)
## all low veg
#allLowVegClassSc, allLowVegNumClassSc = defineLowVegClass(vegClassSc)
#
##open places
#notOpenNumClassSc = list(set(allLowVegNumClassSc).union(set(treeNumClassSc)))
## "******************open no snow class*******************"
#allOpenClassSc, allOpenNumClassSc = differenceBetwee2classes (vegClassSc,notOpenNumClassSc)
#
##%% snow classification for sagehen 26 March
#vegSnowClassSc26M, vegSnowNumClassSc26M = defineSpecificClassGreater (upGroundPointsSnwVegSc26M, -1)
#
##no snow on the ground
#nosnowClassSc26M, noSnowNumClassSc26M = defineSpecificClassLess (vegSnowClassSc26M, 0.15)
##snow on the ground or on the trees
#allSnowClassSc26M, allSnowNumClassSc26M = differenceBetwee2classes (vegSnowClassSc26M,noSnowNumClassSc26M)
## "******************tall canopy snow class*******************"
##snow on the tall canopy >2m
#treeSnowClassSc26M, treeSnowNumClassSc26M = defineSpecificClassGreater (allSnowClassSc26M, 2)
## "******************open snow class**************************"
##snow on the ground
#groundSnowClassSc26M, groundSnowNumClassSc26M = differenceBetwee2classes (allSnowClassSc26M,treeSnowNumClassSc26M)
#
##%% snow classification for sagehen 18 May
#vegSnowClassSc18M, vegSnowNumClassSc18M = defineSpecificClassGreater (upGroundPointsSnwVegSc18M, -1)
#
##no snow on the ground
#nosnowClassSc18M, noSnowNumClassSc18M = defineSpecificClassLess (vegSnowClassSc18M, 0.15)
##snow on the ground or on the trees
#allSnowClassSc18M, allSnowNumClassSc18M = differenceBetwee2classes (vegSnowClassSc18M,noSnowNumClassSc18M)
## "******************tall canopy snow class*******************"
##snow on the tall canopy >2m
#treeSnowClassSc18M, treeSnowNumClassSc18M = defineSpecificClassGreater (allSnowClassSc18M, 2)
## "******************open snow class**************************"
##snow on the ground
#groundSnowClassSc18M, groundSnowNumClassSc18M = differenceBetwee2classes (allSnowClassSc18M,treeSnowNumClassSc18M)
#
##%% ploting
#fig3 = plt.figure(figsize=(20,15))
#ax3 = Axes3D(fig3)
#ax3.scatter(dem_groundPointsVegSc_sp[:, 0], dem_groundPointsVegSc_sp[:, 1], dem_groundPointsVegSc_sp[:, 2])
##ax3.scatter(coordsVegSc_sp[:, 0], coordsVegSc_sp[:, 1], coordsVegSc_sp[:, 2])
##ax3.scatter(coordsSnwSc26M[:, 0], coordsSnwSc26M[:, 1], coordsSnwSc26M[:, 2])
#ax3.scatter(coordsSnwSc18M[:, 0], coordsSnwSc18M[:, 1], coordsSnwSc18M[:, 2])
#
#ax3.legend()
#plt.title('snow Lidar data sagehen May 18',fontsize=30)
#
##for flcl in failclass2:
##    ax3.scatter([x[0] for x in classes_rplc2[flcl]], [x[1] for x in classes_rplc2[flcl]])#, [x[2] for x in classesnow[flcl]])
#plt.savefig('H:\Chapter2_snow_forest\LIDAR_analysis\Sagehen\dem_groundPoints_sc_Veg&dem.png')
#
##%%
#
#def choosingFirstArray0fEachClassAndSettingValueForThatClass(classList,value):
#    class_1stArr = []
#    for arr in classList:
#        class_1stArr.append(arr[0][0])
#    class_1stArr_df = pd.DataFrame(class_1stArr, columns = ['x','y','z'])
#    class_1stArr_df['z'] = value
#    return class_1stArr_df
#
#Nosnow0nGroundNr1_df = choosingFirstArray0fEachClassAndSettingValueForThatClass(nosnowClassVcm,0)
#Nosnow0nGroundNr1_df_int = Nosnow0nGroundNr1_df.astype(int)
#
#Snow0nTallTreeNr1_df = choosingFirstArray0fEachClassAndSettingValueForThatClass(treeSnowClassVcm,0.2)
#Snow0nTallTreeNr1_df_int = Snow0nTallTreeNr1_df.astype(int)
#
#Snow0nGroundNr1_df = choosingFirstArray0fEachClassAndSettingValueForThatClass(groundSnowClassVcm,0.4)
#Snow0nGroundNr1_df_int = Snow0nGroundNr1_df.astype(int)
#
#figNr1 = plt.figure(figsize=(20,15))
#axNr1 = Axes3D(figNr1)
#axNr1.scatter(Nosnow0nGroundNr1_df_int['x'], Nosnow0nGroundNr1_df_int['y'], Nosnow0nGroundNr1_df_int['z'])
#axNr1.scatter(Snow0nGroundNr1_df_int['x'], Snow0nGroundNr1_df_int['y'], Snow0nGroundNr1_df_int['z'])
#axNr1.scatter(Snow0nTallTreeNr1_df_int['x'], Snow0nTallTreeNr1_df_int['y'], Snow0nTallTreeNr1_df_int['z'])
#plt.legend()
#plt.title('snow on the ground vs. on the tall canopy in Jemez_VCM',fontsize=30)
#
#plt.savefig("H:/Chapter2_snow_forest/LIDAR_analysis/Sagehen/snowDepthSc26M.png")
#
##SnowDemNR1 = pd.concat([Nosnow0nGroundNr1_df_int, Snow0nGroundNr1_df_int, Snow0nTallTreeNr1_df_int], axis=0, ignore_index=True)
##SnowDemNR1.sort_values(by=['x'])
##SnowDemNR1_drop = SnowDemNR1.drop_duplicates(subset=['x','y','z'], keep = "first")
##SnowDemNR1_drop.set_index([np.arange(0,len(SnowDemNR1_drop))],inplace=True)
##
##DemNr1 = dem_groundPointsVegNr1_sp0.copy()
##DemNr1_int = DemNr1.astype(int)
##DemNr1_int.set_index([np.arange(0,len(DemNr1_int))],inplace=True)
#
##for demPnt in range(len(DemNr1_int)):
##    for snwPnt in range(len(SnowDemNR1_drop)):
##        if (DemNr1_int['x'][demPnt] == SnowDemNR1_drop['X'][snwPnt]) & (DemNr1_int['y'][demPnt] == SnowDemNR1_drop['Y'][snwPnt]):
##            DemNr1_int['z'][demPnt] = SnowDemNR1_drop['Z'][snwPnt]
##        else:DemNr1_int['z'][demPnt] = -1
#
##Check = DemNr1_int['z'].replace(SnowDemNR1_drop['z'])
#
##DemNr1_int.loc[SnowDemNR1_drop['x'].isin(DemNr1_int.x) & SnowDemNR1_drop['y'].isin(DemNr1_int.y),'z']=SnowDemNR1_drop['z']
#
##DemNr1_int.loc[SnowDemNR1_drop['x']==DemNr1_int.x & SnowDemNR1_drop['y']==DemNr1_int.y,'z']=SnowDemNR1_drop['z']
#
##SameIndexNr1 = DemNr1_int[~(DemNr1_int['x'].isin(SnowDemNR1_drop['x']))].index# | (df1['A'].isnull())
##a = SnowDemNR1_drop[SnowDemNR1_drop['x'].isin(DemNr1_int['x'])].index
##df3['D'][a]=3
#
#
#
#
#
#
#
#
#
#
#
#
#
#
