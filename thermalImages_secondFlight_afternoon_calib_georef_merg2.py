import laspy as ls
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import csv
import os
import rasterio
from matplotlib import cm
from matplotlib import colors
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cv2
import utm
from sklearn.preprocessing import Imputer

#% functions
def readThermalImages(thrTreeFilePath,lstMissNo):
    demset = gdal.Open(thrTreeFilePath)
    band = demset.GetRasterBand(1)
    lst = band.ReadAsArray()
    lst[lst == -9999.] = lstMissNo
    lst[lst > 100000] = lstMissNo
    lst[lst < -20] = lstMissNo
        
    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = lst.shape
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows
    extent=[x0, x1, y1, y0]
    
    return [lst,extent]

#def readThermalImages(thrTreeFilePath,lstMissNo):
#    demset = gdal.Open(thrTreeFilePath)
#    band = demset.GetRasterBand(1)
#    lst = band.ReadAsArray()
#    lst[lst == -9999.] = lstMissNo
#    lst[lst > 100000] = lstMissNo
#    #lst[lst < -20] = lstMissNo
#        
#    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
#    nrows, ncols = lst.shape
#    x1 = x0 + dx * ncols
#    y1 = y0 + dy * nrows
#    extent=[x0, x1, y1, y0]
#    
#    return [lst,extent]
    
def readLasFile (lasFilePath):
    infile = ls.file.File(lasFilePath, mode="r")
    coords = np.vstack((infile.x, infile.y, infile.z)).T
    coords_df0 = pd.DataFrame(coords,columns=['x','y','z'])
    coords_df = pd.concat([coords_df0[['x','y']].astype(int),coords_df0['z']], axis=1)
    coords_df.sort_values(by=['x','y'],inplace=True)
    coords_df.index=np.arange(0,len(coords_df)) 
    return coords_df
    
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

#def plotThermalMap_lowResol(filename,lstMissNo,ouputPath,T): 
#    
#    demset = gdal.Open(filename)
#    band = demset.GetRasterBand(1)
#    lst = band.ReadAsArray()
#    lst[lst == -9999.] = lstMissNo
#    lst[lst > 100000] = lstMissNo
#    lst[lst < -20] = lstMissNo
#    
#    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
#    nrows, ncols = lst.shape
#    x1 = x0 + dx * ncols
#    y1 = y0 + dy * nrows
#    extent1=[x0, x1, y1, y0]
#    
#    plt.figure(figsize=(60,40))
#    #plt.imshow(lst, cmap='gist_earth', extent=extent1)
#    
#    #plotting
#    elseColor = cm.get_cmap('Greens_r', 128)
#    #bottom = cm.get_cmap('Blues', 128)
#    cmap = colors.ListedColormap(['white', elseColor])
#    labels = {lst<3:'snow',lst>3:'else'}
#    colors_cmap = cmap(np.arange(cmap.N))
#    cmap = {-99:colors_cmap[0], -7:colors_cmap[1], -4:colors_cmap[2], 44:colors_cmap[3], 77:colors_cmap[4]}
#
#    arrayShow = np.array([[cmap[i] for i in j] for j in lst])    
#    patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]
#
##    elev_2d = (np.reshape(snow_temp_vegdens_allExtent[:,2],[ncols,nrows])).astype(int)
##    contours = plt.contour(elev_2d, 15, colors='black', linewidths= 4, extent=lat_long_nad83)#
##    plt.clabel(contours, inline=True, fontsize=30, fmt = '%1.0f')
#    
#    plt.imshow(arrayShow, extent=extent1)#'gist_earth'
#
#    plt.xticks(fontsize=60, ha = 'left')#
#    plt.yticks(fontsize=60)
#
#    plt.legend(handles=patches, loc=2, fontsize=40, borderaxespad=0.)
#    plt.title('Snow presence map under forest canopy and open areas in %s' %T, fontsize=70, y=1.015)
#
#    plt.savefig(ouputPath)

def tempLatLon_from_thermalGeorefTiff(georef_calibrated_tiff):
    demset = gdal.Open(georef_calibrated_tiff)
    band = demset.GetRasterBand(1)
    temp = band.ReadAsArray()

    x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
    nrows, ncols = temp.shape
    #x1 = x0 + dx * ncols
    #y1 = y0 + dy * nrows
    #extent1=[x0, x1, y1, y0]
    
    latitude =[]
    for x in range (ncols):
        latitude.append(dx*x+x0)
    longitude = []
    for y in range (nrows):
        longitude.append(y0-dy*y)
   
    latitude_rp = (np.tile(latitude, nrows))
    longitude_rp = (np.repeat(longitude, ncols))
    temp_rp = np.reshape(temp,(nrows*ncols)).T
    temp_lat_lon1 = np.vstack([latitude_rp,longitude_rp,temp_rp]).T

    temp_df = pd.DataFrame(temp_lat_lon1,columns=['x','y','temp'])
    temp_df.sort_values(by=['x','y'],inplace=True)
    temp_df.index=np.arange(0,len(temp_df))
    
    return temp_df

def georeferencingThermalImagesAndGeneratingTIFFfiles (geolocationFile_degree,path_georefImage,lst_images):
    
    with open(geolocationFile_degree) as geoloc_dgr:
        reader = csv.reader(geoloc_dgr)
        geoloc_dgr_col = [r for r in reader]
    geoloc_degree_col = []
    for csv_counter1 in range (len (geoloc_dgr_col)):
        for csv_counter2 in range (4):
            geoloc_degree_col.append(geoloc_dgr_col[csv_counter1][csv_counter2])
    geoloc_degree=np.reshape(geoloc_degree_col,(len (geoloc_dgr_col),4))
    
    utm_lat = []
    utm_lon = []
    for dgr in range (len(geoloc_degree)):
        nrth = geoloc_degree[dgr][1].astype(np.float)
        east = geoloc_degree[dgr][2].astype(np.float)
        utm_geoloc = utm.from_latlon(nrth,east)
        utm_lat.append(utm_geoloc[1])
        utm_lon.append(utm_geoloc[0])

    #create tiff files from georeferenced (degree), calibrated images
    PIXEL_SIZE = 0.1 #offNadir20/x_pixels  # size of the pixel...        
    x_pixels = np.size(lst_images[0], axis = 1)
    y_pixels = np.size(lst_images[0], axis = 0)

    #getting projection from an image
    path_dem_scw = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/dem_sagehen.tif"
    demset = gdal.Open(path_dem_scw)
    wkt_projection = demset.GetProjection()

    #create georeference images
    for tif in range (len(utm_lat)):#
        dst_filename = path_georefImage[tif]
    
        x_min = utm_lon[tif]-PIXEL_SIZE*x_pixels/2 
        y_max = utm_lat[tif]+PIXEL_SIZE*y_pixels/2   # x_min & y_max are like the "top left" corner.
    
        driver = gdal.GetDriverByName('GTiff')
    
        dataset = driver.Create(dst_filename,x_pixels,y_pixels,1,gdal.GDT_Float32, )
    
        dataset.SetGeoTransform((x_min,PIXEL_SIZE,0,y_max,0,-PIXEL_SIZE))  
    
        dataset.SetProjection(wkt_projection)
        dataset.GetRasterBand(1).WriteArray(lst_images[tif])
        dataset.FlushCache()  # Write to disk.
        
#%% loading thermal data
lstMissNo = -9999

path2images= "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/Sagehen/11June2019/Afternoon/Thermal/SfM/tiff"
# get the list of images
listOfImgs=os.listdir(path2images)

# load image using p2i
fullPath2imgs = []
for p2i in listOfImgs:
    fullPath2img=os.path.join(path2images,p2i)
    fullPath2imgs.append(fullPath2img)

lst_hillslope1st = []
for path in fullPath2imgs:
    lst_tree = readThermalImages (path,lstMissNo)[0]
    lst_hillslope1st.append(lst_tree)

#plt.figure(figsize=(30,20))
#plt.hist(lst_hillslope1st[12], facecolor='navy')
#plt.xticks(fontsize=30)#
#plt.yticks(fontsize=30)
#plt.savefig("C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/hist_afternoon-0013.png")

#%% georeferencing thermal images and creating georeference images
geolocationFile_degree = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/geoloc_2flight.txt"
counter = np.arange(1000,1000+len(lst_hillslope1st)) #lst_hillslope_isoCalib1st
lst_hillslope_2flightAfternoon_georef_path = ["C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/thermalImages_afternoon_georef/lst_hillslope_2flightAternoon{}.tif".format(i) for i in counter]
lst_hillslope_2flightAfternoon_georef_name0nly = ["lst_hillslope_2flightAternoon{}.tif".format(i) for i in counter]

georeferencingThermalImagesAndGeneratingTIFFfiles (geolocationFile_degree,lst_hillslope_2flightAfternoon_georef_path,lst_hillslope1st)

#%% read georeferenced images and calibration with isothermal snow
lst_2flight_afternoon = []
for pth in range(len(lst_hillslope_2flightAfternoon_georef_path)):
    lst_tree = readThermalImages(lst_hillslope_2flightAfternoon_georef_path[pth],lstMissNo)[0]
    lst_2flight_afternoon.append(lst_tree)

lst_2flight_afternoon1 = readThermalImages(lst_hillslope_2flightAfternoon_georef_path[174],lstMissNo)

temp0b4_min = np.zeros((len(lst_2flight_afternoon),1))
for tm in range (len(lst_2flight_afternoon)):
    temp0b7s = lst_2flight_afternoon[tm][(lst_2flight_afternoon[tm]<=4)&(lst_2flight_afternoon[tm]>0.39)]
    if np.size(temp0b7s) == 0:
        temp0b4_min[tm] = np.nan
    else: temp0b4_min[tm] = np.min(temp0b7s)

#temp0b7_df = pd.DataFrame(temp0b7)
#temp0b7_df[0][0:100].fillna( method ='bfill', inplace = True, limit = 10) ##ffill  
temp0b4_min[0]=temp0b4_min[53]
temp0b4_min[90]=temp0b4_min[53]
temp0b4_min[145]=temp0b4_min[122]

imputer = Imputer(strategy = "median")
lngth = 2+len(temp0b4_min)/10
imp = range (0,len(temp0b4_min),lngth)
transformed_values = np.zeros((len(temp0b4_min),1))
for counter in imp:  #0,len(temp0b7),31
    values = temp0b4_min[counter:counter+lngth]
    transformed_values0 = imputer.fit_transform(values)
    transformed_values[counter:counter+lngth] = transformed_values0

temp_calib_ice = transformed_values.copy()
lst_2flight_afternoon_calib = []
for imgs in range (len(temp_calib_ice)):
    image_calib = lst_2flight_afternoon[imgs]-temp_calib_ice[imgs,0]
    image_calib[image_calib<=0]=0
    lst_2flight_afternoon_calib.append(image_calib)
    
#    unmixed_snow_pixels1st.append(np.ma.masked_where(lst_hillslope_isoCalib1st[lst] > 1, lst_hillslope_isoCalib1st[lst]))
#    unmixed_canopy_pixels1st.append(np.ma.masked_where(lst_hillslope_isoCalib1st[lst] < 10, lst_hillslope_isoCalib1st[lst]))
#
#unmixed_snow_canopy = []
#for callst in range (len(lst_hillslope_isoCalib1st)):
#    unmixed_snow_canopy_pixels = lst_hillslope_isoCalib1st[callst].copy()
#    unmixed_snow_canopy_pixels[(unmixed_snow_canopy_pixels>1)&(unmixed_snow_canopy_pixels<10)]=np.nan
#    unmixed_snow_canopy.append(unmixed_snow_canopy_pixels)

#creating tiff files of calibrated imagess  
counter2 = np.arange(1001,1001+len(lst_2flight_afternoon_calib)) #lst_hillslope_isoCalib1st
lst_hillslope_2flightAfternoon_georef_calib_path = ['C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/thermalImages_afternoon_georef/calibrated/calibWithIce4/lst_hillslope_2flightAternoon_calibIce4cmin_{}.tif'.format(i) for i in counter2]
georeferencingThermalImagesAndGeneratingTIFFfiles (geolocationFile_degree,lst_hillslope_2flightAfternoon_georef_calib_path,lst_2flight_afternoon_calib)
                           
#%% merged images ????? 
#calibrated by ice; not 20_0ffnadir merged images by pix4d 4.5 c
path_tiff1flight_afternoon_uc1 = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/stick_2flight_Afternoon_calibIce4cmin_3d/3_dsm_ortho/2_mosaic/stick_2flight_Afternoon_calibIce4cmin_3d_transparent_mosaic_grayscale.tif"
tiff1flight_afternoon_uc1_df = tempLatLon_from_thermalGeorefTiff(path_tiff1flight_afternoon_uc1)
tiff1flight_afternoon_uc1 = readThermalImages(path_tiff1flight_afternoon_uc1,-10)
tiff1flight_afternoon_uc10 = tiff1flight_afternoon_uc1[0]#[:3550,:]
tiff1flight_afternoon_uc10[tiff1flight_afternoon_uc10==-0]=np.nan
#tiff1flight_afternoon_uc10[:,2500:]=np.nan
extent = [737406.24736, 737661.91897, 4368096.8227, 4368532.629310001] #tiff1flight_afternoon_uc1[1]

cmap_else = cm.get_cmap('nipy_spectral')
fig, ax = plt.subplots(1, 1, figsize=(30,60)) #fig.set_figheight(30) fig.set_figwidth(60)
plt.xticks(fontsize=30, ha = 'left')#
plt.yticks(fontsize=30)
im = plt.imshow(tiff1flight_afternoon_uc10, interpolation='none', cmap=cmap_else, extent = extent)#tiff1flight_afternoon_uc1[1], vmin=2
axDivider = make_axes_locatable(ax)
cax1 = axDivider.append_axes("right", size="5%", pad="1%")
plt.colorbar(im, cax=cax1)
plt.xticks(fontsize=40, ha = 'left')#
plt.yticks(fontsize=40)
tiff1flight_afternoon_uc1_path = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/stick_2flight_Afternoon_calibIce4cmin_3d/June_Afternoon_calibIce4cmin.png"
plt.savefig(tiff1flight_afternoon_uc1_path)
  
#%% validation with rad tree temp and river temp                             
radTree_st =pd.DataFrame([["2019-06-11 14:50:00",28.74,30.75,30.9,27.58,32.25],
                          ["2019-06-11 14:55:00",28.29,30.41,30.44,27.33,31.83],
                          ["2019-06-11 15:00:00",27.66,29.88,29.8,27.82,31.46],
                          ["2019-06-11 15:05:00",28.5,30.81,31.21,27.73,32.28],
                          ["2019-06-11 15:10:00",27.73,29.8,30.45,26.31,31.63],
                          ["2019-06-11 15:15:00",27.58,30.15,30.79,26.66,31.8],
                          ["2019-06-11 15:20:00",27.7,30.09,30.14,25.93,31.43],
                          ["2019-06-11 15:25:00",27.62,30.3,30.22,26.29,31.77],
                          ["2019-06-11 15:30:00",27.3,30.48,30.28,27.04,32.46],
                          ["2019-06-11 15:35:00",27.09,31.25,31.7,28.1,33.32],
                          ["2019-06-11 15:40:00",26.76,30.75,31.1,27.69,32.82],
                          ["2019-06-11 15:45:00",26.34,30.51,30.86,26.66,32.49],
                          ["2019-06-11 15:50:00",25.15,29.57,29.74,25.46,31.3]], columns = ["TIMESTAMP","TT_C_N_Avg","TT_C_W_Avg","TT_C_S_Avg","TT_C_E_Avg","TT_C_T_Avg"])

radTree_calib_img =pd.DataFrame([["245",24.39,0.62],
                                 ["246",24,1.03]], columns = ["#image","TT_C_T_Avg","iceTemp_besideradTree"])
    
river_temp =pd.DataFrame([["2019-06-11 14:30",12.0],
                          ["2019-06-11 14:45",12.2],
                          ["2019-06-11 15:00",12.3],
                          ["2019-06-11 15:15",12.4],
                          ["2019-06-11 15:30",12.5],
                          ["2019-06-11 15:45",12.5]])
#%%calibration with view angle 10 degree resolution 
offNadir20 = 80 * np.tan((10* np.pi)/180) * 2 #12.8
offNadir45 = 80 * np.tan((22.5* np.pi)/180) * 2
x_pixels_num = np.size(lst_2flight_afternoon_calib[0], axis = 1)
y_pixels_num = np.size(lst_2flight_afternoon_calib[0], axis = 0)

x_pixels = int(x_pixels_num * offNadir20 / offNadir45)
y_pixels = int(y_pixels_num * offNadir20 / offNadir45)

x_indx1 = x_pixels_num/2 - x_pixels/2
x_indx2 = x_indx1 + x_pixels
y_indx1 = y_pixels_num/2 - y_pixels/2
y_indx2 = y_indx1 + y_pixels

lst_hillslope_2stAfternoon_0ffnadir20 = []
for lst in range (len(lst_2flight_afternoon_calib)):
    lst_ofnad20 = lst_2flight_afternoon_calib[lst][y_indx1:y_indx2, x_indx1:x_indx2]
    lst_hillslope_2stAfternoon_0ffnadir20.append(lst_ofnad20)   
    
#%% create all calibrated/20offnadir images in tiff and png format
counter_0ff20n = np.arange(1,len(lst_hillslope_2stAfternoon_0ffnadir20)+1) #lst_hillslope_isoCalib1st
lst_2flightAfternoon_georef_calib0ffnadir_path = ['C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/thermalImages_afternoon_georef/calibrated/0ffnadir20/lst_2stAfternoon_0ffnadir20_{}.tif'.format(i) for i in counter_0ff20n]
georeferencingThermalImagesAndGeneratingTIFFfiles (geolocationFile_degree,lst_2flightAfternoon_georef_calib0ffnadir_path,lst_hillslope_2stAfternoon_0ffnadir20)

lst_hillslope_1st0ffnadir20_nameImage = ['C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/thermalImages_afternoon_georef/calibrated/0ffnadir20/lst_2stAfternoon_0ffnadir20_{}.png'.format(i) for i in counter_0ff20n]
cmap_else = cm.get_cmap('nipy_spectral')

for maps in range (len(lst_hillslope_2stAfternoon_0ffnadir20)):
    fig, ax = plt.subplots(1, 1)
    fig.set_figheight(70)
    fig.set_figwidth(70)
    im = plt.imshow(lst_hillslope_2stAfternoon_0ffnadir20[maps], interpolation='none', cmap=cmap_else)#, vmin=2
    axDivider = make_axes_locatable(ax)
    cax1 = axDivider.append_axes("right", size="5%", pad="1%")
    plt.colorbar(im, cax=cax1)
    
    plt.xticks(fontsize=50, ha = 'left')#
    plt.yticks(fontsize=50)
    
    plt.savefig(lst_hillslope_1st0ffnadir20_nameImage[maps])
    
#%% ploting 2 maps side by side 
    
path_tiff1flight_afternoon_jun = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/second_flight/stick_2flight_Afternoon_calibIce4cmin_3d/3_dsm_ortho/2_mosaic/stick_2flight_Afternoon_calibIce4cmin_3d_transparent_mosaic_grayscale.tif"
tiff1flight_afternoon_jun1_df = tempLatLon_from_thermalGeorefTiff(path_tiff1flight_afternoon_jun)
tiff1flight_afternoon_jun1 = readThermalImages(path_tiff1flight_afternoon_jun,-10)
tiff1flight_afternoon_jun10 = tiff1flight_afternoon_jun1[0]#[:3550,:]
tiff1flight_afternoon_jun10[tiff1flight_afternoon_jun10==-0]=np.nan
#tiff1flight_afternoon_jun10[:,2500:]=np.nan
extentJ = tiff1flight_afternoon_jun1[1] #[737406.24736, 737661.91897, 4368096.8227, 4368532.629310001] 

path_tiff1flight_afternoon_apr = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/first_flight/stick_1flight_Afternoon_45iso_invrs/3_dsm_ortho/2_mosaic/stick_1flight_Afternoon_45iso_invrs_transparent_mosaic_grayscale.tif"
tiff1flight_afternoon_apr1_df = tempLatLon_from_thermalGeorefTiff(path_tiff1flight_afternoon_apr)
tiff1flight_afternoon_apr1 = readThermalImages(path_tiff1flight_afternoon_apr,-10)
tiff1flight_afternoon_apr10 = tiff1flight_afternoon_apr1[0][350:4400,300:2800]
tiff1flight_afternoon_apr10[tiff1flight_afternoon_apr10==-0]=np.nan
#tiff1flight_afternoon_apr10[:,2500:]=np.nan
lst_afternoon = [tiff1flight_afternoon_apr10,tiff1flight_afternoon_jun10]

fig, axs = plt.subplots(1, 2, figsize=(120,60)) #fig.set_figheight(30) fig.set_figwidth(60)

images = []
#for i in range(Nr):
for j in range(2): #Nc
    images.append(axs[j].imshow(lst_afternoon[j], cmap=cmap_else, extent = extentJ))
    axs[j].label_outer()
    axs[j].tick_params(labelsize = 40)
# Find the min and max of all colors for use in setting the color scale.
vmin = min(image.get_array().min() for image in images)
vmax = max(image.get_array().max() for image in images)
norm = colors.Normalize(vmin=vmin, vmax=vmax)
for im in images:
    im.set_norm(norm)

cbar = fig.colorbar(images[0], ax=axs)#, orientation='horizontal', fraction=.1)
cbar.ax.tick_params(labelsize=50) 

tiff1flight_afternoon_junapr_path = "C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/June_Afternoon_calibIce_aprjun.png"
plt.savefig(tiff1flight_afternoon_junapr_path)

#%%
#indxManip = []
#for tif in range (len(utm_lat)):
#    thrmImg = temp_lat_lon [tif]
#    indx = (thrmImg['x']-737000)+(thrmImg['y']-4367500)
#    indxManip.append(indx)
#
#test = []
#for smpl in range (len(indxManip[0])):#
#    indxMLWR = indxManip[1]-indxManip[0][smpl]
#    test.append(indxMLWR)
#lasFile1 = readLasFile (lasFilePath)   

#path_tree1 = lst_hillslope_1st0ffnadir20_nameImage_tif_georef[65]
#temp_df_65 = tempLatLon_from_thermalGeorefTiff(path_tree1)
#temp_df_65_x = (temp_df_65['x']*10).astype(int); temp_df_65_y = (temp_df_65['y']*10).astype(int); 
#
#path_tree2 = lst_hillslope_1st0ffnadir20_nameImage_tif_georef[66]
#temp_df_66 = tempLatLon_from_thermalGeorefTiff(path_tree2)
#temp_df_66_x = (temp_df_66['x']*10).astype(int); temp_df_66_y = (temp_df_66['y']*10).astype(int); 
##temp_df_66_x.values
#for ii in range (len(temp_df_66_x)):
#    aaa_test = [temp_df_66_x.values[1] == temp_df_65_x.values]
##lst_tree_conc = np.concatenate((lst_tree1,lst_tree2))

#%% 
#temp_threshold = {0: (min_max_temp[53][0]+min_max_temp[52][0])/2,
#                  58: (min_max_temp[64][0]+min_max_temp[65][0]+min_max_temp[66][0]+min_max_temp[67][0])/4,
#                  92:(min_max_temp[118][0]+min_max_temp[119][0])/2,
#                  #'180':(((min_max_temp[118][0]+min_max_temp[119][0])/(-2))-((min_max_temp[118][0]+min_max_temp[119][0])/(-2))/21),
#                  159:0, #(min_max_temp[184][0]),
#                  187:0, #(min_max_temp[184][0]+min_max_temp[199][0])/2,
#                  195:0, #(min_max_temp[199][0]),
#                  223:0, #(min_max_temp[228][0]),
#                  239:0, #(min_max_temp[248][0]),
#                  279:0} #(min_max_temp[298][0])

#iceTemp_uc_img =pd.DataFrame([["1031",6.9,7.1,np.nan],
#                          ["1041",7.05,7.6,8.7],
#                          ["1042",7.2,5.9,8.2],
#                          ["1043*",6.09,6.6,7.3],
#                          ["1048",6.7,6.03,6.9],
#                          ["1049",5,5.5],
#                          ["1050",6.2,6.4],
#                          ["1051",6.5,5.7],
#                          ["1052*",5.2,4.2,5.7],
#                          ["1053*",4.5,5.5,3.6],
#                          ["1056",5.7,6.05],
#                          ["1057",6.5,5.4],
#                          ["1064*",6.7,6.8,6.9,5.2,6.1,6.8],
#                          ["1065",6.1,6.2,5.1],
#                          ["1066*",5.3,4.8,6.6],
#                          ["1067",6.9,8.7,6.7],
#                          ["115",6.8,6.1],
#                          ["1116*",5.3,5.4],
#                          ["1117*",5.2,5.3],
#                          ["1118*",4.6,4.3],
#                          ["1119*",5.3,4.1,4.9],
#                          ["1120",5.5,4.2,4.8],
#                          ["1121",5,4.1],
#                          ["1122*",3.8,3.3],
#                          ["1123",4.1,5.3],
#                          ["1181",2.3,2.7],
#                          ["1182",2.3,2,2.2],
#                          ["1183",2.4,2.6,2.8],
#                          ["1184",2.2,3.1,3.8],
#                          ["1185",1.8,2.4,2.5],
#                          ["1242*",1.2,1.5,1.8],
#                          ["1243*",2.7,2.1,2.3],
#                          ["1244*",0.5,0.9,1.6,1.2],
#                          ["1245*",1.2,1.6,0.62,101],
#                          ["1246*",1.2,1.8],
#                          ["1262",5.9],
#                          ["1263",5.3,5.6],columns = ["#image","iceTemp"])
###plt.figure(figsize=(20,15))
##plt.hist(lst_radTree, bins = 1000, facecolor='navy')
##plt.xticks(fontsize=30)#
##plt.yticks(fontsize=30)


## calibration cheshmi
#ice_calib ={34:np.min(lst_2flight_afternoon[31]),
#            58:(np.min(lst_2flight_afternoon[49])+np.min(lst_2flight_afternoon[51])+np.min(lst_2flight_afternoon[52]))/3,
#            100:np.mean([6.7,6.8,6.9,5.2,6.1,6.8,6.1,6.2,5.1,5.3,4.8,6.6]), #1064, 1065, 1066
#            159:np.mean([6.8,6.1,5.3,5.4,5.2,5.3,4.6,4.3,5.3,4.1,4.9,5.5,4.2,4.8,5,4.1,3.8,3.3,4.1,5.3]), #1115-1123
#            198:np.mean([2.3,2.7,2.3,2,2.2,2.4,2.6,2.8,2.2,3.1,3.8,1.8,2.4,2.5]), #1181-1185
#            250:np.mean([1.2,1.5,1.8,2.7,2.1,2.3,0.5,0.9,1.6,1.2,1.2,1.6,0.62,1.1,1.2,1.8]), #1242-1245
#            304:np.mean([5.9,5.3,5.6])} #1262-1263
##np.argpartition(lst_2flight_afternoon1[0], 10)
##print(np.sort(lst_2flight_afternoon1[0].flatten())[:10]) # gives 10 min in an array
#
#temp_calib = np.zeros((len(lst_2flight_afternoon1),1))
#for td in range (len(ice_calib)):
#    key = np.sort(ice_calib.keys())[td]
#    temp_calib0 = ice_calib[key]
#    temp_calib[key:,0]=temp_calib0
    
    
    
#lst_tree = lst_hillslope_1st0ffnadir20[66]
#
#lst_radTree = lst_tree[77:140,87:145]
##unmixed_lst_radTree1st = np.ma.masked_where(lst_radTree < 5, lst_radTree)
#
#canopyTemp_min = 10 #7.5 unmixed canopy temperature
#lst_radTree_image = lst_radTree.copy()
#unmixed_radTemp = lst_radTree_image[(lst_radTree_image>canopyTemp_min)]
##unmixed_canopyTemp = np.ma.masked_where(lst_radTree_image < canopyTemp_min, lst_radTree_image)
#
#bin_radTemp = np.histogram(unmixed_radTemp, bins= 50)
#indx_max = np.argmax(bin_radTemp[0])  
#lst_hot1 = bin_radTemp[1][indx_max]
#lst_hot2 = bin_radTemp[1][indx_max+1]
#lst_radTree_medianHot = np.mean(lst_radTree_image[(lst_radTree_image>lst_hot1)&(lst_radTree_image<lst_hot2)])
#lst_radTree_median = np.mean(unmixed_radTemp)
#
#lst_cal1 = lst_radTree_image[lst]-mostFreqHot

#%% make a video
#image_folder = 'C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/isothermal_calib'
#video_name = 'C:/1UNRuniversityFolder/Dissertation/Chapter3&4/CTEM_flights/analysis/isothermal_calib/lst_afternoon1st.mp4'
#
#images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
#frame = cv2.imread(os.path.join(image_folder, images[0]))
#height, width, layers = frame.shape
#
##sort images
#import re
#convert = lambda text: int(text) if text.isdigit() else text
#alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
#images_sort = sorted(images, key = alphanum_key)
#
#video = cv2.VideoWriter(video_name, 0, 1, (width,height))
#
#for image in images_sort:
#    video.write(cv2.imread(os.path.join(image_folder, image)))
#
##cv2.destroyAllWindows()
#video.release()


