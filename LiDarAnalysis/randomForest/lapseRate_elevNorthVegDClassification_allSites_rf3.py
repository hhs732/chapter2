import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import scipy
import pandas as pd
import itertools
import csv
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import os
import rasterio
from scipy import stats
from matplotlib import colors
from matplotlib import cm
import matplotlib.patches as mpatches
import seaborn as sns
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
import merf
from lapseRate_elevNorthVegDClassification_allSites_rf_utility import *

#%% 
### check projection ##########################################################################
#demset = gdal.Open("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/20091201.tif")
#radiation_proj = demset.GetProjection()
#swr1 = readTiff_creatDF("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/20091201.tif")
#
#demset2 = gdal.Open("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/demSnw0ffKrew.tif")
#dem_proj = demset2.GetProjection()
#dem = readTiff_creatDF("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/demSnw0ffKrew.tif")
#
#demset.SetProjection(dem_proj)
#demset.GetProjection()
#
#band = demset.GetRasterBand(1)
#lwr = band.ReadAsArray()
#    
#x0, dx, dxdy, y0, dydx, dy = demset.GetGeoTransform()
#nrows, ncols = lwr.shape

#'PROJCS["NAD83 / UTM zone 11N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AUTHORITY["EPSG","26911"]]'
#'PROJCS["NAD83 / UTM zone 10N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AUTHORITY["EPSG","26910"]]'
################################################################################################
#load data KREW 2010
#ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/ascii_index_pr_krew.npy')
#
#veg_density_kr = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/veg_density_kr.npy')
#
#elev_krew = ascii_grid_veg_indexKrew[:,2]
#tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339 #R² = 0.5074
#
#veg_snow_temp_density_krew = np.vstack([ascii_grid_veg_indexKrew[:,0],ascii_grid_veg_indexKrew[:,1].astype(int),ascii_grid_veg_indexKrew[:,2],
#                                        ascii_grid_veg_indexKrew[:,3],ascii_grid_veg_indexKrew[:,4],ascii_grid_veg_indexKrew[:,5],
#                                        ascii_grid_veg_indexKrew[:,6],tempDJF_gm_krew,veg_density_kr[:,2],tempDJF_gm_krew,tempDJF_gm_krew]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_krew',veg_snow_temp_density_krew)

## slope less than 30
#snow_temp_vegdens_index_sL30_krew = veg_snow_temp_density_krew[(veg_snow_temp_density_krew[:,5]>-50)&(veg_snow_temp_density_krew[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_krew',snow_temp_vegdens_index_sL30_krew)

#loading LWR and creating 100b100 tif files
#path2images_LWR= "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/LRad_in_Ground_day/res1b1/"
counter = np.arange(1,122) #lst_hillslope_isoCalib1st
path2LWR_tiff_100 = ['C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/LRad_in_Ground_day/LWR_100_{}.tif'.format(i) for i in counter]
#latLon_lwr,path_lwr = changeResolution0fMultipleTiffsAndCreateDF(path2images_LWR,path2LWR_tiff_100)

#%## loading SWR and creating 100b100 tif files
#path2images_SWR = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/SRad_in_Ground_day/res1b1/"
path2SWR_tiff_100 = ['C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/SRad_in_Ground_day/SWR_100_{}.tif'.format(i) for i in counter]
#latLon_swr,path_swr = changeResolution0fMultipleTiffsAndCreateDF(path2images_SWR,path2SWR_tiff_100)

#test projection
#lwr_path_krew = "C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/LRad_in_Ground_day/20091201.tif"
#lwr_df_krew1 = readTiff_creatDF(lwr_path_krew)
#[min(lwr_df_krew1['x']),max(lwr_df_krew1['x']),min(lwr_df_krew1['y']),max(lwr_df_krew1['y'])]

#loading 100b100 lwr&swr tif files
latLon_lwr_krew = creatingMeanRadiationfrom100b100tifFiles(path2LWR_tiff_100)
lwr_ls_krew1_limlLatLong = pd.DataFrame(latLon_lwr_krew[(latLon_lwr_krew[:,1]>=4100400)&(latLon_lwr_krew[:,1]<=4105999)&
                                       (latLon_lwr_krew[:,0]>=303500)&(latLon_lwr_krew[:,0]<=306000)], columns = ['x','y','lwr'])
lwr_ls_krew1_limlLatLong2 = lwr_ls_krew1_limlLatLong.dropna()

latLon_swr_krew = creatingMeanRadiationfrom100b100tifFiles(path2SWR_tiff_100)
swr_ls_krew1_limlLatLong = pd.DataFrame(latLon_swr_krew[(latLon_swr_krew[:,1]>=4100400)&(latLon_swr_krew[:,1]<=4105999)&
                                       (latLon_swr_krew[:,0]>=303500)&(latLon_swr_krew[:,0]<=306000)], columns = ['x','y','swr'])
swr_ls_krew1_limlLatLong2 = swr_ls_krew1_limlLatLong.dropna()

# loading other feature from fusion
snow_temp_vegdens_index_sL30_krew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_krew.npy')

#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(snow_temp_vegdens_index_sL30_krew[:,0]),np.max(snow_temp_vegdens_index_sL30_krew[:,0]),np.min(snow_temp_vegdens_index_sL30_krew[:,1]),np.max(snow_temp_vegdens_index_sL30_krew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

elev_krew = snow_temp_vegdens_index_sL30_krew[:,2]
tempDJF_gm_krew = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339 #R² = 0.5074
elev_class_inx = [101,102,103,104,105,106,107,108,109,110]

snow_temp_vegdens_index_sL30_krew_df = pd.DataFrame(snow_temp_vegdens_index_sL30_krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_krew = np.arange(min(snow_temp_vegdens_index_sL30_krew_df['z']),max(snow_temp_vegdens_index_sL30_krew_df['z']),67)
snow_temp_vegdens_sL30_elevCls_krew = classificationElevation(snow_temp_vegdens_index_sL30_krew,elevationBand_krew,elev_class_inx)

meanTemp_elevClass_krew = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_krew, fsca_0p_krew = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)
#fsca_ut_samp_krew, fsca_0p_samp_krew = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx,sample_size)
##statistical test
#statistic_krew1,pvalue_krew1,average_sample_exp_krew1,average_sample_shl_krew1,significant_krew1 = wilcoxonTest(fsca_ut_samp_krew, fsca_0p_samp_krew)

# under canopy and 0pen fsca----exposed vs. sheltered
fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_Krew, fSCA_0u_shl_Krew = fsca0pn_fscaUnTr(fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew, elev_class_inx)
#statistical test
#sample_size2 = 100      
#fSCA_0u_exp_Krew_samp, fSCA_0u_shl_Krew_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx,sample_size2)
#statistic_Krew2,pvalue_Krew2,average_sample_exp_Krew2,average_sample_shl_Krew2,significant_krew2 = wilcoxonTest(fSCA_0u_exp_Krew_samp, fSCA_0u_shl_Krew_samp)
#
##y = -0.003300x + 7.287576
tempDec_krew = -0.0033*elevationBand_krew + 7.287576 

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_krew = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx)
#fSCA_ut_exp_l_sam_krew, fSCA_ut_exp_h_sam_krew, fSCA_ut_shl_l_sam_krew, fSCA_ut_shl_h_sam_krew = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_krew,elev_class_inx,sample_size)
#statistic_ex_krew,pvalue_ex_krew,average_sample_l_ex_krew,average_sample_h_ex_krew,significant_ex_krew = wilcoxonTest(fSCA_ut_exp_l_sam_krew, fSCA_ut_exp_h_sam_krew)
#statistic_sh_krew,pvalue_sh_krew,average_sample_l_sh_krew,average_sample_h_sh_krew,significant_sh_krew = wilcoxonTest(fSCA_ut_shl_l_sam_krew, fSCA_ut_shl_h_sam_krew)
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/significant_krew',[significant_krew1,significant_krew2,significant_ex_krew,significant_sh_krew])

### precipitation dec to march
elev_krew = snow_temp_vegdens_sL30_elevCls_krew[:,2]
precip_krew1b1 = np.reshape((0.0596 * elev_krew - 13.978),(len(snow_temp_vegdens_sL30_elevCls_krew),1))
#adding precip column
snow_temp_vegdens_sL30_elevCls_krew = snow_temp_vegdens_sL30_elevCls_krew[:,:10]

snow_temp_vegdens_sL30_elevCls_precip_krew = np.append(snow_temp_vegdens_sL30_elevCls_krew, precip_krew1b1, 1)
aaaa_test1 = snow_temp_vegdens_sL30_elevCls_precip_krew[0:105710,:]

## snow maping ##############################################################################
#snow_lwr_vegdens_allExtent_krew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_krew.npy')
#ouputPath_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_contour_krew_lr3.png'
#T_krew = 'KREW'
#plotSnowMap_lowResol(snow_temp_vegdens_allExtent_krew,ouputPath_krew,T_krew)
#
##lat_long_nad83_krew = [np.min(snow_temp_vegdens_allExtent_krew[:,0]),np.max(snow_temp_vegdens_allExtent_krew[:,0]),np.min(snow_temp_vegdens_allExtent_krew[:,1]),np.max(snow_temp_vegdens_allExtent_krew[:,1])]
#
#ouputPath_nrt_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_nrth_contour_krew_lr3.png'
#plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent_krew,ouputPath_nrt_krew,T_krew)

# change resolution to 100b100 #######################################################################
snow_temp_vegdens_sL30_elevCls_precip_cls_krew = addClassifier2changeResolutionTo100B100 (lat_long_nad83_krew,snow_temp_vegdens_sL30_elevCls_precip_krew)
classifier_krew = snow_temp_vegdens_sL30_elevCls_precip_cls_krew [:,13]

in_out_rf_krew0 = define_rf_features_fscaLatLonElevNorthTempPrecip_100B100 (classifier_krew,snow_temp_vegdens_sL30_elevCls_precip_cls_krew)

#define random forest AI features ###################################################################
# defining index for adding lwr and swr
in_out_rf_krew = in_out_rf_krew0.copy()
roundingCoordinates (in_out_rf_krew) #rounding coodinates
in_out_rf_krew_indx = 100*((in_out_rf_krew ['x']).astype(int))+(in_out_rf_krew ['y']).astype(int) #(in_out_rf_krew ['x']-min(in_out_rf_krew ['x']))+(in_out_rf_krew ['y']-min(in_out_rf_krew ['y']))
in_out_rf_krew2 = pd.concat([in_out_rf_krew,in_out_rf_krew_indx],axis =1)
in_out_rf_krew3 = pd.concat([in_out_rf_krew2,in_out_rf_krew_indx],axis =1)
in_out_rf_krew_lswr = pd.concat([in_out_rf_krew3,in_out_rf_krew_indx],axis =1)
in_out_rf_krew_lswr.columns = ['x','y','elev','nrth','temp','precip','VegD','vegD_ut','fsca_0p','fsca_ut','fsca','indx','lwr','swr']

#adding lwr and swr to the features
swr_ls_krew1_limlLatLong20 = swr_ls_krew1_limlLatLong2.copy()
swr_ls_krew1_limlLatLong20.index = np.arange(len(swr_ls_krew1_limlLatLong20))
roundingCoordinates (swr_ls_krew1_limlLatLong20) #rounding coodinates
swr_ls_krew1_limLatLong_indx = 100*((swr_ls_krew1_limlLatLong20 ['x']).astype(int))+(swr_ls_krew1_limlLatLong20 ['y']).astype(int)
swr_ls_krew1_limLatLong3 = pd.concat([swr_ls_krew1_limlLatLong20,swr_ls_krew1_limLatLong_indx],axis =1)

lwr_ls_krew1_limlLatLong20 = lwr_ls_krew1_limlLatLong2.copy()
lwr_ls_krew1_limlLatLong20.index = np.arange(len(lwr_ls_krew1_limlLatLong20))
roundingCoordinates (lwr_ls_krew1_limlLatLong20) #rounding coodinates
lwr_ls_krew1_limLatLong_indx = 100*((lwr_ls_krew1_limlLatLong20 ['x']).astype(int))+(lwr_ls_krew1_limlLatLong20 ['y']).astype(int)
lwr_ls_krew1_limLatLong3 = pd.concat([lwr_ls_krew1_limlLatLong20,lwr_ls_krew1_limLatLong_indx],axis =1)

in_out_rf_krew_lswr_ls = (in_out_rf_krew_lswr.values).copy()
lwr_ls_krew1_limLatLong_ls = lwr_ls_krew1_limLatLong3.values
swr_ls_krew1_limLatLong_ls = swr_ls_krew1_limLatLong3.values
for smpl in range (len(lwr_ls_krew1_limLatLong_indx)):#
    indxMLWR = abs(in_out_rf_krew_lswr_ls[:,11]-lwr_ls_krew1_limLatLong_ls[smpl,3])
    indxMSWR = abs(in_out_rf_krew_lswr_ls[:,11]-swr_ls_krew1_limLatLong_ls[smpl,3])
    in_out_rf_krew_lswr_ls[:,12][indxMLWR<1]=lwr_ls_krew1_limLatLong_ls[smpl,2]
    in_out_rf_krew_lswr_ls[:,13][indxMSWR<1]=swr_ls_krew1_limLatLong_ls[smpl,2]
    
in_out_rf_krew_lswr_ls_drp = in_out_rf_krew_lswr_ls[in_out_rf_krew_lswr_ls[:,11] != in_out_rf_krew_lswr_ls[:,12]]    
in_out_rf_krew_lswr_ls_drp_df = pd.DataFrame(in_out_rf_krew_lswr_ls_drp, columns = ['x','y','elev','nrth','temp','precip','VegD','vegD_ut','fsca_0p','fsca_ut','fsca','indx','lwr','swr'])

# add classifier for fsca
in_out_rf_krew_lswr_ls_drp_df2 = pd.concat([in_out_rf_krew_lswr_ls_drp_df,in_out_rf_krew_lswr_ls_drp_df['fsca']],axis =1)
in_out_rf_krew_lswr_ls_drp_ls = in_out_rf_krew_lswr_ls_drp_df2.values
in_out_rf_krew_lswr_ls_drp_ls[:,14][in_out_rf_krew_lswr_ls_drp_ls[:,10]<=0.3] = 1
in_out_rf_krew_lswr_ls_drp_ls[:,14][(in_out_rf_krew_lswr_ls_drp_ls[:,10]<=0.55) & (in_out_rf_krew_lswr_ls_drp_ls[:,10]>=0.3)] = 2
in_out_rf_krew_lswr_ls_drp_ls[:,14][(in_out_rf_krew_lswr_ls_drp_ls[:,10]<=0.8) & (in_out_rf_krew_lswr_ls_drp_ls[:,10]>=0.55)] = 3
in_out_rf_krew_lswr_ls_drp_ls[:,14][in_out_rf_krew_lswr_ls_drp_ls[:,10]>0.8] = 4

in_out_rf_krew_lswr_fcsa_drp_df = pd.DataFrame(in_out_rf_krew_lswr_ls_drp_ls, columns = ['x','y','elev','nrth','temp','precip','VegD','vegD_ut','fsca_0p','fsca_ut','fsca','indx','lwr','swr','fsca_classifier'])

# random forest modeling #################################################################################
## random forest modeling for different fsca bins for fsca_open--fsca_ut

feature_importances_error_0putBin = randomForestModel44FscaBins(in_out_rf_krew_lswr_fcsa_drp_df)

plt.subplots(figsize=(20,15))
for countFeat in range (len(feature_importances_error_0putBin)): 
    plt.plot(feature_importances_error_0putBin[countFeat][0].values[0], marker='s', markersize=20, linewidth=3.5)

meanAbsError = [feature_importances_error_0putBin[i][1] for i in range (4)]
legend = ['bin{}; MAE = {}'.format(i+1,meanAbsError[i]) for i in range (4)]
plt.legend(legend,fontsize = 40)#, loc = 'upper center'
plt.title('Importance of features in 4 fscaBins in Krew', fontsize=40)#, position=(0.04, 0.88), ha='left'

ax = np.reshape(np.arange(0,len(feature_importances_error_0putBin[0][0].columns)),(1,5))
sa_xticks = feature_importances_error_0putBin[0][0].columns
plt.xticks(ax[0], sa_xticks, fontsize=40)#, rotation=25 
plt.yticks(fontsize=40)
plt.xlabel('Features', fontsize=40)
plt.ylabel('Importance', fontsize=40)

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/importance_0put_krew_bin.png')


#features and target for random forest model for fsca_open ###########################################
in_0p_rf_krew = in_out_rf_krew_lswr_ls_drp_df[['temp','precip','lwr','swr']] #''nrth',elev','x','y', 0.065 random_state=50
out_0p_rf_krew = in_out_rf_krew_lswr_ls_drp_df['fsca_0p']
[feature_importances_0p_mean,mean_abs_error_0p] = randomForestModel4Fsca0pen (in_0p_rf_krew,out_0p_rf_krew,50,200)

pathName_0p = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/importance_0p_krew_lswr_non_150502 (error='
plot_rfm_importance (mean_abs_error_0p, feature_importances_0p_mean, pathName_0p, 'red', 'fSCA_open')

#### 2d plotting fsca_open ###############################################################################
# 2d plotting by color map temp & nrth
savepath_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/temp_nrth_fSCA0pen3.png'
fscaPlottingByColorMap4tempNrth(in_out_rf_krew_lswr_ls_drp_df,out_0p_rf_krew,'OrRd',"fSCA_open",savepath_krew)
   
# 2d plotting by color map LWR & SWR
savepath_lswr_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/swr_lwr_fSCA0pen3.png'
fscaPlottingByColorMap4swrLwr(in_out_rf_krew_lswr_ls_drp_df,out_0p_rf_krew,'OrRd',"fSCA_open",savepath_lswr_krew)
    
#####random forest AI modeling for fsca_underTree  ####################################################
in_ut_rf_krew = in_out_rf_krew_lswr_ls_drp_df[['temp','precip','VegD','lwr','swr']] #vegD_ut,'nrth','elev','x','y', 0.065 random_state=50
out_ut_rf_krew = in_out_rf_krew_lswr_ls_drp_df['fsca_ut']

[feature_importances_ut_mean,mean_abs_error_ut] = randomForestModel4FscaUnderTree0r40p_ut(in_ut_rf_krew,out_ut_rf_krew,50,200)
pathName_ut = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/importance_ut_krew_non_lswr2 (error='
plot_rfm_importance (mean_abs_error_ut, feature_importances_ut_mean, pathName_ut, 'darkorange', 'fSCA_underTree')

#### 2d plotting fsca_ut ###############################################################################
# 2d plotting by color map temp & nrth
savepath_tempNrth_ut_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/temp_nrth_fSCAunderTree2.png'
fscaPlottingByColorMap4tempNrth(in_out_rf_krew_lswr_ls_drp_df,out_ut_rf_krew,'BuGn',"fSCA_underTree",savepath_tempNrth_ut_krew)
# 2d plotting by color map LWR & SWR
savepath_slwr_ut_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/swr_lwr_fSCAunderTree2.png'
fscaPlottingByColorMap4swrLwr(in_out_rf_krew_lswr_ls_drp_df,out_ut_rf_krew,'BuGn',"fSCA_underTree",savepath_slwr_ut_krew)
    
#####random forest AI modeling for fSCA_0pen -- fsca_underTree  ##########################################
in_ut_0p_rf_krew = in_out_rf_krew_lswr_ls_drp_df[['temp','precip','VegD','lwr','swr']] #'nrth','elev','x','y', 0.065 random_state=50
out_ut_0p_rf_krew = pd.DataFrame(in_out_rf_krew_lswr_ls_drp_df['fsca_ut']-in_out_rf_krew_lswr_ls_drp_df['fsca_0p'], columns = ['fSCA_op_M_fSCA_ut']).iloc[:,0]
[feature_importances_0p_ut_mean2,mean_abs_error_oput2] = randomForestModel4FscaUnderTree0r40p_ut(in_ut_0p_rf_krew,out_ut_0p_rf_krew,59,200)

pathName_0p_ut = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/importance_0p_ut_krew_lswr_ne_non (error='
plot_rfm_importance (mean_abs_error_oput, feature_importances_0p_ut_mean, pathName_0p_ut, 'goldenrod', 'fSCAop-fSCAut')

#### 2d plotting fsca_open_ut ###############################################################################
# 2d plotting by color map temp & nrth
savepath_tempNrth_ut0p_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/temp_nrth_fSCA0p_fSCAut2.png'
fscaPlottingByColorMap4tempNrth(in_out_rf_krew_lswr_ls_drp_df,out_ut_0p_rf_krew,cm.Purples_r,"fSCA0p_fSCAut",savepath_tempNrth_ut0p_krew)
# 2d plotting by color map LWR & SWR
savepath_slwr_ut0p_krew = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/randomForest/Krew/swr_lwr_fSCA0p_fSCAut2.png'
fscaPlottingByColorMap4swrLwr(in_out_rf_krew_lswr_ls_drp_df,out_ut_0p_rf_krew,cm.Purples_r,"fSCA0p_fSCAut",savepath_slwr_ut0p_krew)

#%% load data nrc 2010

#ascii_grid_veg_indexnrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/ascii_index_pr_nrc.npy')
#veg_density_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/veg_density_bl.npy')

#lat_long_nad83_nrc = [np.min(snow_temp_vegdens_allExtent_nrc[:,0]),np.max(snow_temp_vegdens_allExtent_nrc[:,0]),np.min(snow_temp_vegdens_allExtent_nrc[:,1]),np.max(snow_temp_vegdens_allExtent_nrc[:,1])]
#lat_long_degre_nrc = [39.978183, 40.025319, -105.58556840332032, -105.52620245393197]

##y = -0.007512x + 14.168475  # R² = 0.969091
#elev_nrc = ascii_grid_veg_indexnrc[:,2]
#tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475 
#
#veg_snow_temp_density_nrc = np.vstack([ascii_grid_veg_indexnrc[:,0],ascii_grid_veg_indexnrc[:,1].astype(int),ascii_grid_veg_indexnrc[:,2],
#                                       ascii_grid_veg_indexnrc[:,3],ascii_grid_veg_indexnrc[:,4],ascii_grid_veg_indexnrc[:,5],
#                                       ascii_grid_veg_indexnrc[:,6],tempDJF_gm_nrc,veg_density_nrc[:,2],tempDJF_gm_nrc,tempDJF_gm_nrc]).T
# slope less than 30
#snow_temp_vegdens_index_sL30_nrc = veg_snow_temp_density_nrc[(veg_snow_temp_density_nrc[:,5]>-90)&(veg_snow_temp_density_nrc[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_nrc',veg_snow_temp_density_nrc)

snow_temp_vegdens_index_sL30_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_nrc.npy')
aaaa_test1 = snow_temp_vegdens_index_sL30_nrc[0:19710,:]

elev_nrc = snow_temp_vegdens_index_sL30_nrc[:,2]
tempDJF_gm_nrc = -0.007512*elev_nrc + 14.168475

snow_temp_vegdens_index_sL30_nrc_df = pd.DataFrame(snow_temp_vegdens_index_sL30_nrc, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elev_class_inx = [101,102,103,104,105,106,107,108,109,110]
elevationBand_nrc = np.arange(min(snow_temp_vegdens_index_sL30_nrc_df['z']),max(snow_temp_vegdens_index_sL30_nrc_df['z']),65)
snow_temp_vegdens_sL30_elevCls_nrc = classificationElevation(snow_temp_vegdens_index_sL30_nrc,elevationBand_nrc,elev_class_inx)
aaaa_test2 = snow_temp_vegdens_sL30_elevCls_nrc[0:710,:]

meanTemp_elevClass_nrc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)   

#sampling
sample_size = 100
fsca_ut_nrc, fsca_0p_nrc = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)
#fsca_ut_samp_nrc, fsca_0p_samp_nrc = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx,sample_size)
#statistic_nrc1,pvalue_nrc1,average_sample_exp_nrc1,average_sample_shl_nrc1,significant_nrc1 = wilcoxonTest(fsca_ut_samp_nrc, fsca_0p_samp_nrc)

# under canopy and 0pen fsca // exposed vs. sheltered
fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_nrc, fSCA_0u_shl_nrc = fsca0pn_fscaUnTr(fSCA_ut_exp_nrc, fSCA_ut_shl_nrc, fSCA_0p_exp_nrc, fSCA_0p_shl_nrc, elev_class_inx)
#sampling
#fSCA_0u_exp_nrc_samp, fSCA_0u_shl_nrc_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx,sample_size)
#statistic_nrc2,pvalue_nrc2,average_sample_exp_nrc2,average_sample_shl_nrc2,significant_nrc2 = wilcoxonTest(fSCA_0u_exp_nrc_samp, fSCA_0u_shl_nrc_samp)

## December
##y = -0.00777x + 17.81629
tempDec_nrc = -0.00777*elevationBand_nrc + 17.81629

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_nrc = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx)
#fSCA_ut_exp_l_sam_nrc, fSCA_ut_exp_h_sam_nrc, fSCA_ut_shl_l_sam_nrc, fSCA_ut_shl_h_sam_nrc = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_nrc,elev_class_inx,sample_size)
#statistic_ex_nrc,pvalue_ex_nrc,average_sample_l_ex_nrc,average_sample_h_ex_nrc,significant_ex_nrc = wilcoxonTest(fSCA_ut_exp_l_sam_nrc, fSCA_ut_exp_h_sam_nrc)
#statistic_sh_nrc,pvalue_sh_nrc,average_sample_l_sh_nrc,average_sample_h_sh_nrc,significant_sh_nrc = wilcoxonTest(fSCA_ut_shl_l_sam_nrc, fSCA_ut_shl_h_sam_nrc)
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/significant_nrc',[significant_nrc1,significant_nrc2,significant_ex_nrc,significant_sh_nrc])

##calculating topo_dimension from cut northness to plot snow on
#snow_temp_vegdens_allExtent_nrc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_nrc.npy')
#ouputPath_nrc = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_contour_nrc_lr2.png'
#T_nrc = 'NRC 2010'
#plotSnowMap_lowResol(snow_temp_vegdens_allExtent_nrc,ouputPath_nrc,T_nrc)
#
#ouputPath_nrt_nrc = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_nrth_contour_nrc_lr2.png'
#plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent_nrc,ouputPath_nrt_nrc,T_nrc)

#%% load data Jemez 2010
#ascii_grid_veg_indexJmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/ascii_index_pr_jmz.npy')
#veg_density_jmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/veg_density_jm.npy')
#
#elev_jmz = ascii_grid_veg_indexJmz[:,2]
#tempDJF_gm_jmz = -0.0046*elev_jmz + 7.7058 
#
#veg_snow_temp_density_jmz = np.vstack([ascii_grid_veg_indexJmz[:,0],ascii_grid_veg_indexJmz[:,1].astype(int),ascii_grid_veg_indexJmz[:,2],
#                                       ascii_grid_veg_indexJmz[:,3],ascii_grid_veg_indexJmz[:,4],ascii_grid_veg_indexJmz[:,5],
#                                       ascii_grid_veg_indexJmz[:,6],tempDJF_gm_jmz,veg_density_jmz[:,2],tempDJF_gm_jmz,tempDJF_gm_jmz]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_jmz',veg_snow_temp_density_jmz)

## slope less than 30
#snow_temp_vegdens_index_sL30_jmz = veg_snow_temp_density_jmz[(veg_snow_temp_density_jmz[:,5]>-50)&(veg_snow_temp_density_jmz[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_jmz',snow_temp_vegdens_index_sL30_jmz)

snow_temp_vegdens_index_sL30_jmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_jmz.npy')
aaaa_test1 = snow_temp_vegdens_index_sL30_jmz[0:710,:]

elev_jmz = snow_temp_vegdens_index_sL30_jmz[:,2]
tempDJF_gm_jmz = -0.0046*elev_jmz + 7.7058 

#calculating topo_dimension from cut northness
lat_long_nad83_Jmz = [np.min(snow_temp_vegdens_index_sL30_jmz[:,0]),np.max(snow_temp_vegdens_index_sL30_jmz[:,0]),np.min(snow_temp_vegdens_index_sL30_jmz[:,1]),np.max(snow_temp_vegdens_index_sL30_jmz[:,1])]
lat_long_degre_Jmz = [35.863504, 35.918451, -106.60600014573875, -106.54061446340079]

snow_temp_vegdens_index_sL30_jmz_df = pd.DataFrame(snow_temp_vegdens_index_sL30_jmz, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_jmz = np.arange(min(snow_temp_vegdens_index_sL30_jmz_df['z']),max(snow_temp_vegdens_index_sL30_jmz_df['z']),90)
snow_temp_vegdens_sL30_elevCls_jmz = classificationElevation(snow_temp_vegdens_index_sL30_jmz,elevationBand_jmz,elev_class_inx)

meanTemp_elevClass_jmz = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_jmz, fsca_0p_jmz = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)
#fsca_ut_samp_jmz, fsca_0p_samp_jmz = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx,sample_size)
##statistical test
#statistic_jmz1,pvalue_jmz1,average_sample_exp_jmz1,average_sample_shl_jmz1,significant_jmz1 = wilcoxonTest(fsca_ut_samp_jmz, fsca_0p_samp_jmz)

# under canopy and 0pen fsca exposed vs. open
fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz = fsca0pn_fscaUnTr(fSCA_ut_exp_Jmz, fSCA_ut_shl_Jmz, fSCA_0p_exp_Jmz, fSCA_0p_shl_Jmz, elev_class_inx)
##statistical test
#fSCA_0u_exp_Jmz_samp, fSCA_0u_shl_Jmz_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx,sample_size)
#statistic_jmz2,pvalue_jmz2,average_sample_exp_jmz2,average_sample_shl_jmz2,significant_jmz2 = wilcoxonTest(fSCA_0u_exp_Jmz_samp, fSCA_0u_shl_Jmz_samp)

# y = -0.0043x + 5.6707
tempDec_jmz = -0.0048*elevationBand_jmz + 8.6707

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_jmz = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx)
#fSCA_ut_exp_l_sam_jmz, fSCA_ut_exp_h_sam_jmz, fSCA_ut_shl_l_sam_jmz, fSCA_ut_shl_h_sam_jmz = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_jmz,elev_class_inx,sample_size)
#statistic_ex_jmz,pvalue_ex_jmz,average_sample_l_ex_jmz,average_sample_h_ex_jmz,significant_ex_jmz = wilcoxonTest(fSCA_ut_exp_l_sam_jmz, fSCA_ut_exp_h_sam_jmz)
#statistic_sh_jmz,pvalue_sh_jmz,average_sample_l_sh_jmz,average_sample_h_sh_jmz,significant_sh_jmz = wilcoxonTest(fSCA_ut_shl_l_sam_jmz, fSCA_ut_shl_h_sam_jmz)
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/significant_jmz',[significant_jmz1,significant_jmz2,significant_ex_jmz,significant_sh_jmz])
#
#calculating topo_dimension from cut northness to plot snow on
#snow_temp_vegdens_allExtent_jmz = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_jmz.npy')
#ouputPath_jmz = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_contour_jmz_lr2.png'
#T_jmz = 'JRBN'
#plotSnowMap_lowResol(snow_temp_vegdens_allExtent_jmz,ouputPath_jmz,T_jmz)
#
##lat_long_nad83_jmz = [np.min(snow_temp_vegdens_allExtent_jmz[:,0]),np.max(snow_temp_vegdens_allExtent_jmz[:,0]),np.min(snow_temp_vegdens_allExtent_jmz[:,1]),np.max(snow_temp_vegdens_allExtent_jmz[:,1])]
#
#ouputPath_nrt_jmz = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_nrth_contour_jmz_lr2.png'
#plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent_jmz,ouputPath_nrt_jmz,T_jmz)

#%% load data sagehen creek 26 March 2016
#ascii_grid_veg_index26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/ascii_index_pr_sc26mar_retile.npy')
##calculating topo_dimension from cut northness
#veg_density_sc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/veg_density_sc.npy')
#
#elev_sc26m = ascii_grid_veg_index26m[:,2]
#tempDJF_sc26m = -0.0016 * elev_sc26m + 1.802
#
#veg_snow_temp_density_sc26m = np.vstack([ascii_grid_veg_index26m[:,0],ascii_grid_veg_index26m[:,1].astype(int),ascii_grid_veg_index26m[:,2],
#                                         ascii_grid_veg_index26m[:,3],ascii_grid_veg_index26m[:,4],ascii_grid_veg_index26m[:,5],
#                                         ascii_grid_veg_index26m[:,6],tempDJF_sc26m,veg_density_sc[:,2],tempDJF_sc26m,tempDJF_sc26m]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sc26m_retile',veg_snow_temp_density_sc26m)
## slope less than 30
#snow_temp_vegdens_index_sL30_sc26m = veg_snow_temp_density_sc26m[(veg_snow_temp_density_sc26m[:,5]>-50)&(veg_snow_temp_density_sc26m[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_sc26m_retile',snow_temp_vegdens_index_sL30_sc26m)

snow_temp_vegdens_index_sL30_sc26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_sc26m_retile.npy')

elev_sc26m = snow_temp_vegdens_index_sL30_sc26m[:,2]
tempDJF_sc26m = -0.0016 * elev_sc26m + 1.802

snow_temp_vegdens_index_sL30_sc26m_df = pd.DataFrame(snow_temp_vegdens_index_sL30_sc26m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_sc = np.arange(min(snow_temp_vegdens_index_sL30_sc26m_df['z']),max(snow_temp_vegdens_index_sL30_sc26m_df['z']),67)
snow_temp_vegdens_sL30_elevCls_sc26m = classificationElevation(snow_temp_vegdens_index_sL30_sc26m,elevationBand_sc,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)   

aaaa_test1 = snow_temp_vegdens_sL30_elevCls_sc26m[0:710,:]

# under canopy and 0pen fsca
fsca_ut_sc26m, fsca_0p_sc26m = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)
#fsca_ut_samp_sc26m, fsca_0p_samp_sc26m = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx,sample_size)
##statistical test
#statistic_sc26m1,pvalue_sc26m1,average_sample_exp_sc26m1,average_sample_shl_sc26m1,significant_sc26m1 = wilcoxonTest(fsca_ut_samp_sc26m, fsca_0p_samp_sc26m)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc26m, fSCA_ut_shl_sc26m, fSCA_0p_exp_sc26m, fSCA_0p_shl_sc26m = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc26m, fSCA_0u_shl_sc26m = fsca0pn_fscaUnTr(fSCA_ut_exp_sc26m, fSCA_ut_shl_sc26m, fSCA_0p_exp_sc26m, fSCA_0p_shl_sc26m, elev_class_inx)
#statistical test
#fSCA_0u_exp_sc26m_samp, fSCA_0u_shl_sc26m_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx,sample_size)
#statistic_sc26m2,pvalue_sc26m2,average_sample_exp_sc26m2,average_sample_shl_sc26m2,significant_sc26m2 = wilcoxonTest(fSCA_0u_exp_sc26m_samp, fSCA_0u_shl_sc26m_samp)
#
##y = y = -0.0013047442x + 2.0480940252
tempDec_sc26m = -0.0013047442 * elevationBand_sc + + 2.0480940252

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_sc26m = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx)
#fSCA_ut_exp_l_sam_sc26m, fSCA_ut_exp_h_sam_sc26m, fSCA_ut_shl_l_sam_sc26m, fSCA_ut_shl_h_sam_sc26m = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_sc26m,elev_class_inx,sample_size)
#statistic_ex_sc26m,pvalue_ex_sc26m,average_sample_l_ex_sc26m,average_sample_h_ex_sc26m,significant_ex_sc26m = wilcoxonTest(fSCA_ut_exp_l_sam_sc26m, fSCA_ut_exp_h_sam_sc26m)
#statistic_sh_sc26m,pvalue_sh_sc26m,average_sample_l_sh_sc26m,average_sample_h_sh_sc26m,significant_sh_sc26m = wilcoxonTest(fSCA_ut_shl_l_sam_sc26m, fSCA_ut_shl_h_sam_sc26m)
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/significant_sc26m_retile',[significant_sc26m1,significant_sc26m2,significant_ex_sc26m,significant_sh_sc26m])
#
##calculating topo_dimension from cut northness to plot snow on
#snow_temp_vegdens_allExtent_sc26m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sc26m_retile.npy')
#ouputPath_sc26m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_contour_sc26m_retile_lr2.png'
#T_sc26m = 'SCWC_26March'
#plotSnowMap_lowResol(snow_temp_vegdens_allExtent_sc26m,ouputPath_sc26m,T_sc26m)
#
##lat_long_nad83_sc26m = [np.min(snow_temp_vegdens_allExtent_sc26m[:,0]),np.max(snow_temp_vegdens_allExtent_sc26m[:,0]),np.min(snow_temp_vegdens_allExtent_sc26m[:,1]),np.max(snow_temp_vegdens_allExtent_sc26m[:,1])]
#
#ouputPath_nrt_sc26m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_nrth_contour_sc26m_retile_lr2.png'
#plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent_sc26m,ouputPath_nrt_sc26m,T_sc26m)

#%% load data sagehen creek 17 April 2016
#ascii_grid_veg_index17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/ascii_index_pr_sc17A.npy')
##calculating topo_dimension from cut northness
#veg_density_sc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/veg_density_sc.npy')
#
#elev_sc17a = ascii_grid_veg_index17a[:,2]
#tempDJF_sc17a = -0.0016 * elev_sc17a + 1.802
#
#veg_snow_temp_density_sc17a = np.vstack([ascii_grid_veg_index17a[:,0],ascii_grid_veg_index17a[:,1].astype(int),ascii_grid_veg_index17a[:,2],
#                                         ascii_grid_veg_index17a[:,3],ascii_grid_veg_index17a[:,4],ascii_grid_veg_index17a[:,5],
#                                         ascii_grid_veg_index17a[:,6],tempDJF_sc17a,veg_density_sc[:,2],tempDJF_sc17a,tempDJF_sc17a]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sc17a',veg_snow_temp_density_sc17a)

## slope less than 30
#snow_temp_vegdens_index_sL30_sc17a = veg_snow_temp_density_sc17a[(veg_snow_temp_density_sc17a[:,5]>-50)&(veg_snow_temp_density_sc17a[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_sc17a',snow_temp_vegdens_index_sL30_sc17a)

snow_temp_vegdens_index_sL30_sc17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_sc17a.npy')
aaaa_test2 = snow_temp_vegdens_index_sL30_sc17a[0:710,:]

elev_sc17a = snow_temp_vegdens_index_sL30_sc17a[:,2]
tempDJF_sc17a = -0.0016 * elev_sc17a + 1.802

snow_temp_vegdens_index_sL30_sc17a_df = pd.DataFrame(snow_temp_vegdens_index_sL30_sc17a, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])
aaaa_test = snow_temp_vegdens_index_sL30_sc17a[700780:710057,:]

elevationBand_sc17a = np.arange(min(snow_temp_vegdens_index_sL30_sc17a_df['z']),max(snow_temp_vegdens_index_sL30_sc17a_df['z']),67)
snow_temp_vegdens_sL30_elevCls_sc17a = classificationElevation(snow_temp_vegdens_index_sL30_sc17a,elevationBand_sc17a,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc17a, fsca_0p_sc17a = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)
#fsca_ut_samp_sc17a, fsca_0p_samp_sc17a = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx,sample_size)
##statistical test
#statistic_sc17a1,pvalue_sc17a1,average_sample_exp_sc17a1,average_sample_shl_sc17a1,significant_sc17a1 = wilcoxonTest(fsca_ut_samp_sc17a, fsca_0p_samp_sc17a)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a = fsca0pn_fscaUnTr(fSCA_ut_exp_sc17a, fSCA_ut_shl_sc17a, fSCA_0p_exp_sc17a, fSCA_0p_shl_sc17a, elev_class_inx)
##statistical test
#fSCA_0u_exp_sc17a_samp, fSCA_0u_shl_sc17a_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx,sample_size)
#statistic_sc17a2,pvalue_sc17a2,average_sample_exp_sc17a2,average_sample_shl_sc17a2,significant_sc17a2 = wilcoxonTest(fSCA_0u_exp_sc17a_samp, fSCA_0u_shl_sc17a_samp)
#
##y = -0.0013x + 2.7982
tempDec_sc17a = -0.0013 * elevationBand_sc17a + 2.7982

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_sc17a = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx)
#fSCA_ut_exp_l_sam_sc17a, fSCA_ut_exp_h_sam_sc17a, fSCA_ut_shl_l_sam_sc17a, fSCA_ut_shl_h_sam_sc17a = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_sc17a,elev_class_inx,sample_size)
#statistic_ex_sc17a,pvalue_ex_sc17a,average_sample_l_ex_sc17a,average_sample_h_ex_sc17a,significant_ex_sc17a = wilcoxonTest(fSCA_ut_exp_l_sam_sc17a, fSCA_ut_exp_h_sam_sc17a)
#statistic_sh_sc17a,pvalue_sh_sc17a,average_sample_l_sh_sc17a,average_sample_h_sh_sc17a,significant_sh_sc17a = wilcoxonTest(fSCA_ut_shl_l_sam_sc17a, fSCA_ut_shl_h_sam_sc17a)
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/significant_sc17a',[significant_sc17a1,significant_sc17a2,significant_ex_sc17a,significant_sh_sc17a])
#
#snow_temp_vegdens_allExtent_sc17a = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sc17a.npy')
#ouputPath_sc17a = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_contour_sc17a_lr2.png'
#T_sc17a = 'SCWC_17April'
#plotSnowMap_lowResol(snow_temp_vegdens_allExtent_sc17a,ouputPath_sc17a,T_sc17a)
#
##lat_long_nad83_sc17a = [np.min(snow_temp_vegdens_allExtent_sc17a[:,0]),np.max(snow_temp_vegdens_allExtent_sc17a[:,0]),np.min(snow_temp_vegdens_allExtent_sc17a[:,1]),np.max(snow_temp_vegdens_allExtent_sc17a[:,1])]
#
#ouputPath_nrt_sc17a = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_nrth_contour_sc17a_lr2.png'
#plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent_sc17a,ouputPath_nrt_sc17a,T_sc17a)

#%% load data sagehen creek 18 May 2016
#ascii_grid_veg_index18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/ascii_index_pr_sc24May.npy')
##calculating topo_dimension from cut northness
#
#veg_density_sc = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/veg_density_sc.npy')
#
#elev_sc18m = ascii_grid_veg_index18m[:,2]
#tempDJF_sc18m = -0.0016 * elev_sc18m + 1.802
#
#veg_snow_temp_density_sc18m = np.vstack([ascii_grid_veg_index18m[:,0],ascii_grid_veg_index18m[:,1].astype(int),ascii_grid_veg_index18m[:,2],
#                                         ascii_grid_veg_index18m[:,3],ascii_grid_veg_index18m[:,4],ascii_grid_veg_index18m[:,5],
#                                         ascii_grid_veg_index18m[:,6],tempDJF_sc18m,veg_density_sc[:,2],tempDJF_sc18m,tempDJF_sc18m]).T
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sc18m',veg_snow_temp_density_sc18m)

## slope less than 30
#snow_temp_vegdens_index_sL30_sc18m = veg_snow_temp_density_sc18m[(veg_snow_temp_density_sc18m[:,5]>-50)&(veg_snow_temp_density_sc18m[:,4]<30)]
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_sc18m',snow_temp_vegdens_index_sL30_sc18m)

snow_temp_vegdens_index_sL30_sc18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sL30_sc18m.npy')
aaaa_test3 = snow_temp_vegdens_index_sL30_sc18m[0:710,:]

elev_sc18m = snow_temp_vegdens_index_sL30_sc18m[:,2]
tempDJF_sc18m = -0.0016 * elev_sc18m + 1.802

snow_temp_vegdens_index_sL30_sc18m_df = pd.DataFrame(snow_temp_vegdens_index_sL30_sc18m, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','vegDens','indx1','indx2'])

elevationBand_sc18m = np.arange(min(snow_temp_vegdens_index_sL30_sc18m_df['z']),max(snow_temp_vegdens_index_sL30_sc18m_df['z']),67)
snow_temp_vegdens_sL30_elevCls_sc18m = classificationElevation(snow_temp_vegdens_index_sL30_sc18m,elevationBand_sc18m,elev_class_inx)

meanTemp_elevClass_sc = calculationMeanTempforEachElevClassification(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)   

# under canopy and 0pen fsca
fsca_ut_sc18m, fsca_0p_sc18m = fSCAcalculation_elevClassific(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)
#fsca_ut_samp_sc18m, fsca_0p_samp_sc18m = fSCAcalculation_elevRandom(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx,sample_size)
##statistical test
#statistic_sc18m1,pvalue_sc18m1,average_sample_exp_sc18m1,average_sample_shl_sc18m1,significant_sc18m1 = wilcoxonTest(fsca_ut_samp_sc18m, fsca_0p_samp_sc18m)

# under canopy and 0pen fsca---exposed vs. open
fSCA_ut_exp_sc18m, fSCA_ut_shl_sc18m, fSCA_0p_exp_sc18m, fSCA_0p_shl_sc18m = fSCAcalculation_elevNorthClassific(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)
#fsca open -minus- fsca under canopy
fSCA_0u_exp_sc18m, fSCA_0u_shl_sc18m = fsca0pn_fscaUnTr(fSCA_ut_exp_sc18m, fSCA_ut_shl_sc18m, fSCA_0p_exp_sc18m, fSCA_0p_shl_sc18m, elev_class_inx)
##statistical test
#fSCA_0u_exp_sc18m_samp, fSCA_0u_shl_sc18m_samp = fSCAcalculation_elevNorthRandom(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx,sample_size)
#statistic_sc18m2,pvalue_sc18m2,average_sample_exp_sc18m2,average_sample_shl_sc18m2,significant_sc18m2 = wilcoxonTest(fSCA_0u_exp_sc18m_samp, fSCA_0u_shl_sc18m_samp)
#
##y = -0.0020x + 5.2481
tempDec_sc18m = -0.002 * elevationBand_sc18m + 5.2481

# under canopy fsca // exposed vs. sheltered  //  dense vs. spare veg density
fSCA_uT_rad_vegDens_sc18m = fSCAcalculation_elevNorthVegDensClassific(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx)
#fSCA_ut_exp_l_sam_sc18m, fSCA_ut_exp_h_sam_sc18m, fSCA_ut_shl_l_sam_sc18m, fSCA_ut_shl_h_sam_sc18m = fSCAcalculation_elevNorthVegDensRandom(snow_temp_vegdens_sL30_elevCls_sc18m,elev_class_inx,sample_size)
#statistic_ex_sc18m,pvalue_ex_sc18m,average_sample_l_ex_sc18m,average_sample_h_ex_sc18m,significant_ex_sc18m = wilcoxonTest(fSCA_ut_exp_l_sam_sc18m, fSCA_ut_exp_h_sam_sc18m)
#statistic_sh_sc18m,pvalue_sh_sc18m,average_sample_l_sh_sc18m,average_sample_h_sh_sc18m,significant_sh_sc18m = wilcoxonTest(fSCA_ut_shl_l_sam_sc18m, fSCA_ut_shl_h_sam_sc18m)
#np.save('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/significant_sc18m',[significant_sc18m1,significant_sc18m2,significant_ex_sc18m,significant_sh_sc18m])
#
#snow_temp_vegdens_allExtent_sc18m = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snow_temp_vegdens_index_sc18m.npy')
#ouputPath_sc18m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_contour_sc18m_lr2.png'
#T_sc18m = 'SCWC_18May'
#plotSnowMap_lowResol(snow_temp_vegdens_allExtent_sc18m,ouputPath_sc18m,T_sc18m)
#
##lat_long_nad83_sc18m = [np.min(snow_temp_vegdens_allExtent_sc18m[:,0]),np.max(snow_temp_vegdens_allExtent_sc18m[:,0]),np.min(snow_temp_vegdens_allExtent_sc18m[:,1]),np.max(snow_temp_vegdens_allExtent_sc18m[:,1])]
#
#ouputPath_nrt_sc18m = 'C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/snowMap_nrth_contour_sc18m_lr2.png'
#plotSnowMap_northness_lowResol(snow_temp_vegdens_allExtent_sc18m,ouputPath_nrt_sc18m,T_sc18m)


#%% ploting fSCA vs. Dec temp lapse rate
fSCA_0u = [fsca_ut_sc26m, fsca_0p_sc26m, fsca_ut_sc17a, fsca_0p_sc17a,
           fsca_ut_sc18m, fsca_0p_sc18m, fsca_ut_krew, fsca_0p_krew,
           fsca_ut_jmz, fsca_0p_jmz, fsca_ut_nrc, fsca_0p_nrc]

meanTempDec = [tempDec_sc26m,tempDec_sc26m,tempDec_sc17a,tempDec_sc17a,
               tempDec_sc18m,tempDec_sc18m,tempDec_krew,tempDec_krew,
               tempDec_jmz,tempDec_jmz,tempDec_nrc,tempDec_nrc]

label = ['SCWC 26MAR2016-UnderTree','SCWC 26MAR2016-0pen','SCWC 17APR2016-UnderTree','SCWC 17APR2016-0pen',
         'SCWC 18MAY2016-UnderTree','SCWC 18MAY2016-0pen','Krew 2010-UnderTree','Krew 2010-0pen',
         'JMZ 2010-UnderTree','JMZ 2010 - 0pen','NRC 2010-UnderTree','NRC 2010 - 0pen']
color = ['darkred','gold','darkred','gold','darkred','gold',
         'red','orange','darkgreen','lightgreen','navy','lightblue']
marker = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

location = ['lower left','lower left','lower left','lower left','lower left','lower left']
xlimit_fsca = [(-1.5,1.5),(-1.5,1.5),(-1.5,1.5),(0.5,3),(-8,-1.5),(-8,-1.5)]

plt.subplots(figsize=(30,20)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
    for pnt in range (10):
        plt.scatter(meanTempDec[2*i][pnt],fSCA_0u[2*i][pnt], s=25**2, color = color[2*i], marker = marker[i][pnt]) #label = label[2*i],
        plt.scatter(meanTempDec[2*i+1][pnt],fSCA_0u[2*i+1][pnt],  s=25**2, color = color[2*i+1], marker = marker[i][pnt]) #label = label[2*i+1],

    plt.ylabel('fSCA%', fontsize=25)
    plt.yticks(fontsize=15)
    plt.xlabel('average temp from Dec till date of flights(C)', fontsize=20)
    plt.xticks(fontsize=15)
    
    plt.legend([label[2*i], label[2*i+1]], fontsize=20, loc = location[i])
    plt.ylim((0,110))  # adjust the top leaving bottom unchanged
    plt.xlim(xlimit_fsca[i])

plt.title('fSCA under canopy and in open sites in different temp lapse rate in 4 sites', fontsize=35, y=2.25, x=-1) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/tempLsr_dec_fsca_all4.png')
#%% ploting fsca vs. DJF temp lapse rate
fSCA_0u = [fsca_ut_sc26m, fsca_0p_sc26m, fsca_ut_sc17a, fsca_0p_sc17a,
           fsca_ut_sc18m, fsca_0p_sc18m, fsca_ut_krew, fsca_0p_krew,
           fsca_ut_jmz, fsca_0p_jmz, fsca_ut_nrc, fsca_0p_nrc]

meanTemp = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
            meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
            meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label0 = ['SCWC 26MAR2016-UnderTree','SCWC 26MAR2016-0pen','SCWC 17APR2016-UnderTree','SCWC 17APR2016-0pen',
          'SCWC 18MAY2016-UnderTree','SCWC 18MAY2016-0pen','KREW 2010-UnderTree','KREW 2010-0pen',
          'JRBN 2010-UnderTree','JRBN 2010 - 0pen','NRC 2010-UnderTree','NRC 2010 - 0pen']
color0 = ['plum','purple','plum','purple','plum','purple',
          'gold','darkred','lightgreen','darkgreen','deepskyblue','navy']#'olive',
markerS = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

marker0 = ['^','o','^','o','^','o','^','o','^','o','^','o']

location0 = ['lower left','lower left','uper right','lower left','lower left','lower left']
xlimit_fsca_DJF = [(-2.4,-1.2),(-2.4,-1.2),(-2.4,-1.2),(0.4,2.8),(-8.3,-3.5),(-10.5,-5)]

plt.subplots(figsize=(60,30)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
    for pnt in range (10):
        if fSCA_0u[2*i][pnt]>fSCA_0u[2*i+1][pnt]:
            plt.scatter(meanTemp[2*i][pnt],fSCA_0u[2*i][pnt]+7, s=25**2, color = 'k', marker = markerS[i][pnt]) #label = label[2*i],
        else:
            plt.scatter(meanTemp[2*i+1][pnt],fSCA_0u[2*i+1][pnt]+7,  s=25**2, color = 'k', marker = markerS[i][pnt]) #label = label[2*i+1],

    plt.plot(meanTemp[2*i],fSCA_0u[2*i], color = color0[2*i], linewidth=4, marker = marker0[2*i], markersize=40) #, 
    plt.plot(meanTemp[2*i+1],fSCA_0u[2*i+1], color = color0[2*i+1], linewidth=4, marker = marker0[2*i+1], markersize=40) #

    plt.ylabel('fSCA%', fontsize=50)
    plt.yticks(fontsize=40)
    plt.xlabel('average temp of DJF (C)', fontsize=40)
    plt.xticks(fontsize=40)
    
    plt.legend([label0[2*i], label0[2*i+1]], fontsize=35, loc = location0[i])
    plt.ylim((0,110))  # adjust the top leaving bottom unchanged
    plt.xlim(xlimit_fsca_DJF[i])

plt.title('fSCA under canopy and in open sites in different temp lapse rate in 4 sites', fontsize=70, y=2.25, x=-0.8) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/tempLsr_jfd_fsca_all_retile2.png')

#%% ploting DJF temp vs terrain features delta fsca
fSCA_d0u = [fSCA_0u_exp_sc26m, fSCA_0u_shl_sc26m,fSCA_0u_exp_sc17a, fSCA_0u_shl_sc17a,  
            fSCA_0u_exp_sc18m, fSCA_0u_shl_sc18m,fSCA_0u_exp_Krew, fSCA_0u_shl_Krew,  
            fSCA_0u_exp_Jmz, fSCA_0u_shl_Jmz,fSCA_0u_exp_nrc, fSCA_0u_shl_nrc]
meanTemp2 = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
             meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label = ['SCWC 26Mar2016 - exposed','SCWC 26Mar2016 - sheltered',
         'SCWC 17Apr2016 - exposed','SCWC 17Apr2016 - sheltered',
         'SCWC 18May2016 - exposed','SCWC 18May2016 - sheltered',
         'JRBN 2010 - exposed','JRBN 2010 - sheltered',         
         'KREW 2010 - exposed','KREW 2010 - sheltered',
         'NRC 2010 - exposed','NRC 2010 - sheltered']

color1 = ['plum','purple','plum','purple','plum','purple','gold','darkred','lightgreen','darkgreen','deepskyblue','navy']#'olive',
#marker = [significant_sc26m1,significant_sc17a1,significant_sc18m1,significant_krew1,significant_jmz1,significant_nrc1]

marker1 = ['o','^','o','^','o','^','o','^','o','^','o','^']

markerS2 = [significant_sc26m2,significant_sc17a2,significant_sc18m2,significant_krew2,significant_jmz2,significant_nrc2]

location = ['lower left','lower left','lower left','lower left','lower left','lower left']
xlimit_deltAfsca = [(-2.4,-1.2),(-2.4,-1.2),(-2.4,-1.2),(0.4,2.8),(-8.3,-3.5),(-10.5,-5)]
ylimit = [(-1.2,0.8),(-1.2,0.8),(-1.2,0.8),(-0.9,0.8),(-1.2,0.8),(-1.2,0.8)]

plt.subplots(figsize=(60,40)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
#    for pnt in range (10):
#        plt.scatter(meanTemp2[2*i][pnt],fSCA_d0u[2*i][pnt], s=50**2, color = color1[2*i], marker = marker1[2*i]) #label = label[2*i],
#        plt.scatter(meanTemp2[2*i+1][pnt],fSCA_d0u[2*i+1][pnt],  s=50**2, color = color1[2*i+1], marker = marker1[2*i+1]) #label = label[2*i+1],
    for pnt in range (10):
        if fSCA_d0u[2*i][pnt]>fSCA_d0u[2*i+1][pnt]:
            plt.scatter(meanTemp2[2*i][pnt],fSCA_d0u[2*i][pnt]+0.1, s=25**2, color = 'k', marker = markerS2[i][pnt]) #label = label[2*i],
        else:
            plt.scatter(meanTemp2[2*i+1][pnt],fSCA_d0u[2*i+1][pnt]+0.1,  s=25**2, color = 'k', marker = markerS2[i][pnt]) #label = label[2*i+1],

    plt.plot(meanTemp2[2*i],fSCA_d0u[2*i], color = color1[2*i], linewidth=4, marker = marker1[2*i], markersize=40) #, 
    plt.plot(meanTemp2[2*i+1],fSCA_d0u[2*i+1], color = color1[2*i+1], linewidth=4, marker = marker1[2*i+1], markersize=40) #

    plt.ylabel('(fSCA_open - fSCA_underTree) / fSCA_open', fontsize=45)
    plt.yticks(fontsize=30)
    plt.xlabel('average temp of DJF (C)', fontsize=45)
    plt.xticks(fontsize=30)
    
    plt.legend([label[2*i], label[2*i+1]], fontsize=40, loc = location[i])
    plt.ylim(ylimit[i])
    plt.xlim(xlimit_deltAfsca[i])
    
    x = [-11,-8,-5,-2,0,2,4]
    y = [0,0,0,0,0,0,0]
    plt.plot(x,y, color = 'black') #line_temp_ls[2*i],
    
    #plt.arrow(arrow_loc[i], 0.05, 0, 0.5, fc="k", ec="k", head_width=0.1, head_length=0.1)

plt.title('(fSCA_open - fSCA_underTree)/fSCA_open based on northness in 4 sites', fontsize=60, y=2.24, x=-0.9) # loc = 'right', 

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/tempLsr_nrth_deltaFsca_all8_retile.png')

#%% ploting DJF temp vs terrain features delta fsca
#fSCA_uT_rad_vegDens_krew = [fSCA_unTr_exp_l, fSCA_unTr_shl_l, fSCA_unTr_exp_h, fSCA_unTr_shl_h]

fSCA_0uV = [fSCA_uT_rad_vegDens_sc26m[0], fSCA_uT_rad_vegDens_sc26m[1], fSCA_uT_rad_vegDens_sc26m[2], fSCA_uT_rad_vegDens_sc26m[3], 
            fSCA_uT_rad_vegDens_sc17a[0], fSCA_uT_rad_vegDens_sc17a[1], fSCA_uT_rad_vegDens_sc17a[2], fSCA_uT_rad_vegDens_sc17a[3],
            fSCA_uT_rad_vegDens_sc18m[0], fSCA_uT_rad_vegDens_sc18m[1], fSCA_uT_rad_vegDens_sc18m[2], fSCA_uT_rad_vegDens_sc18m[3],
            fSCA_uT_rad_vegDens_krew[0], fSCA_uT_rad_vegDens_krew[1], fSCA_uT_rad_vegDens_krew[2], fSCA_uT_rad_vegDens_krew[3],
            fSCA_uT_rad_vegDens_jmz[0], fSCA_uT_rad_vegDens_jmz[1], fSCA_uT_rad_vegDens_jmz[2], fSCA_uT_rad_vegDens_jmz[3], 
            fSCA_uT_rad_vegDens_nrc[0], fSCA_uT_rad_vegDens_nrc[1], fSCA_uT_rad_vegDens_nrc[2], fSCA_uT_rad_vegDens_nrc[3]]

meanTemp2 = [meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,meanTemp_elevClass_sc,
             meanTemp_elevClass_krew,meanTemp_elevClass_krew,meanTemp_elevClass_krew,meanTemp_elevClass_krew,
             meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,meanTemp_elevClass_jmz,
             meanTemp_elevClass_nrc,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc,meanTemp_elevClass_nrc]

label = ['SCWC26Mar2016-exposed-lowVD','SCWC26Mar2016-shelter-lowVD',
         'SCWC26Mar2016-exposed-highVD','SCWC26Mar2016-shelter-highVD',
         'SCWC17Apr2016-exp-lowVD','SCWC17Apr2016-shel-lowVD',
         'SCWC17Apr2016-exp-highVD','SCWC17Apr2016-shel-highVD',
         'SCWC18May2016-exposed-lowVD','SCWC18May2016-shelter-lowVD',
         'SCWC18May2016-exposed-highVD','SCWC18May2016-shelter-highVD',
         'KREW2010-exposed-lowVD','KREW2010-shelter-lowVD',
         'KREW2010-exposed-highVD','KREW2010-shelter-highVD',
         'JRBN2010-exposed-lowVD','JRBN2010-shelter-lowVD',
         'JRBN2010-exposed-highVD','JRBN2010-shelter-highVD',
         'NRC2010-exposed-lowVD','NRC2010-shelter-lowVD',
         'NRC2010-exposed-highVD','NRC2010-shelter-highVD']

color = ['orchid','purple','plum','darkorchid',
         'orchid','purple','plum','darkorchid',
         'orchid','purple','plum','darkorchid',
         'red','darkred','gold','orange',
         'lightgreen','darkgreen','olive','green',
         'deepskyblue','navy','lightblue','blue']

marker = ['o','s','^','*',
          'o','s','^','*',
          'o','s','^','*',
          'o','s','^','*',
          'o','s','^','*',
          'o','s','^','*']

#markerS3 = [significant_ex_sc26m,significant_ex_sc17a,significant_ex_sc18m,significant_ex_krew,significant_ex_nrc,significant_ex_jmz]
#markerS4 = [significant_sh_sc26m,significant_sh_sc17a,significant_sh_sc18m,significant_sh_krew,significant_sh_nrc,significant_sh_jmz]


location = ['lower left','lower left','uper right','lower left','lower left','lower left']
xlimit_deltAfsca = [(-2.4,-1.2),(-2.4,-1.2),(-2.4,-1.2),(0.4,2.8),(-8.3,-3.5),(-10.5,-5)]

plt.subplots(figsize=(60,40)) #fig, ax = plt.subplots(figsize=(20,15))

for i in range (6):
    
    plt.subplot(231+i)
    
#    for pnt in range (10):
#        plt.scatter(meanTemp2[4*i][pnt],fSCA_d0u[4*i][pnt], s=30**2, color = color[4*i], marker = marker[4*i]) #label = label[2*i],
#        plt.scatter(meanTemp2[4*i+1][pnt],fSCA_d0u[4*i+1][pnt],  s=25**2, color = color[4*i+1], marker = marker[4*i+1]) #label = label[2*i+1],
#        plt.scatter(meanTemp2[4*i+2][pnt],fSCA_d0u[4*i+2][pnt],  s=35**2, color = color[4*i+2], marker = marker[4*i+2]) #label = label[2*i+1],
#        plt.scatter(meanTemp2[4*i+3][pnt],fSCA_d0u[4*i+3][pnt],  s=40**2, color = color[4*i+3], marker = marker[4*i+3]) #label = label[2*i+1],

#    for pnt in range (10):
#        if fSCA_0uV[4*i][pnt]>fSCA_0uV[4*i+2][pnt]:
#            plt.scatter(meanTemp2[4*i][pnt],fSCA_0uV[4*i][pnt]+7, s=25**2, color = 'k', marker = markerS3[i][pnt]) #label = label[2*i],
#        else:
#            plt.scatter(meanTemp2[4*i+2][pnt],fSCA_0uV[4*i+2][pnt]+7,  s=25**2, color = 'k', marker = markerS3[i][pnt]) #label = label[2*i+1],

    plt.plot(meanTemp2[4*i],fSCA_0uV[4*i], color = color[4*i], linewidth=4, marker = marker[4*i], markersize=30) #, 
    plt.plot(meanTemp2[4*i+1],fSCA_0uV[4*i+1], color = color[4*i+1], linewidth=4, marker = marker[4*i+1], markersize=25) #
    plt.plot(meanTemp2[4*i+2],fSCA_0uV[4*i+2], color = color[4*i+2], linewidth=4, marker = marker[4*i+2], markersize=35) #
    plt.plot(meanTemp2[4*i+3],fSCA_0uV[4*i+3], color = color[4*i+3], linewidth=4, marker = marker[4*i+3], markersize=40) #

    plt.ylabel('fSCA%', fontsize=50)
    plt.yticks(fontsize=50)
    plt.xlabel('average temp of DJF (C)', fontsize=50)
    plt.xticks(fontsize=40)
    
    plt.legend([label[4*i], label[4*i+1], label[4*i+2], label[4*i+3]], fontsize=30, loc = location[i])
    plt.ylim((0,110))
    plt.xlim(xlimit_deltAfsca[i])
    

plt.title('fSCA classification based on northness and  vegetation density (VD)', fontsize=60, y=2.24, x=-0.65) # loc = 'right', 

#plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/tempLsr_nrth_vegDens_all8_retile.png')

#%%
#fsca_0pUt_shl2550 = fsca0p_fscaUt_fscaLimit (fSCA_0p_shl_sc26m,fSCA_0p_shl_sc17a,fSCA_0p_shl_sc18m,fSCA_0p_shl_Krew,fSCA_0p_shl_Jmz,fSCA_0p_shl_nrc,
#                                             fSCA_ut_shl_sc26m,fSCA_ut_shl_sc17a,fSCA_ut_shl_sc18m,fSCA_ut_shl_Krew,fSCA_ut_shl_Jmz,fSCA_ut_shl_nrc,25,50)
#fsca_0pUt_exp2550 = fsca0p_fscaUt_fscaLimit (fSCA_0p_exp_sc26m,fSCA_0p_exp_sc17a,fSCA_0p_exp_sc18m,fSCA_0p_exp_Krew,fSCA_0p_exp_Jmz,fSCA_0p_exp_nrc,
#                                             fSCA_ut_exp_sc26m,fSCA_ut_exp_sc17a,fSCA_ut_exp_sc18m,fSCA_ut_exp_Krew,fSCA_ut_exp_Jmz,fSCA_ut_exp_nrc,25,50)

fsca_0t2550 = fsca0p_fscaUt_fscaLimit (fsca_0p_sc26m,fsca_0p_sc17a,fsca_0p_sc18m,fsca_0p_krew,fsca_0p_jmz,fsca_0p_nrc,
                                       fsca_ut_sc26m,fsca_ut_sc17a,fsca_ut_sc18m,fsca_ut_krew,fsca_ut_jmz,fsca_ut_nrc,25,50)
fsca_0t5075 = fsca0p_fscaUt_fscaLimit (fsca_0p_sc26m,fsca_0p_sc17a,fsca_0p_sc18m,fsca_0p_krew,fsca_0p_jmz,fsca_0p_nrc,
                                       fsca_ut_sc26m,fsca_ut_sc17a,fsca_ut_sc18m,fsca_ut_krew,fsca_ut_jmz,fsca_ut_nrc,50,75)
fsca_0t7599 = fsca0p_fscaUt_fscaLimit (fsca_0p_sc26m,fsca_0p_sc17a,fsca_0p_sc18m,fsca_0p_krew,fsca_0p_jmz,fsca_0p_nrc,
                                       fsca_ut_sc26m,fsca_ut_sc17a,fsca_ut_sc18m,fsca_ut_krew,fsca_ut_jmz,fsca_ut_nrc,75,99)

fsca_shEx_0p2550 = fsca0p_fscaUt_fscaLimit (fSCA_0p_shl_sc26m,fSCA_0p_shl_sc17a,fSCA_0p_shl_sc18m,fSCA_0p_shl_Krew,fSCA_0p_shl_Jmz,fSCA_0p_shl_nrc,
                                            fSCA_0p_exp_sc26m,fSCA_0p_exp_sc17a,fSCA_0p_exp_sc18m,fSCA_0p_exp_Krew,fSCA_0p_exp_Jmz,fSCA_0p_exp_nrc,25,50)
fsca_shEx_ut2550 = fsca0p_fscaUt_fscaLimit (fSCA_ut_shl_sc26m,fSCA_ut_shl_sc17a,fSCA_ut_shl_sc18m,fSCA_ut_shl_Krew,fSCA_ut_shl_Jmz,fSCA_ut_shl_nrc,
                                            fSCA_ut_exp_sc26m,fSCA_ut_exp_sc17a,fSCA_ut_exp_sc18m,fSCA_ut_exp_Krew,fSCA_ut_exp_Jmz,fSCA_ut_exp_nrc,25,50)
fsca_shEx_0p5075 = fsca0p_fscaUt_fscaLimit (fSCA_0p_shl_sc26m,fSCA_0p_shl_sc17a,fSCA_0p_shl_sc18m,fSCA_0p_shl_Krew,fSCA_0p_shl_Jmz,fSCA_0p_shl_nrc,
                                            fSCA_0p_exp_sc26m,fSCA_0p_exp_sc17a,fSCA_0p_exp_sc18m,fSCA_0p_exp_Krew,fSCA_0p_exp_Jmz,fSCA_0p_exp_nrc,50,75)
fsca_shEx_ut5075 = fsca0p_fscaUt_fscaLimit (fSCA_ut_shl_sc26m,fSCA_ut_shl_sc17a,fSCA_ut_shl_sc18m,fSCA_ut_shl_Krew,fSCA_ut_shl_Jmz,fSCA_ut_shl_nrc,
                                            fSCA_ut_exp_sc26m,fSCA_ut_exp_sc17a,fSCA_ut_exp_sc18m,fSCA_ut_exp_Krew,fSCA_ut_exp_Jmz,fSCA_ut_exp_nrc,50,75)
fsca_shEx_0p7599 = fsca0p_fscaUt_fscaLimit (fSCA_0p_shl_sc26m,fSCA_0p_shl_sc17a,fSCA_0p_shl_sc18m,fSCA_0p_shl_Krew,fSCA_0p_shl_Jmz,fSCA_0p_shl_nrc,
                                            fSCA_0p_exp_sc26m,fSCA_0p_exp_sc17a,fSCA_0p_exp_sc18m,fSCA_0p_exp_Krew,fSCA_0p_exp_Jmz,fSCA_0p_exp_nrc,75,99)
fsca_shEx_ut7599 = fsca0p_fscaUt_fscaLimit (fSCA_ut_shl_sc26m,fSCA_ut_shl_sc17a,fSCA_ut_shl_sc18m,fSCA_ut_shl_Krew,fSCA_ut_shl_Jmz,fSCA_ut_shl_nrc,
                                            fSCA_ut_exp_sc26m,fSCA_ut_exp_sc17a,fSCA_ut_exp_sc18m,fSCA_ut_exp_Krew,fSCA_ut_exp_Jmz,fSCA_ut_exp_nrc,75,99)

fsca_hlvd_exp2550 = fsca0p_fscaUt_fscaLimit (fSCA_uT_rad_vegDens_sc26m[2],fSCA_uT_rad_vegDens_sc17a[2],fSCA_uT_rad_vegDens_sc18m[2],
                                             fSCA_uT_rad_vegDens_krew[2],fSCA_uT_rad_vegDens_jmz[2],fSCA_uT_rad_vegDens_nrc[2],
                                             fSCA_uT_rad_vegDens_sc26m[0],fSCA_uT_rad_vegDens_sc17a[0],fSCA_uT_rad_vegDens_sc18m[0],
                                             fSCA_uT_rad_vegDens_krew[0],fSCA_uT_rad_vegDens_jmz[0],fSCA_uT_rad_vegDens_nrc[0],25,50)

fsca_hlvd_shl2550 = fsca0p_fscaUt_fscaLimit (fSCA_uT_rad_vegDens_sc26m[3],fSCA_uT_rad_vegDens_sc17a[3],fSCA_uT_rad_vegDens_sc18m[3],
                                             fSCA_uT_rad_vegDens_krew[3],fSCA_uT_rad_vegDens_jmz[3],fSCA_uT_rad_vegDens_nrc[3],
                                             fSCA_uT_rad_vegDens_sc26m[1],fSCA_uT_rad_vegDens_sc17a[1],fSCA_uT_rad_vegDens_sc18m[1],
                                             fSCA_uT_rad_vegDens_krew[1],fSCA_uT_rad_vegDens_jmz[1],fSCA_uT_rad_vegDens_nrc[1],25,50)

fsca_hlvd_exp5075 = fsca0p_fscaUt_fscaLimit (fSCA_uT_rad_vegDens_sc26m[2],fSCA_uT_rad_vegDens_sc17a[2],fSCA_uT_rad_vegDens_sc18m[2],
                                             fSCA_uT_rad_vegDens_krew[2],fSCA_uT_rad_vegDens_jmz[2],fSCA_uT_rad_vegDens_nrc[2],
                                             fSCA_uT_rad_vegDens_sc26m[0],fSCA_uT_rad_vegDens_sc17a[0],fSCA_uT_rad_vegDens_sc18m[0],
                                             fSCA_uT_rad_vegDens_krew[0],fSCA_uT_rad_vegDens_jmz[0],fSCA_uT_rad_vegDens_nrc[0],50,75)

fsca_hlvd_shl5075 = fsca0p_fscaUt_fscaLimit (fSCA_uT_rad_vegDens_sc26m[3],fSCA_uT_rad_vegDens_sc17a[3],fSCA_uT_rad_vegDens_sc18m[3],
                                             fSCA_uT_rad_vegDens_krew[3],fSCA_uT_rad_vegDens_jmz[3],fSCA_uT_rad_vegDens_nrc[3],
                                             fSCA_uT_rad_vegDens_sc26m[1],fSCA_uT_rad_vegDens_sc17a[1],fSCA_uT_rad_vegDens_sc18m[1],
                                             fSCA_uT_rad_vegDens_krew[1],fSCA_uT_rad_vegDens_jmz[1],fSCA_uT_rad_vegDens_nrc[1],50,75)

fsca_hlvd_exp7599 = fsca0p_fscaUt_fscaLimit (fSCA_uT_rad_vegDens_sc26m[2],fSCA_uT_rad_vegDens_sc17a[2],fSCA_uT_rad_vegDens_sc18m[2],
                                             fSCA_uT_rad_vegDens_krew[2],fSCA_uT_rad_vegDens_jmz[2],fSCA_uT_rad_vegDens_nrc[2],
                                             fSCA_uT_rad_vegDens_sc26m[0],fSCA_uT_rad_vegDens_sc17a[0],fSCA_uT_rad_vegDens_sc18m[0],
                                             fSCA_uT_rad_vegDens_krew[0],fSCA_uT_rad_vegDens_jmz[0],fSCA_uT_rad_vegDens_nrc[0],75,99)

fsca_hlvd_shl7599 = fsca0p_fscaUt_fscaLimit (fSCA_uT_rad_vegDens_sc26m[3],fSCA_uT_rad_vegDens_sc17a[3],fSCA_uT_rad_vegDens_sc18m[3],
                                             fSCA_uT_rad_vegDens_krew[3],fSCA_uT_rad_vegDens_jmz[3],fSCA_uT_rad_vegDens_nrc[3],
                                             fSCA_uT_rad_vegDens_sc26m[1],fSCA_uT_rad_vegDens_sc17a[1],fSCA_uT_rad_vegDens_sc18m[1],
                                             fSCA_uT_rad_vegDens_krew[1],fSCA_uT_rad_vegDens_jmz[1],fSCA_uT_rad_vegDens_nrc[1],75,99)

meanTempDec = [tempDec_sc26m,tempDec_sc26m,tempDec_sc17a,tempDec_sc17a,
               tempDec_sc18m,tempDec_sc18m,tempDec_krew,tempDec_krew,
               tempDec_jmz,tempDec_jmz,tempDec_nrc,tempDec_nrc]

fsca_tempDec_sc26m = np.array([np.array((tempDec_sc26m[0:10]+tempDec_sc26m[1:11])/2.),(np.array(fsca_0p_sc26m)/2+np.array(fsca_ut_sc26m)/2)])
tempDec_75_99_sc26m = np.mean(fsca_tempDec_sc26m[0][(fsca_tempDec_sc26m[1]<=99)&(fsca_tempDec_sc26m[1]>=75)])
tempDec_50_75_sc26m = np.mean(fsca_tempDec_sc26m[0][(fsca_tempDec_sc26m[1]<=75)&(fsca_tempDec_sc26m[1]>=50)])
tempDec_25_50_sc26m = np.mean(fsca_tempDec_sc26m[0][(fsca_tempDec_sc26m[1]<=50)&(fsca_tempDec_sc26m[1]>=25)])

fsca_tempDec_sc17a = np.array([np.array((tempDec_sc17a[0:10]+tempDec_sc17a[1:11])/2.),(np.array(fsca_0p_sc17a)/2+np.array(fsca_ut_sc17a)/2)])
tempDec_75_99_sc17a = np.mean(fsca_tempDec_sc17a[0][(fsca_tempDec_sc17a[1]<=99)&(fsca_tempDec_sc17a[1]>=75)])
tempDec_50_75_sc17a = np.mean(fsca_tempDec_sc17a[0][(fsca_tempDec_sc17a[1]<=75)&(fsca_tempDec_sc17a[1]>=50)])
tempDec_25_50_sc17a = np.mean(fsca_tempDec_sc17a[0][(fsca_tempDec_sc17a[1]<=50)&(fsca_tempDec_sc17a[1]>=25)])

fsca_tempDec_sc18m = np.array([np.array((tempDec_sc18m[0:10]+tempDec_sc18m[1:11])/2.),(np.array(fsca_0p_sc18m)/2+np.array(fsca_ut_sc18m)/2)])
tempDec_75_99_sc18m = np.mean(fsca_tempDec_sc18m[0][(fsca_tempDec_sc18m[1]<=99)&(fsca_tempDec_sc18m[1]>=75)])
tempDec_50_75_sc18m = np.mean(fsca_tempDec_sc18m[0][(fsca_tempDec_sc18m[1]<=75)&(fsca_tempDec_sc18m[1]>=50)])
tempDec_25_50_sc18m = np.mean(fsca_tempDec_sc18m[0][(fsca_tempDec_sc18m[1]<=50)&(fsca_tempDec_sc18m[1]>=25)])

fsca_tempDec_krew = np.array([np.array((tempDec_krew[0:10]+tempDec_krew[1:11])/2.),(np.array(fsca_0p_krew)/2+np.array(fsca_ut_krew)/2)])
tempDec_75_99_krew = np.mean(fsca_tempDec_krew[0][(fsca_tempDec_krew[1]<=99)&(fsca_tempDec_krew[1]>=75)])
tempDec_50_75_krew = np.mean(fsca_tempDec_krew[0][(fsca_tempDec_krew[1]<=75)&(fsca_tempDec_krew[1]>=50)])
tempDec_25_50_krew = np.mean(fsca_tempDec_krew[0][(fsca_tempDec_krew[1]<=50)&(fsca_tempDec_krew[1]>=25)])

fsca_tempDec_jmz = np.array([np.array((tempDec_jmz[0:10]+tempDec_jmz[1:11])/2.),(np.array(fsca_0p_jmz)/2+np.array(fsca_ut_jmz)/2)])
tempDec_75_99_jmz = np.mean(fsca_tempDec_jmz[0][(fsca_tempDec_jmz[1]<=99)&(fsca_tempDec_jmz[1]>=75)])
tempDec_50_75_jmz = np.mean(fsca_tempDec_jmz[0][(fsca_tempDec_jmz[1]<=75)&(fsca_tempDec_jmz[1]>=50)])
tempDec_25_50_jmz = np.mean(fsca_tempDec_jmz[0][(fsca_tempDec_jmz[1]<=50)&(fsca_tempDec_jmz[1]>=25)])

fsca_tempDec_nrc = np.array([np.array((tempDec_nrc[0:10]+tempDec_nrc[1:11])/2.),(np.array(fsca_0p_nrc)/2+np.array(fsca_ut_nrc)/2)])
tempDec_75_99_nrc = np.mean(fsca_tempDec_nrc[0][(fsca_tempDec_nrc[1]<=99)&(fsca_tempDec_nrc[1]>=75)])
tempDec_50_75_nrc = np.mean(fsca_tempDec_nrc[0][(fsca_tempDec_nrc[1]<=75)&(fsca_tempDec_nrc[1]>=50)])
tempDec_25_50_nrc = np.mean(fsca_tempDec_nrc[0][(fsca_tempDec_nrc[1]<=50)&(fsca_tempDec_nrc[1]>=25)])

#%%bar plot open-under

ind = np.arange(0,6)  # the x locations for the groups
width = 0.2  # the width of the bars
plt.subplots(figsize=(60,40))
plt.bar(ind - width, fsca_0t7599, width, #yerr=women_std,
        color='skyblue', label='fSCA_75_99') #palegoldenrod  goldenrod darkgoldenrod
plt.bar(ind , fsca_0t5075, width, #yerr=men_std,
         color='forestgreen', label='fSCA_50_75') # palevioletred mediumvioletred firebrick
plt.bar(ind + width, fsca_0t2550, width, #yerr=women_std,
        color='rosybrown', label='fSCA_25_50') #palegoldenrod  goldenrod darkgoldenrod

plt.ylabel('fSCA_open - fSCA_underTree', fontsize=80) #(fSCA_op - fSCA_ut) in sheltered
plt.yticks(fontsize=80)
plt.xlabel('Sites', fontsize=90) #(fSCA_op - fSCA_ut) in exposed
plt.xticks(ind, ('SCWC_Mar', 'SCWC_Mar', 'SCWC_Mar', 'KREW', 'JRBN', 'NRC'),fontsize=70)

x = [-1,-8,-5,-2,0,2,10]
y = [0,0,0,0,0,0,0]
plt.plot(x,y, color = 'black',linewidth=4)

marker = ['-1C','-1C','-1C','+2C','-4C','-6C']#
ytemp = [4.5,2.3,7,2.7,6.7,4.1]
xtemp = ind-0.2
for t in range (6):
    plt.text(xtemp[t],ytemp[t],str(marker[t]),fontsize=90)

nanX = []
for na in range (6):
    if np.isnan(fsca_0t2550[na]):
        nanX.append(ind[na] + width)

plt.bar(nanX[0], 1.5, width, #yerr=women_std,
        color='white', edgecolor='black', hatch='/', label='nan') 
plt.bar(nanX[1], 1.5, width, #yerr=women_std,
        color='white', edgecolor='black', hatch='/')
plt.bar(nanX[0], -1.5, width, #yerr=women_std,
        color='white', edgecolor='black', hatch='/') 
plt.bar(nanX[1], -1.5, width, #yerr=women_std,
        color='white', edgecolor='black', hatch='/') 
  
plt.ylim((-9,9))
plt.xlim((-0.9,6))
plt.legend(fontsize=55, loc = 'upper left')

#plt.title('50% to 75% fSCA', fontsize=60) #, y=2.24, x=-0.65 loc = 'right', 
 
#plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_opUt_25507599.png')

#%%bar plot open-under exp-shl
width = 0.1  # the width of the bars

plt.subplots(figsize=(60,40))

plt.bar(ind - 3*width, fsca_shEx_0p7599, width,color='darkred', label='open_75_99') 
plt.bar(ind - 2*width, fsca_shEx_0p5075, width,color='mediumvioletred', label='open_50_75') 
plt.bar(ind - width, fsca_shEx_0p2550, width, color='pink', label='open_25_50') 

plt.bar(ind , fsca_shEx_ut7599, width, color='darkgreen', label='underTree_75_99') 
plt.bar(ind + width, fsca_shEx_ut5075, width, color='darkseagreen', label='underTree_50_75') 
plt.bar(ind + 2*width, fsca_shEx_ut2550, width, color='lightgreen', label='underTree_25_50') 

plt.ylabel('fSCA_sheltered - fSCA_exposed', fontsize=80) #(fSCA_op - fSCA_ut) in sheltered
plt.yticks(fontsize=80)
plt.xlabel('Sites', fontsize=90) #(fSCA_op - fSCA_ut) in exposed
plt.xticks(ind, ('SCWC_Mar', 'SCWC_Mar', 'SCWC_Mar', 'KREW', 'JRBN', 'NRC'),fontsize=70)

x = [-1,-8,-5,-2,0,2,10]
y = [0,0,0,0,0,0,0]
plt.plot(x,y, color = 'black',linewidth=4)

marker = ['-1C','-1C','-1C','+2C','-4C','-6C']#
ytemp = [-13,-13,-13,-13,-13,-13]
xtemp = ind-0.2
for t in range (6):
    plt.text(xtemp[t],ytemp[t],str(marker[t]),fontsize=80)

nanX1 = [-0.2,-0.1,0.1,0.2,0.9,1.7,2,2.8,3.9,4.2]

plt.bar(nanX1[0], 1.5, width, color='white', edgecolor='black', hatch='/', label='nan') 
for bars in range (len(nanX1)):
    plt.bar(nanX1[bars], 1.5, width, color='white', edgecolor='black', hatch='/') 
for bars in range (len(nanX1)):
    plt.bar(nanX1[bars], -1.5, width, color='white', edgecolor='black', hatch='/') 
    
plt.ylim((-14,18))
plt.xlim((-0.9,6))
plt.legend(fontsize=55, loc = 'upper right')

#plt.title('50% to 75% fSCA', fontsize=60) #, y=2.24, x=-0.65 loc = 'right', 
 
#plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_opUt_shlExp_25507599.png')

#%%bar plot hvd-lvd exp-shl
fsca_hlvd_shl5075 = np.array([np.nan, -11.95366529, 4.14264728, np.nan, 11.16302138, 0.2])
fsca_hlvd_exp2550 = np.array([-8.20500649, 8.99780804, np.nan, 8.77211629, np.nan, 0.2])

plt.subplots(figsize=(60,40))
plt.bar(ind - 3*width, fsca_hlvd_exp7599, width,color='darkgreen', label='exposed_75_99')
plt.bar(ind - 2*width, fsca_hlvd_exp5075, width,color='darkseagreen', label='exposed_50_75')
plt.bar(ind - width, fsca_hlvd_exp2550, width, color='lightgreen', label='exposed_25_50') 

plt.bar(ind , fsca_hlvd_shl7599, width,color='navy', label='exposed_75_99')
plt.bar(ind + width, fsca_hlvd_shl5075, width, color='deepskyblue', label='sheltered_50_75') 
plt.bar(ind + 2*width, fsca_hlvd_shl2550, width, color='powderblue', label='sheltered_25_50') 

plt.ylabel('fSCA_highVD - fSCA_lowVD', fontsize=80) #(fSCA_op - fSCA_ut) in sheltered
plt.yticks(fontsize=80)
plt.xlabel('Sites', fontsize=90) #(fSCA_op - fSCA_ut) in exposed
plt.xticks(ind, ('SCWC_Mar', 'SCWC_Mar', 'SCWC_Mar', 'KREW', 'JRBN', 'NRC'),fontsize=70)

x = [-1,-8,-5,-2,0,2,10]
y = [0,0,0,0,0,0,0]
plt.plot(x,y, color = 'black',linewidth=2)

marker = ['-1C','-1C','-1C','+2C','-4C','-6C']#
ytemp = [-13,-13,-13,-13,-13,-13]
xtemp = ind-0.2
for t in range (6):
    plt.text(xtemp[t],ytemp[t],str(marker[t]),fontsize=80)

nanX2 = [0.1,0.2,1.2,1.7,1.9,2.8,3.1,3.9,4.2,4.8]

plt.bar(nanX2[0], 1.5, width, color='white', edgecolor='black', hatch='/', label='nan') 
for bars in range (len(nanX1)):
    plt.bar(nanX2[bars], 1.5, width, color='white', edgecolor='black', hatch='/') 
for bars in range (len(nanX1)):
    plt.bar(nanX2[bars], -1.5, width, color='white', edgecolor='black', hatch='/') 

plt.ylim((-14,18))
plt.xlim((-0.9,6))
plt.legend(fontsize=55, loc = 'upper left')

#plt.title('50% to 75% fSCA', fontsize=60) #, y=2.24, x=-0.65 loc = 'right', 
 
#plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_hlvd_shlExp_25507599.png')

#%% scatter plot

tempDec7599 = [tempDec_75_99_sc26m,tempDec_75_99_sc17a,tempDec_75_99_sc18m,
               tempDec_75_99_krew,tempDec_75_99_jmz,tempDec_75_99_nrc]
tempDec5075 = [tempDec_50_75_sc26m,tempDec_50_75_sc17a,tempDec_50_75_sc18m,
               tempDec_50_75_krew,tempDec_50_75_jmz,tempDec_50_75_nrc]
tempDec2550 = [tempDec_25_50_sc26m,tempDec_25_50_sc17a,tempDec_25_50_sc18m,
               tempDec_25_50_krew,tempDec_25_50_jmz,tempDec_25_50_nrc]
A0t7599 = np.array([tempDec7599,fsca_0t7599])
A0t5075 = np.array([tempDec5075,fsca_0t5075])
A0t2550 = np.array([tempDec2550,fsca_0t2550])

temp_dec = []
temp_dec.extend(tempDec7599); temp_dec.extend(tempDec5075); temp_dec.extend(tempDec2550);
fsca_0t = []
fsca_0t.extend(fsca_0t7599); fsca_0t.extend(fsca_0t5075); fsca_0t.extend(fsca_0t2550);

sns_0t = pd.DataFrame(np.array([temp_dec,fsca_0t]).T,columns=['temp','fsca'])
sns_0t.dropna(axis=0, inplace = True)
sns_0t.index = np.arange(0,len(sns_0t))

plt.subplots(figsize=(60,40))
plt.plot(A0t7599[0,A0t7599[0,:].argsort()],A0t7599[1,A0t7599[0,:].argsort()],linestyle='None',
         color='navy', label=r'$/Delta$fSCA 75_99',linewidth=20, marker = "o", markersize = 60)#, markerfacecolor = 'black')
plt.plot(A0t5075[0,A0t5075[0,:].argsort()],A0t5075[1,A0t5075[0,:].argsort()],linestyle='None',
         color='navy', label=r'$/Delta$fSCA 50_75',linewidth=20, marker = "s", markersize = 60)#, markerfacecolor = 'black')
plt.plot(A0t2550[0,A0t2550[0,:].argsort()],A0t2550[1,A0t2550[0,:].argsort()],linestyle='None',
         color='navy', label=r'$/Delta$fSCA 25_50',linewidth=20, marker = "X", markersize = 60)#, markerfacecolor = 'black')

lowess = sm.nonparametric.lowess
Z2 = lowess(sns_0t['fsca'],sns_0t['temp'],frac=0.7,it=2)
plt.plot(Z2[:,0],Z2[:,1],'navy',lw=20);
#Z2 = lowess(sns_hlvd_shl['fsca'],sns_hlvd_shl['temp'],frac=0.7,it=2)
#plt.plot(Z2[:,0],Z2[:,1],'darkcyan',lw=20);
#
# confidence band exposed
n2=len(sns_0t['temp'])
linreg2=lin_reg(sns_0t['temp'],sns_0t['fsca'])
Yhat2=linreg2['Yhat'];
sse2=linreg2['sse']; 
MSE2=sse2/np.float(n2-2)
barX2=np.mean(sns_0t['temp'])
XminusbarX2=sns_0t['temp']-barX2

s_of_yh_hat2=np.sqrt(MSE2*(1.0/n2+(sns_0t['temp']-barX2)**2/sum(XminusbarX2**2)))
W2=np.sqrt(2.0*stats.f.ppf(0.95,2,n2-2))
cb_upper2=Yhat2+W2*s_of_yh_hat2
cb_lower2=Yhat2-W2*s_of_yh_hat2

X_ci2 = []
cb_lower_pl2 = []
cb_upper_pl2 = []
for idx2 in np.argsort(sns_0t['temp']):
    X_ci2.append(sns_0t['temp'][idx2])
    cb_lower_pl2.append(cb_lower2[idx2])
    cb_upper_pl2.append(cb_upper2[idx2])

plt.fill_between(X_ci2, cb_lower_pl2, cb_upper_pl2, facecolor = 'lavender', edgecolor = 'k')

plt.ylabel('fSCA_open - fSCA_underTree', fontsize=80) #(fSCA_op - fSCA_ut) in sheltered
plt.yticks(fontsize=80)
plt.xlabel('average Dec to flights date temperature (C)', fontsize=90) #(fSCA_op - fSCA_ut) in exposed
plt.xticks(fontsize=70)
plt.xlim((-7,3))

x = [-8,-5,-2,0,2,10]
y = [0,0,0,0,0,0]
plt.plot(x,y, color = 'black',linewidth=5)

plt.legend(fontsize=80, loc = 'lower left')
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_0t_point_25507599_bound.png')
#%%scatter plot
Bse0p7599 = np.array([tempDec7599,fsca_shEx_0p7599])
Bse0p5075 = np.array([tempDec5075,fsca_shEx_0p5075])
Bse0p2550 = np.array([tempDec2550,fsca_shEx_0p2550])
BseUt7599 = np.array([tempDec7599,fsca_shEx_ut7599])
BseUt5075 = np.array([tempDec5075,fsca_shEx_ut5075])
BseUt2550 = np.array([tempDec2550,fsca_shEx_ut2550])

fsca_shEx_0p = []
fsca_shEx_0p.extend(fsca_shEx_0p7599); fsca_shEx_0p.extend(fsca_shEx_0p5075); fsca_shEx_0p.extend(fsca_shEx_0p2550);
fsca_shEx_ut = []
fsca_shEx_ut.extend(fsca_shEx_ut7599); fsca_shEx_ut.extend(fsca_shEx_ut5075); fsca_shEx_ut.extend(fsca_shEx_ut2550);

sns_shEx_ut = pd.DataFrame(np.array([temp_dec,fsca_shEx_ut]).T,columns=['temp','fsca'])
sns_shEx_ut.dropna(axis=0, inplace = True)
sns_shEx_ut.index = np.arange(0,len(sns_shEx_ut))
sns_shEx_0p = pd.DataFrame(np.array([temp_dec,fsca_shEx_0p]).T,columns=['temp','fsca'])
sns_shEx_0p.dropna(axis=0, inplace = True)
sns_shEx_0p.index = np.arange(0,len(sns_shEx_0p))

plt.subplots(figsize=(60,40))
#plt.plot(Bse0p7599[0,Bse0p7599[0,:].argsort()],Bse0p7599[1,Bse0p7599[0,:].argsort()],linestyle='None',
#         color='mediumvioletred', label=r'$/Delta$fSCA open_75_99',linewidth=20, marker = "o", markersize = 70)#, markerfacecolor = 'black')
#plt.plot(Bse0p5075[0,Bse0p5075[0,:].argsort()],Bse0p5075[1,Bse0p5075[0,:].argsort()],linestyle='None',
#        color='mediumvioletred', label=r'$/Delta$fSCA open_50_75',linewidth=20, marker = "s", markersize = 70)#, markerfacecolor = 'black')
#plt.plot(Bse0p2550[0,Bse0p2550[0,:].argsort()],Bse0p2550[1,Bse0p2550[0,:].argsort()],linestyle='None',
#        color='mediumvioletred', label=r'$/Delta$fSCA open_25_50',linewidth=20, marker = "X", markersize = 70)#, markerfacecolor = 'black')
plt.plot(BseUt7599[0,BseUt7599[0,:].argsort()],BseUt7599[1,BseUt7599[0,:].argsort()],linestyle='None',
         color='deepskyblue', label=r'$/Delta$fSCA underTree_75_99',linewidth=20, marker = "o", markersize = 70)#, markerfacecolor = 'black')
plt.plot(BseUt5075[0,BseUt5075[0,:].argsort()],BseUt5075[1,BseUt5075[0,:].argsort()],linestyle='None',
        color='deepskyblue', label=r'$/Delta$fSCA underTree_50_75',linewidth=20, marker = "s", markersize = 70)#, markerfacecolor = 'black')
plt.plot(BseUt2550[0,BseUt2550[0,:].argsort()],BseUt2550[1,BseUt2550[0,:].argsort()],linestyle='None',
        color='deepskyblue', label=r'$/Delta$fSCA underTree_25_50',linewidth=20, marker = "X", markersize = 70)#, markerfacecolor = 'black')

lowess = sm.nonparametric.lowess
#Z3 = lowess(sns_shEx_0p['fsca'],sns_shEx_0p['temp'],frac=0.7,it=2)
#plt.plot(Z3[:,0],Z3[:,1],'mediumvioletred',lw=20);
Z4 = lowess(sns_shEx_ut['fsca'],sns_shEx_ut['temp'],frac=0.6,it=2)
plt.plot(Z4[:,0],Z4[:,1],'deepskyblue',lw=20);

## confidence band exposed
#n3=len(sns_shEx_0p['temp'])
#linreg3=lin_reg(sns_shEx_0p['temp'],sns_shEx_0p['fsca'])
#Yhat3=linreg3['Yhat'];
#sse3=linreg3['sse']; 
#MSE3=sse3/np.float(n3-2)
#barX3=np.mean(sns_shEx_0p['temp'])
#XminusbarX3=sns_shEx_0p['temp']-barX3
#
#s_of_yh_hat3=np.sqrt(MSE3*(1.0/n3+(sns_shEx_0p['temp']-barX3)**2/sum(XminusbarX3**2)))
#W3=np.sqrt(2.0*stats.f.ppf(0.95,2,n3-2))
#cb_upper3=Yhat3+W3*s_of_yh_hat3
#cb_lower3=Yhat3-W3*s_of_yh_hat3
#
#X_ci3 = []
#cb_lower_pl3 = []
#cb_upper_pl3 = []
#for idx3 in np.argsort(sns_shEx_0p['temp']):
#    X_ci3.append(sns_shEx_0p['temp'][idx3])
#    cb_lower_pl3.append(cb_lower3[idx3])
#    cb_upper_pl3.append(cb_upper3[idx3])
#
#plt.fill_between(X_ci3, cb_lower_pl3, cb_upper_pl3, facecolor = 'lavenderblush', edgecolor = 'k')

# confidence band exposed
n4=len(sns_shEx_ut['temp'])
linreg4=lin_reg(sns_shEx_ut['temp'],sns_shEx_ut['fsca'])
Yhat4=linreg4['Yhat'];
sse4=linreg4['sse']; 
MSE4=sse4/np.float(n4-2)
barX4=np.mean(sns_shEx_ut['temp'])
XminusbarX4=sns_shEx_ut['temp']-barX4

s_of_yh_hat4=np.sqrt(MSE4*(1.0/n4+(sns_shEx_ut['temp']-barX4)**2/sum(XminusbarX4**2)))
W4=np.sqrt(2.0*stats.f.ppf(0.95,2,n4-2))
cb_upper4=Yhat4+W4*s_of_yh_hat4
cb_lower4=Yhat4-W4*s_of_yh_hat4

X_ci4 = []
cb_lower_pl4 = []
cb_upper_pl4 = []
for idx4 in np.argsort(sns_shEx_ut['temp']):
    X_ci4.append(sns_shEx_ut['temp'][idx4])
    cb_lower_pl4.append(cb_lower4[idx4])
    cb_upper_pl4.append(cb_upper4[idx4])

plt.fill_between(X_ci4, cb_lower_pl4, cb_upper_pl4, facecolor = 'aliceblue', edgecolor = 'k')

plt.ylabel('fSCA_sheltered - fSCA_expoxed', fontsize=80) #(fSCA_op - fSCA_ut) in sheltered
plt.yticks(fontsize=80)
plt.xlabel('average Dec to flights date temperature (C)', fontsize=90) #(fSCA_op - fSCA_ut) in exposed
plt.xticks(fontsize=70)
plt.xlim((-7,3))

x = [-8,-5,-2,0,2,10]
y = [0,0,0,0,0,0]
plt.plot(x,y, color = 'black',linewidth=5)

plt.legend(fontsize=60, loc = 'upper left')
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_shEx_ut_point_25507599_bound.png')
#%%scatter plot
ChlE7599 = np.array([tempDec7599,fsca_hlvd_exp7599])
ChlE5075 = np.array([tempDec5075,fsca_hlvd_exp5075])
ChlE2550 = np.array([tempDec2550,fsca_hlvd_exp2550])
ChlS7599 = np.array([tempDec7599,fsca_hlvd_shl7599])
ChlS5075 = np.array([tempDec5075,fsca_hlvd_shl5075])
ChlS2550 = np.array([tempDec2550,fsca_hlvd_shl2550])

fsca_hlvd_shl = []
fsca_hlvd_shl.extend(fsca_hlvd_shl7599); fsca_hlvd_shl.extend(fsca_hlvd_shl5075); fsca_hlvd_shl.extend(fsca_hlvd_shl2550);
fsca_hlvd_exp = []
fsca_hlvd_exp.extend(fsca_hlvd_exp7599); fsca_hlvd_exp.extend(fsca_hlvd_exp5075); fsca_hlvd_exp.extend(fsca_hlvd_exp2550);

sns_hlvd_exp = pd.DataFrame(np.array([temp_dec,fsca_hlvd_exp]).T,columns=['temp','fsca'])
sns_hlvd_exp.dropna(axis=0, inplace = True)
sns_hlvd_exp.index = np.arange(0,len(sns_hlvd_exp))
sns_hlvd_shl = pd.DataFrame(np.array([temp_dec,fsca_hlvd_shl]).T,columns=['temp','fsca'])
sns_hlvd_shl.dropna(axis=0, inplace = True)
sns_hlvd_shl.index = np.arange(0,len(sns_hlvd_shl))

plt.subplots(figsize=(60,40))

plt.plot(ChlE7599[0,ChlE7599[0,:].argsort()],ChlE7599[1,ChlE7599[0,:].argsort()],linestyle='None',
         color='darkorange', label=r'$/Delta$fSCA exposed_75_99',linewidth=20, marker = "o", markersize = 70)#, markerfacecolor = 'black'
plt.plot(ChlE5075[0,ChlE5075[0,:].argsort()],ChlE5075[1,ChlE5075[0,:].argsort()],linestyle='None',
        color='darkorange', label=r'$/Delta$fSCA exposed_50_75',linewidth=20, marker = "s", markersize = 70)#, markerfacecolor = 'black')
plt.plot(ChlE2550[0,ChlE2550[0,:].argsort()],ChlE2550[1,ChlE2550[0,:].argsort()],linestyle='None',
        color='darkorange', label=r'$/Delta$fSCA exposed_25_50',linewidth=20, marker = "X", markersize = 70)#, markerfacecolor = 'black')
#plt.plot(ChlS7599[0,ChlS7599[0,:].argsort()],ChlS7599[1,ChlS7599[0,:].argsort()],linestyle='None',
#         color='darkgreen', label=r'$/Delta$fSCA sheltered_75_99',linewidth=20, marker = "o", markersize = 70)#, markerfacecolor = 'black')
#plt.plot(ChlS5075[0,ChlS5075[0,:].argsort()],ChlS5075[1,ChlS5075[0,:].argsort()],linestyle='None',
#        color='darkgreen', label=r'$/Delta$fSCA sheltered_50_75',linewidth=20, marker = "s", markersize = 70)#, markerfacecolor = 'black')
#plt.plot(ChlS2550[0,ChlS2550[0,:].argsort()],ChlS2550[1,ChlS2550[0,:].argsort()],linestyle='None',
#        color='darkgreen', label=r'$/Delta$fSCA sheltered_25_50',linewidth=20, marker = "X", markersize = 70)#, markerfacecolor = 'black')

#sns.regplot(x="temp", y="fsca", data = sns_data, lowess = True)#, fit_reg = True,  order=2, ci=90 
lowess = sm.nonparametric.lowess
Z = lowess(sns_hlvd_exp['fsca'],sns_hlvd_exp['temp'],frac=0.6,it=2)
plt.plot(Z[:,0],Z[:,1],'tomato',lw=20);
#Z2 = lowess(sns_hlvd_shl['fsca'],sns_hlvd_shl['temp'],frac=0.7,it=2)
#plt.plot(Z2[:,0],Z2[:,1],'darkcyan',lw=20);
#
# confidence band exposed
n=len(sns_hlvd_exp['temp'])
linreg=lin_reg(sns_hlvd_exp['temp'],sns_hlvd_exp['fsca'])
Yhat=linreg['Yhat'];
sse=linreg['sse']; 
MSE=sse/np.float(n-2)
barX=np.mean(sns_hlvd_exp['temp'])
XminusbarX=sns_hlvd_exp['temp']-barX

s_of_yh_hat=np.sqrt(MSE*(1.0/n+(sns_hlvd_exp['temp']-barX)**2/sum(XminusbarX**2)))
W=np.sqrt(2.0*stats.f.ppf(0.95,2,n-2))
cb_upper=Yhat+W*s_of_yh_hat
cb_lower=Yhat-W*s_of_yh_hat

X_ci = []
cb_lower_pl = []
cb_upper_pl = []
for idx in np.argsort(sns_hlvd_exp['temp']):
    X_ci.append(sns_hlvd_exp['temp'][idx])
    cb_lower_pl.append(cb_lower[idx])
    cb_upper_pl.append(cb_upper[idx])

plt.fill_between(X_ci, cb_lower_pl, cb_upper_pl, facecolor = 'papayawhip', edgecolor = 'k')

## confidence band sheltered
#n1=len(sns_hlvd_shl['temp'])
#linreg1=lin_reg(sns_hlvd_shl['temp'],sns_hlvd_shl['fsca'])
#Yhat1=linreg1['Yhat'];
#sse1=linreg1['sse']; 
#MSE1=sse1/np.float(n1-2)
#barX1=np.mean(sns_hlvd_shl['temp'])
#XminusbarX1=sns_hlvd_shl['temp']-barX1
#
#s_of_yh_hat1=np.sqrt(MSE1*(1.0/n1+(sns_hlvd_shl['temp']-barX1)**2/sum(XminusbarX1**2)))
#W1=np.sqrt(2.0*stats.f.ppf(0.95,2,n1-2))
#cb_upper1=Yhat1+W1*s_of_yh_hat1
#cb_lower1=Yhat1-W1*s_of_yh_hat1
#
#X_ci1 = []
#cb_lower_pl1 = []
#cb_upper_pl1 = []
#for idx1 in np.argsort(sns_hlvd_shl['temp']):
#    X_ci1.append(sns_hlvd_shl['temp'][idx1])
#    cb_lower_pl1.append(cb_lower1[idx1])
#    cb_upper_pl1.append(cb_upper1[idx1])
#
#plt.fill_between(X_ci1, cb_lower_pl1, cb_upper_pl1, facecolor = 'mintcream', edgecolor = 'k')
#plt.plot(X_ci1,cb_lower_pl1,'k--',linewidth=3)
#plt.plot(X_ci1,cb_upper_pl1,'k--',linewidth=3);

plt.ylabel('fSCA_highVD - fSCA_lowVD', fontsize=80) #(fSCA_op - fSCA_ut) in sheltered
plt.yticks(fontsize=80)
plt.xlabel('average Dec to flights date temperature (C)', fontsize=90) #(fSCA_op - fSCA_ut) in exposed
plt.xticks(fontsize=70)
plt.xlim((-7,3))

x = [-8,-5,-2,0,2,10]
y = [0,0,0,0,0,0]
plt.plot(x,y, color = 'black',linewidth=5)

plt.legend(fontsize=60, loc = 'upper left')
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_lhvd_exp_25507599_point_bound.png')

#%%
#plt.plot(X_ci,cb_lower_pl,'k--',linewidth=3)
#plt.plot(X_ci,cb_upper_pl,'k--',linewidth=3);
#plt.plot(Z[:,0],Z[:,1],'g-',lw=5);

#temp_exp = [-1.0237605523464839, -0.30603496093749977, 0.987795836425781, -5.616034765624999,
#            -6.475201755371095, -0.586671245346484, -0.0011849609374998682, 0.5393539062499997,
#            -3.456034765624999, 0.17301503906250004, 2.203845836425781, -3.9499517553710946]
#
#fsca_hlvd_exp = [-1.3878828787384805, -5.616993888029953, 1.4124858732115229, -1.4123652457854092,
#                 -6.164086339731853, -2.5201319732285725, 14.830789597936466, -1.4641616046594237,
#                 -1.0898873953720738, 8.997808044408217, 8.772116291637879, 0.03628949774000034]


#label4 = ['SCWC 26MAR2016-exposed-25-50','SCWC 17APR2016-exposed-25-50','SCWC 18MAY2016-exposed-25-50','KREW 2010-exposed-25-50','JRBN 2010-exposed-25-50','NRC 2010-exposed-25-50']
#label5 = ['SCWC 26MAR2016-sheltered-25-50-25-50','SCWC 17APR2016-sheltered-25-50','SCWC 18MAY2016-sheltered-25-50','KREW 2010-sheltered-25-50','JRBN 2010-sheltered-25-50','NRC 2010-sheltered-25-50']
#label6 = ['SCWC 26MAR2016-exposed-50-75','SCWC 17APR2016-exposed-50-75','SCWC 18MAY2016-exposed-50-75','KREW 2010-exposed-50-75','JRBN 2010-exposed-50-75','NRC 2010-exposed-50-75']
#label7 = ['SCWC 26MAR2016-sheltered-50-75','SCWC 17APR2016-sheltered-50-75','SCWC 18MAY2016-sheltered-50-75','KREW 2010-sheltered-50-75','JRBN 2010-sheltered-50-75','NRC 2010-sheltered-50-75']
#
##color = ['plum','purple','hotpink','red','darkgreen','deepskyblue']#'olive',
#marker = ['^','>','v','D','s','o']
#
#plt.subplots(figsize=(60,60))
#
#for i in range (1,6):
#    plt.scatter(fsca_hlvd_exp2550[i],fsca_shEx_ut2550[i], s=80**2, color = 'lightblue', marker = marker[i], label =  label4[i])#
#    plt.scatter(fsca_hlvd_shl2550[i],fsca_shEx_ut2550[i], s=80**2, color = 'blue', marker = marker[i], label =  label5[i])#
#    plt.scatter(fsca_hlvd_exp5075[i],fsca_shEx_ut5075[i], s=80**2, color = 'violet', marker = marker[i], label =  label6[i])#
#    plt.scatter(fsca_hlvd_shl5075[i],fsca_shEx_ut5075[i], s=80**2, color = 'navy', marker = marker[i], label =  label7[i])#
#
#plt.ylabel('(fSCA_sheltered - fSCA_exposed)', fontsize=80)
#plt.yticks(fontsize=80)
#plt.xlabel('(fSCA_highVD - fSCA_lowVD)', fontsize=80)
#plt.xticks(fontsize=80)
#plt.ylim((-40,40))
#plt.xlim((-41,40))
#plt.legend(fontsize=50, loc = 'upper left')
#
#x = [-41,-8,-5,-2,0,2,40]
#y = [0,0,0,0,0,0,0]
#plt.plot(x,y, color = 'black') #line_temp_ls[2*i],
#
#y2 = [-40,-8,-5,-2,0,2,40]
#x2 = [0,0,0,0,0,0,0]
#plt.plot(x2,y2, color = 'black') #line_temp_ls[2*i],   
##plt.title('50% to 75% fSCA', fontsize=60) #, y=2.24, x=-0.65 loc = 'right', 
# 
#plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_fusion_lidar/deltaFsca_expShl_hlvd_melt.png')

#%%
#fSCA_0p_shl_2550 = mean_fsca_Limit (fSCA_0p_shl_sc26m,fSCA_0p_shl_sc17a,fSCA_0p_shl_sc18m,fSCA_0p_shl_Krew,fSCA_0p_shl_Jmz,fSCA_0p_shl_nrc,25,50) 
#fSCA_0p_exp_2550 = mean_fsca_Limit (fSCA_0p_exp_sc26m,fSCA_0p_exp_sc17a,fSCA_0p_exp_sc18m,fSCA_0p_exp_Krew,fSCA_0p_exp_Jmz,fSCA_0p_exp_nrc,25,50) 
#fSCA_ut_shl_2550 = mean_fsca_Limit (fSCA_ut_shl_sc26m,fSCA_ut_shl_sc17a,fSCA_ut_shl_sc18m,fSCA_ut_shl_Krew,fSCA_ut_shl_Jmz,fSCA_ut_shl_nrc,25,50) 
#fSCA_ut_exp_2550 = mean_fsca_Limit (fSCA_ut_exp_sc26m,fSCA_ut_exp_sc17a,fSCA_ut_exp_sc18m,fSCA_ut_exp_Krew,fSCA_ut_exp_Jmz,fSCA_ut_exp_nrc,25,50) 
#
#fSCA_0p_shl_5075 = mean_fsca_Limit (fSCA_0p_shl_sc26m,fSCA_0p_shl_sc17a,fSCA_0p_shl_sc18m,fSCA_0p_shl_Krew,fSCA_0p_shl_Jmz,fSCA_0p_shl_nrc,50,75) 
#fSCA_0p_exp_5075 = mean_fsca_Limit (fSCA_0p_exp_sc26m,fSCA_0p_exp_sc17a,fSCA_0p_exp_sc18m,fSCA_0p_exp_Krew,fSCA_0p_exp_Jmz,fSCA_0p_exp_nrc,50,75) 
#fSCA_ut_shl_5075 = mean_fsca_Limit (fSCA_ut_shl_sc26m,fSCA_ut_shl_sc17a,fSCA_ut_shl_sc18m,fSCA_ut_shl_Krew,fSCA_ut_shl_Jmz,fSCA_ut_shl_nrc,50,75) 
#fSCA_ut_exp_5075 = mean_fsca_Limit (fSCA_ut_exp_sc26m,fSCA_ut_exp_sc17a,fSCA_ut_exp_sc18m,fSCA_ut_exp_Krew,fSCA_ut_exp_Jmz,fSCA_ut_exp_nrc,50,75) 
#
#for sites in range (6):
#    plt.scatter(ind[sites] - 3*width/2, fSCA_0p_shl_2550[sites], #width, yerr=men_std,
#                color=colorp[sites], s=60**2, marker = markerp[0]) #, label='0p_shl_2550' palevioletred mediumvioletred firebrick
#    plt.scatter(ind[sites] - width/2, fSCA_0p_exp_2550[sites], #width, yerr=women_std,
#                color=colorp[sites], s=60**2, marker = markerp[1]) #, label='0p_exp_2550'palegoldenrod  goldenrod darkgoldenrod
#    plt.scatter(ind[sites] + width/2, fSCA_ut_shl_2550[sites], #width, yerr=women_std,
#                color=colorp[sites], s=60**2, marker = markerp[2]) #, label='ut_shl_2550'seagreen palegreen lawngreen darkolivegreen
#    plt.scatter(ind[sites] + 3*width/2, fSCA_ut_exp_2550[sites], #width, yerr=women_std,
#                color=colorp[sites], s=60**2, marker = markerp[3]) #, label='ut_exp_2550'powderblue dodgerblue cornflowerblue cadetblue blueviolet
#    plt.scatter(ind[sites] - 3*width/2, fSCA_0p_shl_5075[sites], #width, 
#                color=colorp[sites], s=85**2, marker = markerp[0])#, label='0p_shl_5075'bottom=fSCA_0p_shl_2550, 
#    plt.scatter(ind[sites] - width/2, fSCA_0p_exp_5075[sites], #width, yerr=women_std,
#                color=colorp[sites], s=85**2, marker = markerp[1]) #, label='0p_exp_5075'palegoldenrod  goldenrod darkgoldenrod
#    plt.scatter(ind[sites] + width/2, fSCA_ut_shl_5075[sites], #width, yerr=women_std,
#                color=colorp[sites], s=85**2, marker = markerp[2]) #, label='ut_shl_5075'seagreen palegreen lawngreen darkolivegreen
#    plt.scatter(ind[sites] + 3*width/2, fSCA_ut_exp_5075[sites], #width, yerr=women_std,
#                color=colorp[sites], s=85**2, marker = markerp[3]) #, label='ut_exp_5075'powderblue dodgerblue cornflowerblue cadetblue blueviolet


#marker = ['^','>','v','D','s','o']
#for i in range (0,6):
#        plt.scatter(fsca_shEx_0p2550[i],fsca_shEx_ut2550[i], s=80**2, color = 'goldenrod', marker = marker[i], label =  label0[i])#
#        plt.scatter(fsca_0pUt_shl2550[i],fsca_0pUt_exp2550[i], s=80**2, color = 'darkorange', marker = marker[i], label =  label0[i])#

#    plt.scatter(fsca_shEx_0p2550[i],fsca_0t2550[i], s=80**2, color = 'goldenrod', marker = marker[i], label =  label0[i])#
#    plt.scatter(fsca_shEx_ut2550[i],fsca_0t2550[i], s=80**2, color = 'darkorange', marker = marker[i], label =  label1[i])#
#    plt.scatter(fsca_shEx_0p5070[i],fsca_0t5070[i], s=80**2, color = 'red', marker = marker[i], label =  label0[i])#
#    plt.scatter(fsca_shEx_ut5070[i],fsca_0t5070[i], s=80**2, color = 'darkred', marker = marker[i], label =  label1[i])#

#x = [-41,-8,-5,-2,0,2,40]
#y = [0,0,0,0,0,0,0]
#plt.plot(x,y, color = 'black') #line_temp_ls[2*i],

#y2 = [-40,-8,-5,-2,0,2,40]
#x2 = [0,0,0,0,0,0,0]
#plt.plot(x2,y2, color = 'black') #line_temp_ls[2*i],   


#%%# mixed effect random forest

# mixed effect random forest
#inin_0p_merf_krew = pd.concat([in_0p_rf_krew,in_0p_rf_krew['temp']], axis = 1)
#in_0p_merf_krew_ls = inin_0p_merf_krew.values
#class_len = (max(in_0p_merf_krew_ls[:,2]) - min(in_0p_merf_krew_ls[:,2]))/5
#in_0p_merf_krew_ls[:,5][(in_0p_merf_krew_ls[:,2]>=min(in_0p_merf_krew_ls[:,2]))&
#                  (in_0p_merf_krew_ls[:,2]<=min(in_0p_merf_krew_ls[:,2])+class_len)] = 11
#in_0p_merf_krew_ls[:,5][(in_0p_merf_krew_ls[:,2]>=min(in_0p_merf_krew_ls[:,2])+class_len)&
#                  (in_0p_merf_krew_ls[:,2]<=min(in_0p_merf_krew_ls[:,2])+class_len+class_len)] = 12
#in_0p_merf_krew_ls[:,5][(in_0p_merf_krew_ls[:,2]>=min(in_0p_merf_krew_ls[:,2])+class_len+class_len)&
#                  (in_0p_merf_krew_ls[:,2]<=min(in_0p_merf_krew_ls[:,2])+class_len+class_len+class_len)] = 13                  
#in_0p_merf_krew_ls[:,5][(in_0p_merf_krew_ls[:,2]>=min(in_0p_merf_krew_ls[:,2])+class_len+class_len+class_len)&
#                  (in_0p_merf_krew_ls[:,2]<=min(in_0p_merf_krew_ls[:,2])+class_len+class_len+class_len+class_len)] = 14         
#in_0p_merf_krew_ls[:,5][(in_0p_merf_krew_ls[:,2]>=min(in_0p_merf_krew_ls[:,2])+class_len+class_len+class_len+class_len)] = 15         
#
##features, clusters and Z
#in_0p_merf_krew = pd.DataFrame(in_0p_merf_krew_ls[:,0:5], columns = in_0p_rf_krew.columns)
#z_merf_krew = pd.DataFrame(np.ones((1355,1)), columns = ['z']) #
#cluster_0p_merf_krew = pd.DataFrame(in_0p_merf_krew_ls[:,5], columns = ['cluster'])
#
#in_0p_rf_krew_train1, in_0p_rf_krew_test1, out_0p_rf_krew_train1, out_0p_rf_krew_test1 = train_test_split(in_0p_rf_krew, out_0p_rf_krew, test_size=0.3, random_state=random_state, shuffle= True) # 70% training and 30% test
#
#from merf.utils import MERFDataGenerator
#from merf.merf import MERF
#mrf = MERF(n_estimators=300, max_iterations=100)
#mrf.fit(in_0p_merf_krew, z_merf_krew, cluster_0p_merf_krew.iloc[:,0], out_0p_rf_krew.iloc[:])
#
#mrf.score(in_0p_merf_krew, z_merf_krew, cluster_0p_merf_krew.iloc[:,0], out_0p_rf_krew.iloc[:])
#













