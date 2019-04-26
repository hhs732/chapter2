import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import itertools
import csv
from osgeo import gdal, ogr, osr
from mpl_toolkits.mplot3d import Axes3D
import os
import rasterio
from scipy import stats
#%% functions
# classification based on lat and long in 100x100 pixels
#veg_snow_temp_sc_df_tst = pd.DataFrame(veg_snow_temp_sc[0:1000,:], columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])
#lat_long_df_sc_tst = lat_long_df_sc[['x','y']][0:200]
#cluster100_sc = []
#for row in veg_snow_temp_sc_df.itertuples():
#    indexFinl = np.where(((lat_long_df_sc['x']-row.x<=100) & (lat_long_df_sc['y']-row.y<=100) & (lat_long_df_sc['x']-row.x>=0) & (lat_long_df_sc['y']-row.y>=0)))[0]#
#    if (np.size(indexFinl))>0:
#        cluster100_sc.append([lat_long_df_sc['x'][indexFinl[0]],lat_long_df_sc['y'][indexFinl[0]]])

#def classificationGridCells(points_df_ls,gridcell):#0.25
# 
#    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[0][0][0])&(points_df_ls[:,0]<gridcell[0][2][0])&
#                      (points_df_ls[:,1]>=gridcell[0][0][1])&(points_df_ls[:,1]<gridcell[0][1][1])]=11
#                
#    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[1][0][0])&(points_df_ls[:,0]<gridcell[1][2][0])&
#                      (points_df_ls[:,1]>=gridcell[1][0][1])&(points_df_ls[:,1]<gridcell[1][1][1])]=12
#                
#    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[2][0][0])&(points_df_ls[:,0]<gridcell[2][2][0])&
#                      (points_df_ls[:,1]>=gridcell[2][0][1])&(points_df_ls[:,1]<gridcell[2][1][1])]=13
#                
#    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[3][0][0])&(points_df_ls[:,0]<gridcell[3][2][0])&
#                      (points_df_ls[:,1]>=gridcell[3][0][1])&(points_df_ls[:,1]<gridcell[3][1][1])]=21
#                
#    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[4][0][0])&(points_df_ls[:,0]<gridcell[4][2][0])&
#                      (points_df_ls[:,1]>=gridcell[4][0][1])&(points_df_ls[:,1]<gridcell[4][1][1])]=22
#                
#    points_df_ls[:,6][(points_df_ls[:,0]>=gridcell[5][0][0])&(points_df_ls[:,0]<gridcell[5][2][0])&
#                      (points_df_ls[:,1]>=gridcell[5][0][1])&(points_df_ls[:,1]<gridcell[5][1][1])]=23
#    
#    pointFile_df = pd.DataFrame(points_df_ls,columns=['x','y','z','north','slope','vegClass','gridClass'])
#    
#    return pointFile_df.values#, pointFile_df

def classificationElevation(points_df_ls,elevationBand):#0.25
 
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[0])&(points_df_ls[:,2]<elevationBand[1])]=101
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[1])&(points_df_ls[:,2]<elevationBand[2])]=102
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[2])&(points_df_ls[:,2]<elevationBand[3])]=103
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[3])&(points_df_ls[:,2]<elevationBand[4])]=104
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[4])&(points_df_ls[:,2]<elevationBand[5])]=105
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[5])&(points_df_ls[:,2]<elevationBand[6])]=106
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[6])&(points_df_ls[:,2]<elevationBand[7])]=107
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[7])&(points_df_ls[:,2]<elevationBand[8])]=108
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[8])&(points_df_ls[:,2]<elevationBand[9])]=109
    points_df_ls[:,8][(points_df_ls[:,2]>=elevationBand[9])&(points_df_ls[:,2]<elevationBand[10])]=110
    
    points_df_ls[:,8][(points_df_ls[:,8]<=100)]=110
                
    return points_df_ls#, pointFile_df

def fSCAcalculation_elevNorthClassific(veg_snow_temp_elev_clsf_file): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification
    
    #under tree with snow, exposed
    ascii_ut_snow1011 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1021 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1031 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1041 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1051 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1061 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1071 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1081 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1091 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    
    #under tree no snow, exposed
    ascii_ut_nsnow1011 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1021 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1031 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1041 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1051 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1061 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1071 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1081 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1091 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_nsnow1101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])

    #under tree with snow, sheltered
    ascii_ut_snow1012 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1022 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1032 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1042 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1052 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1062 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1072 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1082 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1092 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_ut_snow1102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==77)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    
    #under tree no snow, sheltered
    ascii_ut_nsnow1012 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1022 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1032 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1042 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1052 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1062 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1072 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1082 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1092 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_ut_nsnow1102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-7)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])

    #fSCA under tree, exposed
    ascii_fSCA_ut1011 = 100.*ascii_ut_snow1011/(ascii_ut_snow1011+ascii_ut_nsnow1011)
    ascii_fSCA_ut1021 = 100.*ascii_ut_snow1021/(ascii_ut_snow1021+ascii_ut_nsnow1021)
    ascii_fSCA_ut1031 = 100.*ascii_ut_snow1031/(ascii_ut_snow1031+ascii_ut_nsnow1031)
    ascii_fSCA_ut1041 = 100.*ascii_ut_snow1041/(ascii_ut_snow1041+ascii_ut_nsnow1041)
    ascii_fSCA_ut1051 = 100.*ascii_ut_snow1051/(ascii_ut_snow1051+ascii_ut_nsnow1051)
    ascii_fSCA_ut1061 = 100.*ascii_ut_snow1061/(ascii_ut_snow1061+ascii_ut_nsnow1061)
    ascii_fSCA_ut1071 = 100.*ascii_ut_snow1071/(ascii_ut_snow1071+ascii_ut_nsnow1071)
    ascii_fSCA_ut1081 = 100.*ascii_ut_snow1081/(ascii_ut_snow1081+ascii_ut_nsnow1081)
    ascii_fSCA_ut1091 = 100.*ascii_ut_snow1091/(ascii_ut_snow1091+ascii_ut_nsnow1091)
    ascii_fSCA_ut1101 = 100.*ascii_ut_snow1101/(ascii_ut_snow1101+ascii_ut_nsnow1101)

    #fSCA under tree, exposed
    ascii_fSCA_ut1012 = 100.*ascii_ut_snow1012/(ascii_ut_snow1012+ascii_ut_nsnow1012)
    ascii_fSCA_ut1022 = 100.*ascii_ut_snow1022/(ascii_ut_snow1022+ascii_ut_nsnow1022)
    ascii_fSCA_ut1032 = 100.*ascii_ut_snow1032/(ascii_ut_snow1032+ascii_ut_nsnow1032)
    ascii_fSCA_ut1042 = 100.*ascii_ut_snow1042/(ascii_ut_snow1042+ascii_ut_nsnow1042)
    ascii_fSCA_ut1052 = 100.*ascii_ut_snow1052/(ascii_ut_snow1052+ascii_ut_nsnow1052)
    ascii_fSCA_ut1062 = 100.*ascii_ut_snow1062/(ascii_ut_snow1062+ascii_ut_nsnow1062)
    ascii_fSCA_ut1072 = 100.*ascii_ut_snow1072/(ascii_ut_snow1072+ascii_ut_nsnow1072)
    ascii_fSCA_ut1082 = 100.*ascii_ut_snow1082/(ascii_ut_snow1082+ascii_ut_nsnow1082)
    ascii_fSCA_ut1092 = 100.*ascii_ut_snow1092/(ascii_ut_snow1092+ascii_ut_nsnow1092)
    ascii_fSCA_ut1102 = 100.*ascii_ut_snow1102/(ascii_ut_snow1102+ascii_ut_nsnow1102)
    
    fSCA_unTr_exp = [ascii_fSCA_ut1011, ascii_fSCA_ut1021, ascii_fSCA_ut1031, ascii_fSCA_ut1041, 
                     ascii_fSCA_ut1051, ascii_fSCA_ut1061, ascii_fSCA_ut1071, ascii_fSCA_ut1081, 
                     ascii_fSCA_ut1091, ascii_fSCA_ut1101]
    fSCA_unTr_shl = [ascii_fSCA_ut1012, ascii_fSCA_ut1022, ascii_fSCA_ut1032, ascii_fSCA_ut1042, 
                     ascii_fSCA_ut1052, ascii_fSCA_ut1062, ascii_fSCA_ut1072, ascii_fSCA_ut1082, 
                     ascii_fSCA_ut1092, ascii_fSCA_ut1102]
    
    #0pen with snow, exposed
    ascii_0p_snow1011 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1021 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1031 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1041 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1051 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1061 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1071 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1081 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1091 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    
    #0pen no snow, exposed
    ascii_0p_nsnow1011 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1021 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1031 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1041 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1051 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1061 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1071 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1081 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1091 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_nsnow1101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])

    #0pen with snow, sheltered
    ascii_0p_snow1012 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1022 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1032 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1042 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1052 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1062 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1072 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1082 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1092 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    ascii_0p_snow1102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==44)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]<-0.1)])
    
    #0pen no snow, sheltered
    ascii_0p_nsnow1012 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==101)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1022 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==102)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1032 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==103)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1042 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==104)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1052 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==105)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1062 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==106)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1072 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==107)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1082 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==108)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1092 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==109)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])
    ascii_0p_nsnow1102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,5]==-4)&(veg_snow_temp_elev_clsf_file[:,8]==110)&(veg_snow_temp_elev_clsf_file[:,3]>0.1)])

    #fSCA 0pen, exposed
    ascii_fSCA_0p1011 = 100.*ascii_0p_snow1011/(ascii_0p_snow1011+ascii_0p_nsnow1011)
    ascii_fSCA_0p1021 = 100.*ascii_0p_snow1021/(ascii_0p_snow1021+ascii_0p_nsnow1021)
    ascii_fSCA_0p1031 = 100.*ascii_0p_snow1031/(ascii_0p_snow1031+ascii_0p_nsnow1031)
    ascii_fSCA_0p1041 = 100.*ascii_0p_snow1041/(ascii_0p_snow1041+ascii_0p_nsnow1041)
    ascii_fSCA_0p1051 = 100.*ascii_0p_snow1051/(ascii_0p_snow1051+ascii_0p_nsnow1051)
    ascii_fSCA_0p1061 = 100.*ascii_0p_snow1061/(ascii_0p_snow1061+ascii_0p_nsnow1061)
    ascii_fSCA_0p1071 = 100.*ascii_0p_snow1071/(ascii_0p_snow1071+ascii_0p_nsnow1071)
    ascii_fSCA_0p1081 = 100.*ascii_0p_snow1081/(ascii_0p_snow1081+ascii_0p_nsnow1081)
    ascii_fSCA_0p1091 = 100.*ascii_0p_snow1091/(ascii_0p_snow1091+ascii_0p_nsnow1091)
    ascii_fSCA_0p1101 = 100.*ascii_0p_snow1101/(ascii_0p_snow1101+ascii_0p_nsnow1101)

    #fSCA 0pen, exposed
    ascii_fSCA_0p1012 = 100.*ascii_0p_snow1012/(ascii_0p_snow1012+ascii_0p_nsnow1012)
    ascii_fSCA_0p1022 = 100.*ascii_0p_snow1022/(ascii_0p_snow1022+ascii_0p_nsnow1022)
    ascii_fSCA_0p1032 = 100.*ascii_0p_snow1032/(ascii_0p_snow1032+ascii_0p_nsnow1032)
    ascii_fSCA_0p1042 = 100.*ascii_0p_snow1042/(ascii_0p_snow1042+ascii_0p_nsnow1042)
    ascii_fSCA_0p1052 = 100.*ascii_0p_snow1052/(ascii_0p_snow1052+ascii_0p_nsnow1052)
    ascii_fSCA_0p1062 = 100.*ascii_0p_snow1062/(ascii_0p_snow1062+ascii_0p_nsnow1062)
    ascii_fSCA_0p1072 = 100.*ascii_0p_snow1072/(ascii_0p_snow1072+ascii_0p_nsnow1072)
    ascii_fSCA_0p1082 = 100.*ascii_0p_snow1082/(ascii_0p_snow1082+ascii_0p_nsnow1082)
    ascii_fSCA_0p1092 = 100.*ascii_0p_snow1092/(ascii_0p_snow1092+ascii_0p_nsnow1092)
    ascii_fSCA_0p1102 = 100.*ascii_0p_snow1102/(ascii_0p_snow1102+ascii_0p_nsnow1102)
    
    fSCA_0pen_exp = [ascii_fSCA_0p1011, ascii_fSCA_0p1021, ascii_fSCA_0p1031, ascii_fSCA_0p1041, 
                     ascii_fSCA_0p1051, ascii_fSCA_0p1061, ascii_fSCA_0p1071, ascii_fSCA_0p1081, 
                     ascii_fSCA_0p1091, ascii_fSCA_0p1101]
    fSCA_0pen_shl = [ascii_fSCA_0p1012, ascii_fSCA_0p1022, ascii_fSCA_0p1032, ascii_fSCA_0p1042, 
                     ascii_fSCA_0p1052, ascii_fSCA_0p1062, ascii_fSCA_0p1072, ascii_fSCA_0p1082, 
                     ascii_fSCA_0p1092, ascii_fSCA_0p1102]
    
    return fSCA_unTr_exp, fSCA_unTr_shl, fSCA_0pen_exp, fSCA_0pen_shl

def fsca0pn_fscaUnTr(fSCA_ut_exp,fSCA_ut_shl,fSCA_0p_exp,fSCA_0p_shl):
    #fSCA in exposed area
    fSCA_0u01 = (fSCA_0p_exp[0] - fSCA_ut_exp[0])/fSCA_0p_exp[0]
    fSCA_0u11 = (fSCA_0p_exp[1] - fSCA_ut_exp[1])/fSCA_0p_exp[1]
    fSCA_0u21 = (fSCA_0p_exp[2] - fSCA_ut_exp[2])/fSCA_0p_exp[2]
    fSCA_0u31 = (fSCA_0p_exp[3] - fSCA_ut_exp[3])/fSCA_0p_exp[3]
    fSCA_0u41 = (fSCA_0p_exp[4] - fSCA_ut_exp[4])/fSCA_0p_exp[4]
    fSCA_0u51 = (fSCA_0p_exp[5] - fSCA_ut_exp[5])/fSCA_0p_exp[5]
    fSCA_0u61 = (fSCA_0p_exp[6] - fSCA_ut_exp[6])/fSCA_0p_exp[6]
    fSCA_0u71 = (fSCA_0p_exp[7] - fSCA_ut_exp[7])/fSCA_0p_exp[7]
    fSCA_0u81 = (fSCA_0p_exp[8] - fSCA_ut_exp[8])/fSCA_0p_exp[8]
    fSCA_0u91 = (fSCA_0p_exp[9] - fSCA_ut_exp[9])/fSCA_0p_exp[9]
    
    #fSCA in shelter area
    fSCA_0u02 = (fSCA_0p_shl[0] - fSCA_ut_shl[0])/fSCA_0p_shl[0]
    fSCA_0u12 = (fSCA_0p_shl[1] - fSCA_ut_shl[1])/fSCA_0p_shl[1]
    fSCA_0u22 = (fSCA_0p_shl[2] - fSCA_ut_shl[2])/fSCA_0p_shl[2]
    fSCA_0u32 = (fSCA_0p_shl[3] - fSCA_ut_shl[3])/fSCA_0p_shl[3]
    fSCA_0u42 = (fSCA_0p_shl[4] - fSCA_ut_shl[4])/fSCA_0p_shl[4]
    fSCA_0u52 = (fSCA_0p_shl[5] - fSCA_ut_shl[5])/fSCA_0p_shl[5]
    fSCA_0u62 = (fSCA_0p_shl[6] - fSCA_ut_shl[6])/fSCA_0p_shl[6]
    fSCA_0u72 = (fSCA_0p_shl[7] - fSCA_ut_shl[7])/fSCA_0p_shl[7]
    fSCA_0u82 = (fSCA_0p_shl[8] - fSCA_ut_shl[8])/fSCA_0p_shl[8]
    fSCA_0u92 = (fSCA_0p_shl[9] - fSCA_ut_shl[9])/fSCA_0p_shl[9]    
    
    fSCA_0p_ut_exp = [fSCA_0u01, fSCA_0u11, fSCA_0u21, fSCA_0u31, fSCA_0u41, fSCA_0u51, fSCA_0u61, 
                      fSCA_0u71, fSCA_0u81, fSCA_0u91] 
    fSCA_0p_ut_shl = [fSCA_0u02, fSCA_0u12, fSCA_0u22, fSCA_0u32, fSCA_0u42, fSCA_0u52, fSCA_0u62, 
                      fSCA_0u72, fSCA_0u82, fSCA_0u92]
    
    return fSCA_0p_ut_exp, fSCA_0p_ut_shl

def calculationMeanTempforEachElevClassification(veg_snow_temp_elev_clsf_file):
    temp101_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==101)][:,7])
    temp102_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==102)][:,7])
    temp103_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==103)][:,7])
    temp104_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==104)][:,7])
    temp105_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==105)][:,7])
    temp106_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==106)][:,7])
    temp107_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==107)][:,7])
    temp108_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==108)][:,7])
    temp109_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==109)][:,7])
    temp110_mean = np.mean(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==110)][:,7])
    
    meanTemp = [temp101_mean, temp102_mean, temp103_mean, temp104_mean, temp105_mean, temp106_mean, temp107_mean, temp108_mean, temp109_mean, temp110_mean]
    
    return meanTemp

def number0fGridsInEachElevationBand(veg_snow_temp_elev_clsf_file): #veg_snow_temp_elev_clsf_file[:,8]==101 contains elevation classification
    
    numGrid101 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==101)])
    numGrid102 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==102)])
    numGrid103 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==103)])
    numGrid104 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==104)])
    numGrid105 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==105)])
    numGrid106 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==106)])
    numGrid107 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==107)])
    numGrid108 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==108)])
    numGrid109 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==109)])
    numGrid110 = len(veg_snow_temp_elev_clsf_file[(veg_snow_temp_elev_clsf_file[:,8]==110)])

    numGrids = [numGrid101, numGrid102, numGrid103, numGrid104, 
                 numGrid105, numGrid106, numGrid107, numGrid108, 
                 numGrid109, numGrid110]
    
    return numGrids
#%% load data sagehen creek 26 march 2010
ascii_grid_veg_indexKrew = np.load('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/ascii_index_pr_krew.npy')
#calculating topo_dimension from cut northness
lat_long_nad83_krew = [np.min(ascii_grid_veg_indexKrew[:,0]),np.max(ascii_grid_veg_indexKrew[:,0]),np.min(ascii_grid_veg_indexKrew[:,1]),np.max(ascii_grid_veg_indexKrew[:,1])]
lat_long_degre_krew = [37.029257, 37.08022, -119.20915406614355, -119.18251130946379]

# slope less than 30
ascii_grid_veg_index_lsKrew = ascii_grid_veg_indexKrew[(ascii_grid_veg_indexKrew[:,5]>-90)&(ascii_grid_veg_indexKrew[:,4]<30)]
aaaa_test = ascii_grid_veg_index_lsKrew[700780:710057,:]

elev_krew = ascii_grid_veg_index_lsKrew[:,2]
tempDJF_gm = -0.0033*elev_krew + 7.2339 #y = -0.0033x + 7.2339
 #R² = 0.5074

veg_snow_temp_Krew = np.vstack([ascii_grid_veg_index_lsKrew[:,0],ascii_grid_veg_index_lsKrew[:,1].astype(int),ascii_grid_veg_index_lsKrew[:,2],
                                ascii_grid_veg_index_lsKrew[:,3],ascii_grid_veg_index_lsKrew[:,4],ascii_grid_veg_index_lsKrew[:,5],
                                ascii_grid_veg_index_lsKrew[:,6],tempDJF_gm,tempDJF_gm]).T

veg_snow_temp_dfKrew = pd.DataFrame(veg_snow_temp_Krew, columns=['x','y','z','nrth','slp','snwIndx','grid','temp','indx'])

elevationBand = np.arange(min(veg_snow_temp_dfKrew['z']),max(veg_snow_temp_dfKrew['z']),67)
veg_snow_temp_clsfKrew = classificationElevation(veg_snow_temp_Krew,elevationBand)

numGrids_krew = number0fGridsInEachElevationBand(veg_snow_temp_clsfKrew)

meanTemp_elevClass = calculationMeanTempforEachElevClassification(veg_snow_temp_clsfKrew)   
# under canopy and 0pen fsca
fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew = fSCAcalculation_elevNorthClassific(veg_snow_temp_clsfKrew)

#fsca open -minus- fsca under canopy
fSCA_0u_exp_Krew, fSCA_0u_shl_Krew = fsca0pn_fscaUnTr(fSCA_ut_exp_Krew, fSCA_ut_shl_Krew, fSCA_0p_exp_Krew, fSCA_0p_shl_Krew)

#%% ploting DJF temp
fig, ax = plt.subplots(figsize=(20,15))
ax.scatter(meanTemp_elevClass,fSCA_0u_exp_Krew, label = 'Krew 2010 --- exposed', s=20**2, color = 'orange')#r2=%.6f'%(r_value_tempFscaKrew**2)
ax.scatter(meanTemp_elevClass,fSCA_0u_shl_Krew, label = 'Krew 2010 --- sheltered', s=20**2, color = 'red')


plt.ylabel('delta fSCA/fSCA_op', fontsize=35)
plt.yticks(fontsize=30)
plt.xlabel('average temp of DJF (c)', fontsize=35)
plt.xticks(fontsize=30)

plt.legend(fontsize=20, loc = 'upper right')

# the line equation:
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/LiDarAnalysis/lapseRate_gridMet/tempLsr_krew3.png')











