###       /bin/bash runTestCases_dockerSC.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
from sklearn.metrics import mean_squared_error
import itertools
import csv

def readSpecificDatafromAllHRUs(variablename,hruname,day):
    dayData = []
    for names in hruname:
        dayData.append(variablename[names][day])
    return dayData

def readVariablefromNcfilesDatasetasDF(NcfilesDataset,variable,hruname):
    variableNameList = []
    for datasets in NcfilesDataset:
        variableNameList.append(pd.DataFrame(datasets[variable][:][:]))
    variableNameDF = pd.concat (variableNameList, axis=1)
    variableNameDF.columns = hruname
    counter = pd.DataFrame(np.arange(0,np.size(variableNameDF[hruname[0]])),columns=['counter'])
    counter.set_index(variableNameDF.index,inplace=True)
    variableNameDF = pd.concat([counter, variableNameDF], axis=1)
    return variableNameDF

def calculateDay0fSnowDissappearance(swe_df,hruname_df):
    av_swe_df4000 = swe_df[:][4000:8784]
    av_swe_df13000 = swe_df[:][13000:17137]
    
    zerosnowdate = []
    for val in hruname_df[0]:
        zerosnowdate.append(np.where(av_swe_df4000[val]==0))
        zerosnowdate.append(np.where(av_swe_df13000[val]==0))

    zerosnowdate_omg = [item[0] for item in zerosnowdate] #change tuple to array

    for zdindx in range (len(zerosnowdate_omg)/2):
        for i,item in enumerate(zerosnowdate_omg[2*zdindx]):
            if np.size(item) == 0:
                zerosnowdate_omg[2*zdindx][i] = 4783
        for i,item in enumerate(zerosnowdate_omg[2*zdindx]):
            zerosnowdate_omg[2*zdindx][i] = zerosnowdate_omg[2*zdindx][i]+4000

        for i,item in enumerate(zerosnowdate_omg[2*zdindx+1]):
            if np.size(item) == 0:
                zerosnowdate_omg[2*zdindx+1][2*zdindx+1] = 13137
        for i,item in enumerate(zerosnowdate_omg[2*zdindx+1]):
            zerosnowdate_omg[2*zdindx+1][i] = zerosnowdate_omg[2*zdindx+1][i]+13000

    dayofsnowDiss = []
    for dosd in range (len(zerosnowdate_omg)/2):
        dayofsnowDiss.append([zerosnowdate_omg[2*dosd][0],zerosnowdate_omg[2*dosd+1][0]])
    
    dayofsnowDiss_df = pd.DataFrame(np.array(dayofsnowDiss))
    dayofsnowDiss_df.columns = ['2016','2017']
    
    return dayofsnowDiss_df

def dateTime(ncFlist): #2 years should be in 2 consecutive ncFile
    ncFdataset = []
    for ncfiles in ncFlist:
        ncFdataset.append(Dataset(ncfiles))
    
    timeFirstYear = ncFdataset[0].variables['time'][:] # get values
    timeSecondYear = ncFdataset[1].variables['time'][:] # get values
    time = np.concatenate((timeFirstYear,timeSecondYear), axis=0)

    t_unit = ncFdataset[0].variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

    try :

        t_cal = ncFdataset[0].variables['time'].calendar

    except AttributeError : # Attribute doesn't exist

        t_cal = u"gregorian" # or standard

    tvalue = num2date(time, units=t_unit, calendar=t_cal)
    date = [i.strftime("%Y-%m-%d %H:%M") for i in tvalue] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")  
    
    return date

def readNcfdatasetF0rEachVariableT0dataframe(ncFlist,variableName,hrunameDF,date):
    ncFdataset = []
    for ncfiles in ncFlist:
        ncFdataset.append(Dataset(ncfiles))
    
    variableList = []
    for ds in ncFdataset:
        variableList.append(pd.DataFrame(ds[variableName][:]))

    variable_2yearcons = []
    for dfs in range (len(variableList)/2):
        variable_2yearcons.append(pd.concat([variableList[2*dfs],variableList[2*dfs+1]], ignore_index=True))
    
    variable_df = pd.concat (variable_2yearcons, axis=1)
    variable_df.columns = hrunameDF[0]

    variable_df.set_index(pd.DatetimeIndex(date),inplace=True)
    counter = pd.DataFrame(np.arange(0,np.size(date)),columns=['counter'])
    counter.set_index(variable_df.index,inplace=True)
    variable_df2 = pd.concat([counter, variable_df], axis=1)
    
    return variable_df2
#%% SWE data
with open("input_SWE.csv") as scvd:
    reader = csv.reader(scvd)
    raw_swe = [r for r in reader]
sc_swe_column = []
for csv_counter1 in range (len (raw_swe)):
    for csv_counter2 in range (2):
        sc_swe_column.append(raw_swe[csv_counter1][csv_counter2])
sc_swe=np.reshape(sc_swe_column,(len (raw_swe),2))
sc_swe = sc_swe[1:]
sc_swe_obs_date = pd.DatetimeIndex(sc_swe[:,0])
sc_swe_obs = [float(value) for value in sc_swe[:,1]]
swe_obs_df = pd.DataFrame(sc_swe_obs, columns = ['observed swe']) 
swe_obs_df.set_index(sc_swe_obs_date,inplace=True)
#%%
p1 = [273.66,273.75] #273.66 tempCritRain	
p2 = [1.05,1.1] # 1.045 frozenPrecipMultip	
p3 = [2,3] #2, 3, 4] #mw_exp exponent for meltwater flow
p4 = [0.89,0.94] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500
p5 = [0.89,0.94] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p6 = [0.68,0.75] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500
p7 = [0.75,0.8] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p8 = [0.35,0.45] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p9 = [0.3,0.5] #0.5albedoSootLoad

def hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9):#, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21):
    ix1 = np.arange(1,len(p1)+1)
    ix2 = np.arange(1,len(p2)+1)
    ix3 = np.arange(1,len(p3)+1)
    ix4 = np.arange(1,len(p4)+1)
    ix5 = np.arange(1,len(p5)+1)
    ix6 = np.arange(1,len(p6)+1)
    ix7 = np.arange(1,len(p7)+1)
    ix8 = np.arange(1,len(p8)+1)
    ix9 = np.arange(1,len(p9)+1)

    c = list(itertools.product(ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9))#,ix10,ix11,ix12,ix13,ix14,ix15,ix16,ix17,ix18,ix19,ix20,ix21))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9)#, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21)
hru_num = np.size(hruidxID)

out_names = ['pt11', 'pt12', 'pt13', 'pt21', 'pt22', 'pt23']
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df0 = pd.DataFrame (hru_names1)

#%% reading output_swe files for open scenario
av_ncfiles = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp11_2015-2016_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp11_2016-2017_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp12_2015-2016_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp12_2016-2017_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp13_2015-2016_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp13_2016-2017_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp21_2015-2016_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp21_2016-2017_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp22_2015-2016_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp22_2016-2017_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp23_2015-2016_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_bp23_2016-2017_senatorVariableDecayRate_1.nc"]

av_all = []
for ncfiles in av_ncfiles:
    av_all.append(Dataset(ncfiles))

av_swe = []
for dfs in av_all:
    av_swe.append(pd.DataFrame(dfs['scalarSWE'][:]))

av_swe_2yearcons = []
for dfs2 in range (len(av_swe)/2):
    av_swe_2yearcons.append(pd.concat([av_swe[2*dfs2],av_swe[2*dfs2+1]], ignore_index=True))
av_swe_df = pd.concat (av_swe_2yearcons, axis=1)
av_swe_df.columns = hru_names_df0[0]

#%% output time step
TimeSc16 = av_all[0].variables['time'][:] # get values
TimeSc17 = av_all[1].variables['time'][:] # get values
TimeSc = np.concatenate((TimeSc16,TimeSc17), axis=0)

t_unitSc = av_all[0].variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = av_all[0].variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueSc = num2date(TimeSc, units=t_unitSc, calendar=t_cal)
DateSc = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueSc] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")    
#%%
av_swe_df.set_index(pd.DatetimeIndex(DateSc),inplace=True)
counter = pd.DataFrame(np.arange(0,np.size(DateSc)),columns=['counter'])
counter.set_index(av_swe_df.index,inplace=True)
av_swe_df2 = pd.concat([counter, av_swe_df], axis=1)
#%%   calculating day of snow disapperance open scenario
dayofsnowDiss_df=calculateDay0fSnowDissappearance(av_swe_df2,hru_names_df0)

dayofsnowDiss_obs = np.array([np.array([4686]),np.array([14430])]).T
dayofsnowDiss_obs_df = pd.DataFrame(dayofsnowDiss_obs,columns=['2016','2017'])

dosd_residual=[]
for years in dayofsnowDiss_df.columns:
    dosd_residual.append(abs(dayofsnowDiss_obs_df[years][0]-dayofsnowDiss_df[years])/24)

dosd_residual_df = pd.DataFrame(np.reshape(np.array(dosd_residual),(2,3072)).T, columns=['dosd2016','dosd2017'])
#%%# scenario 2 (veg1)
p8 = [4.5,6.6,8] #refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)
p9 = [0.3,0.45,0.6] #throughfallScaleSnow
p10 = [700,874,950] #specificHeatVeg   j/kg k         |     874.0000 |     500.0000 |    1500.0000
p11 = [0.02,0.04,0.06] #leafDimension             |       0.0400 |       0.0100 |       0.1000

def hru_ix_ID(p1, p2, p3, p4):#, p5, p6, p7, p8, p9, p10, p11):#, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21):
    ix1 = np.arange(1,len(p1)+1)
    ix2 = np.arange(1,len(p2)+1)
    ix3 = np.arange(1,len(p3)+1)
    ix4 = np.arange(1,len(p4)+1)

    c = list(itertools.product(ix1,ix2,ix3,ix4))#,,ix5,ix6,ix7,ix8,ix9,ix10,ix11ix12,ix13,ix14,ix15,ix16,ix17,ix18,ix19,ix20,ix21))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p8, p9, p10, p11)#,, p5, p6, p7, p8, p9, p10, p11 p12, p13, p14, p15, p16, p17, p18, p19, p20, p21)
#
hru_num = np.size(hruidxID)

out_names = ['ssc','dlb','dlc','dsb','dsc','slb','slc','ssb']
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df1 = pd.DataFrame (hru_names1)
#%% reading output_swe files for veg scenario1
av_ncfiles_vg1 = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1ssc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1ssc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dlb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dlb_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dlc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dlc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dsb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dsb_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dsc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1dsc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1slb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1slb_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1slc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1slc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1ssb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS1ssb_2016-2017_senatorVariableDecayRate_1.nc",
                  ]

av_all_vg1 = []
for ncfiles1 in av_ncfiles_vg1:
    av_all_vg1.append(Dataset(ncfiles1))

av_swe_vg1 = []
for dfs1 in av_all_vg1:
    av_swe_vg1.append(pd.DataFrame(dfs1['scalarSWE'][:]))

av_swe_2yearcons_vg1 = []
for dfs21 in range (len(av_swe_vg1)/2):
    av_swe_2yearcons_vg1.append(pd.concat([av_swe_vg1[2*dfs21],av_swe_vg1[2*dfs21+1]], ignore_index=True))
av_swe_df_vg1 = pd.concat (av_swe_2yearcons_vg1, axis=1)
av_swe_df_vg1.columns = hru_names_df1[0]

av_swe_df_vg1.set_index(pd.DatetimeIndex(DateSc),inplace=True)
av_swe_df_vg12 = pd.concat([counter, av_swe_df_vg1], axis=1)

#%%   calculating day of snow disapperance veg1 scenario
dayofsnowDiss_vg1_df=calculateDay0fSnowDissappearance(av_swe_df_vg12,hru_names_df1)

#%% scalarLWNetGround

av_nlwrG_vg1 = []
for dfs1 in av_all_vg1:
    av_nlwrG_vg1.append(pd.DataFrame(dfs1['scalarLWNetGround'][:]))

av_nlwrG_2yearcons_vg1 = []
for dfs21 in range (len(av_nlwrG_vg1)/2):
    av_nlwrG_2yearcons_vg1.append(pd.concat([av_nlwrG_vg1[2*dfs21],av_nlwrG_vg1[2*dfs21+1]], ignore_index=True))
av_nlwrG_df_vg1 = pd.concat (av_nlwrG_2yearcons_vg1, axis=1)
av_nlwrG_df_vg1.columns = hru_names_df1[0]

av_nlwrG_df_vg1.set_index(pd.DatetimeIndex(DateSc),inplace=True)
av_nlwrG_df_vg12 = pd.concat([counter, av_nlwrG_df_vg1], axis=1)

av_nlwrG_df_vg12016 = av_nlwrG_df_vg12[0:8785]
av_nlwrG_df_vg12017 = av_nlwrG_df_vg12[8784:]

av_nlwrG_df_vg12016SS = av_nlwrG_df_vg12016[0:6217]
av_nlwrG_df_vg12017SS = av_nlwrG_df_vg12017[0:6193]

sumNlwrG_df_vg12016SS = (av_nlwrG_df_vg12016SS.sum(axis=0)[1:])/1000
sumNlwrG_df_vg12017SS = (av_nlwrG_df_vg12017SS.sum(axis=0)[1:])/1000

#%% scalarGroundAbsorbedSolar
av_nswrG_vg1 = []
for dfs11 in av_all_vg1:
    av_nswrG_vg1.append(pd.DataFrame(dfs11['scalarGroundAbsorbedSolar'][:]))

av_nswrG_2yearcons_vg1 = []
for dfs211 in range (len(av_nswrG_vg1)/2):
    av_nswrG_2yearcons_vg1.append(pd.concat([av_nswrG_vg1[2*dfs211],av_nswrG_vg1[2*dfs211+1]], ignore_index=True))
av_nswrG_df_vg1 = pd.concat (av_nswrG_2yearcons_vg1, axis=1)
av_nswrG_df_vg1.columns = hru_names_df1[0]

av_nswrG_df_vg1.set_index(pd.DatetimeIndex(DateSc),inplace=True)
av_nswrG_df_vg12 = pd.concat([counter, av_nswrG_df_vg1], axis=1)

av_nswrG_df_vg12016 = av_nswrG_df_vg12[0:8785]
av_nswrG_df_vg12017 = av_nswrG_df_vg12[8784:]

av_nswrG_df_vg12016SS = av_nswrG_df_vg12016[0:6217]
av_nswrG_df_vg12017SS = av_nswrG_df_vg12017[0:6193]

sumNswrG_df_vg12016SS = (av_nswrG_df_vg12016SS.sum(axis=0)[1:])/1000
sumNswrG_df_vg12017SS = (av_nswrG_df_vg12017SS.sum(axis=0)[1:])/1000

#%% reading output_swe files for veg scenario2
av_ncfiles_vg2 = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2ssc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2ssc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dlb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dlb_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dlc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dlc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dsb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dsb_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dsc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2dsc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2slb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2slb_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2slc_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2slc_2016-2017_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2ssb_2015-2016_senatorVariableDecayRate_1.nc",
                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T1_vegS2ssb_2016-2017_senatorVariableDecayRate_1.nc"]
av_all_vg2 = []
for ncfiles2 in av_ncfiles_vg2:
    av_all_vg2.append(Dataset(ncfiles2))

av_swe_vg2 = []
for dfs2 in av_all_vg2:
    av_swe_vg2.append(pd.DataFrame(dfs2['scalarSWE'][:]))

av_swe_2yearcons_vg2 = []
for dfs22 in range (len(av_swe_vg2)/2):
    av_swe_2yearcons_vg2.append(pd.concat([av_swe_vg2[2*dfs22],av_swe_vg2[2*dfs22+1]], ignore_index=True))
av_swe_df_vg2 = pd.concat (av_swe_2yearcons_vg2, axis=1)
av_swe_df_vg2.columns = hru_names_df1[0]

av_swe_df_vg2.set_index(pd.DatetimeIndex(DateSc),inplace=True)
av_swe_df_vg22 = pd.concat([counter, av_swe_df_vg2], axis=1)

#%%   calculating day of snow disapperance veg2 scenario
dayofsnowDiss_vg2_df=calculateDay0fSnowDissappearance(av_swe_df_vg22,hru_names_df1)

#%% scalarGroundNetNrgFlux
av_nlwrG_vg2 = []
for dfs2 in av_all_vg2:
    av_nlwrG_vg2.append(pd.DataFrame(dfs2['scalarLWNetGround'][:]))

av_nlwrG_2yearcons_vg2 = []
for dfs22 in range (len(av_nlwrG_vg2)/2):
    av_nlwrG_2yearcons_vg2.append(pd.concat([av_nlwrG_vg2[2*dfs22],av_nlwrG_vg2[2*dfs22+1]], ignore_index=True))
av_nlwrG_df_vg2 = pd.concat (av_nlwrG_2yearcons_vg2, axis=1)
av_nlwrG_df_vg2.columns = hru_names_df1[0]

av_nlwrG_df_vg2.set_index(pd.DatetimeIndex(DateSc),inplace=True)
av_nlwrG_df_vg22 = pd.concat([counter, av_nlwrG_df_vg2], axis=1)

av_nlwrG_df_vg22016 = av_nlwrG_df_vg22[0:8785]
av_nlwrG_df_vg22017 = av_nlwrG_df_vg22[8784:]

av_nlwrG_df_vg22016SS = av_nlwrG_df_vg22016[0:6217]
av_nlwrG_df_vg22017SS = av_nlwrG_df_vg22017[0:6193]

sumNlwrG_df_vg22016SS = (av_nlwrG_df_vg22016SS.sum(axis=0)[1:])/1000 #kw/m2
sumNlwrG_df_vg22017SS = (av_nlwrG_df_vg22017SS.sum(axis=0)[1:])/1000

#%% scalarGroundAbsorbedSolar
av_nswrG_vg2 = []
for dfs22 in av_all_vg2:
    av_nswrG_vg2.append(pd.DataFrame(dfs22['scalarGroundAbsorbedSolar'][:]))

av_nswrG_2yearcons_vg2 = []
for dfs222 in range (len(av_nswrG_vg2)/2):
    av_nswrG_2yearcons_vg2.append(pd.concat([av_nswrG_vg2[2*dfs222],av_nswrG_vg2[2*dfs222+1]], ignore_index=True))
av_nswrG_df_vg2 = pd.concat (av_nswrG_2yearcons_vg2, axis=1)
av_nswrG_df_vg2.columns = hru_names_df1[0]

av_nswrG_df_vg2.set_index(pd.DatetimeIndex(DateSc),inplace=True)
av_nswrG_df_vg22 = pd.concat([counter, av_nswrG_df_vg2], axis=1)

av_nswrG_df_vg22016 = av_nswrG_df_vg22[0:8785]
av_nswrG_df_vg22017 = av_nswrG_df_vg22[8784:]

av_nswrG_df_vg22016SS = av_nswrG_df_vg22016[0:6217]
av_nswrG_df_vg22017SS = av_nswrG_df_vg22017[0:6193]

sumNswrG_df_vg22016SS = (av_nswrG_df_vg22016SS.sum(axis=0)[1:])/1000 #kw/m2
sumNswrG_df_vg22017SS = (av_nswrG_df_vg22017SS.sum(axis=0)[1:])/1000
#%%  boxplot  dayofsnowDiss_obs_df
d1 = [dayofsnowDiss_vg1_df['2016']/24,(dayofsnowDiss_vg1_df['2017']-8784)/24,dayofsnowDiss_vg2_df['2016']/24,(dayofsnowDiss_vg2_df['2017']-8784)/24]#,]
fig = plt.subplots(1,1, figsize=(20,15))
bp1 = plt.boxplot(d1, patch_artist=True)
bp1['boxes'][0].set(color='navy', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][1].set(color='blue', linewidth=2, facecolor = 'olive', hatch = '/')
bp1['boxes'][2].set(color='green', linewidth=2, facecolor = 'pink', hatch = '/')
bp1['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][4].set(color='turquoise', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][5].set(color='darkturquoise', linewidth=2, facecolor = 'pink', hatch = '/')


plt.xticks([1,2,3,4], ['2016veg1','2017veg1','2016veg2','2017veg2'], fontsize=20)#'open2016','open2017',
plt.yticks(fontsize=20)
plt.xlabel('Scenarios', fontsize=30)
plt.ylabel('day of snow dissappearance', fontsize=30)
plt.savefig('dosdVeg12.png')


#%%  boxplot  scalarGroundNetNrgFlux
d2 = [sumNlwrG_df_vg12016SS,sumNlwrG_df_vg12017SS,sumNswrG_df_vg12016SS,sumNswrG_df_vg12017SS,sumNlwrG_df_vg22016SS,sumNlwrG_df_vg22017SS,sumNswrG_df_vg22016SS,sumNswrG_df_vg22017SS]#sumNefG_df_vg12016SS,sumNefG_df_vg12017SS,
fig2 = plt.subplots(1,1, figsize=(20,15))
bp2 = plt.boxplot(d2, patch_artist=True)
bp2['boxes'][0].set(color='blue', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][1].set(color='navy', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][2].set(color='red', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][3].set(color='orange', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][4].set(color='blue', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][5].set(color='navy', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][6].set(color='red', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][7].set(color='orange', linewidth=2, facecolor = 'pink', hatch = '/')

plt.xticks([1,2,3,4,5,6,7,8], ['lwr16veg1','lwr17veg1','swr16veg1','swr17veg1','lwr16veg2','lwr17veg2','swr16veg2','swr17veg2'], fontsize=20)#'2016veg1','2017veg1',
plt.yticks(fontsize=25)
plt.xlabel('Scenarios', fontsize=35)
plt.ylabel('Net radiation on the ground (kw/m2)', fontsize=30)
plt.savefig('netRVeg12.png')
#%% plotting
#DateSc_2 = [i.strftime("%Y-%m") for i in tvalueSc]
#sax = np.arange(0,np.size(DateSc))
#sa_xticks = DateSc_2
#safig, saax = plt.subplots(1,1, figsize=(20,15))
#plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
#saax.xaxis.set_major_locator(ticker.AutoLocator())
#plt.yticks(fontsize=20)
#for hru in av_swe_df_vg1.columns:
#    plt.plot(av_swe_df_vg1[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#
#plt.plot(swe_obs_df, 'k', markersize=10)
##
#plt.title('SWE (mm)', position=(0.04, 0.88), ha='left', fontsize=40)
#plt.xlabel('Time 2015-2017', fontsize=30)
#plt.ylabel('SWE(mm)', fontsize=30)
##plt.legend()
##plt.show()
#plt.savefig('initialResults_sc\swe_ScT1_scveg1.png')





















