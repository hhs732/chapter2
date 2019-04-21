###       /bin/bash runTestCases_dockerScT1_calib.sh
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
    
    return date, tvalue

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
#%% SWE observed data T1 open
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
#counter = pd.DataFrame(np.arange(0,len(swe_obs_df)), columns = ['counter']); counter.set_index(sc_swe_obs_date,inplace=True)
#swe_obs_df2 = pd.concat([counter, swe_obs_df], axis=1)

maxSwe1 = swe_obs_df['observed swe'][0:365].max()
maxSwe2 = swe_obs_df['observed swe'][365:].max()
maxSwe_date1 = swe_obs_df['observed swe'][0:365].idxmax()
maxSwe_date2 = swe_obs_df['observed swe'][365:].idxmax()

#%% SWE observed data T1 under canopy
with open("T1hill_snw_under.csv") as scvd2:
    reader2 = csv.reader(scvd2)
    raw_swe2 = [r2 for r2 in reader2]
sc_swe_column2 = []
for csv_counter12 in range (len (raw_swe2)):
    for csv_counter22 in range (2):
        sc_swe_column2.append(raw_swe2[csv_counter12][csv_counter22])
sc_swe_uc=np.reshape(sc_swe_column2,(len (raw_swe2),2))
sc_swe_uc = sc_swe_uc[1:]
sc_swe_obs_date_uc = pd.DatetimeIndex(sc_swe_uc[:,0])
sc_swe_obs_uc = [float(value) for value in sc_swe_uc[:,1]]
swe_obs_uc_df = pd.DataFrame(sc_swe_obs_uc, columns = ['observed swe']) 
swe_obs_uc_df.set_index(sc_swe_obs_date_uc,inplace=True)
swe_obs_uc_day = swe_obs_uc_df.resample('D').mean()


















