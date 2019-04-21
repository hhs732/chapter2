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
#%% SWE observed data T4
with open("hhs_jemezVcm_Swe2.csv") as scvd:
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

swe_obs_hr = swe_obs_df.resample('H').mean()

maxSwe1 = swe_obs_df['observed swe'][0:7250].max()
maxSwe2 = swe_obs_df['observed swe'][7250:].max()
maxSwe_date1 = swe_obs_df['observed swe'][0:7250].idxmax()
maxSwe_date2 = swe_obs_df['observed swe'][7250:].idxmax()

#%% calib1
#p1 = [3,4] #mw_exp exponent for meltwater flow
#p2 = [0.001,0.002] #zsnow
#
#p3 = [0.84,0.86] #albedoMax 
#p4 = [0.87,0.95] #0.80,albedoMaxVisible 
#p5 = [0.7,0.75] #albedoMinVisible 
#p6 = [0.6,0.65] #albedoMaxNearIR 
#p7 = [0.3,0.4] #albedoMinNearIR 
#p8 = [500000,700000,1000000]#albedoDecayRate
#p9 = [1,3] #albedoRefresh 
#p10 = [0.3] #,0.5albedoSootLoad
#     
#p11 = [1] #rootingDepth
#p12 = [0.02] #winterSAI
#p13 = [0.1] #summerLAI
#p14 = [0.01] #LAIMIN
#p15 = [1] #LAIMAX
#p16 = [0.3] #heightCanopyTop
#p17 = [0.03] #heightCanopyBottom
#p18 = [0.02] #0.1#z0Canopy  
#
##p4 = [0.28] #,0.4windReductionParam
##p1 = [273.16] #tempCritRain
##p13 = [0.45] #albedoMinSpring           |       0.5000 |       0.3000 |       1.0000
##p6 = [100] #newSnowDenMin 
##p2 = [1,1.05] # frozenPrecipMultip

#def hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9):
#    ix1 = np.arange(1,len(p1)+1)
#    ix2 = np.arange(1,len(p2)+1)
#    ix3 = np.arange(1,len(p3)+1)
#    ix4 = np.arange(1,len(p4)+1)
#    ix5 = np.arange(1,len(p5)+1)
#    ix6 = np.arange(1,len(p6)+1)
#    ix7 = np.arange(1,len(p7)+1)
#    ix8 = np.arange(1,len(p8)+1)
#    ix9 = np.arange(1,len(p9)+1)
#   
#    c = list(itertools.product(ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9))
#    ix_numlist=[]
#    for tup in c:
#        ix_numlist.append(''.join(map(str, tup)))
#    new_list = [float(i) for i in ix_numlist]
#
#    return(new_list)  
#
#hruidxID = hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9)#
#
hruidxID =[1001]
hru_num = np.size(hruidxID)

out_names = ['senator','senator2']#'llj','llt','lsj','lst','slj','lmj','lmt','slt','smj','smt','ssj','sst',
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df0 = pd.DataFrame (hru_names1)
#%% reading output_swe files for open scenario
av_ncfiles = ["C:\Users\HHS\summaTestCases_2.x\output\Jemez\Jemez_senator_2009-2010_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\Jemez\Jemez_senator_2010-2011_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\Jemez\Jemez_senator2_2009-2010_senatorVariableDecayRate_1.nc",
              "C:\Users\HHS\summaTestCases_2.x\output\Jemez\Jemez_senator2_2010-2011_senatorVariableDecayRate_1.nc",
               ]

Date, tvalue = dateTime(av_ncfiles)

av_swe_df = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles,'scalarSWE',hru_names_df0,Date)
#%% plotting
#sax = np.arange(0,np.size(Date))
#date = [i.strftime("%Y-%m") for i in tvalue]
#sa_xticks = date
#safig, saax = plt.subplots(1,1, figsize=(20,15))
#plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
#saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.subplots(1,1, figsize=(20,15))
plt.xticks(rotation=25, fontsize=20)
plt.yticks(fontsize=20)
for hru in hru_names_df0[0]:
    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#
plt.plot(swe_obs_hr, 'k', markersize=10)
#
plt.title('SWE (mm); Jemez_VCM', position=(0.04, 0.88), ha='left', fontsize=40)
plt.xlabel('Time 2009-2011', fontsize=30)
plt.ylabel('SWE(mm)', fontsize=30)
#plt.legend()
#plt.show()
plt.savefig('swe_vcm_calib.png')

#%%   calculating day of snow disapperance open scenario
dosd0pen_df=calculateDay0fSnowDissappearance(av_swe_df,hru_names_df0)

#dayofsnowDiss_obs = np.array([np.array([4656]),np.array([14424])]).T
#dayofsnowDiss_obs_df = pd.DataFrame(dayofsnowDiss_obs,columns=['2016','2017'])
#
#dofs_residual=[]
#for year in dayofsnowDiss_obs_df.columns:
#    dofs_residual.append((dayofsnowDiss_obs_df[year][0]-dosd0pen_df[year])/24)
#
#dosd_residual_df = pd.DataFrame(np.reshape(np.array(dofs_residual),(2,len(dosd0pen_df))).T, columns=['dosd2016','dosd2017'])

#%% find 2 max swe in each year
#maxSweFY_1obs = maxSwe1.copy() #[294.64] # 15/3/2016 
#maxSweFY_1 = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df0[0],3966)
##maxSweFY_2obs = 1143.2 #31/3/2016 8am
##maxSweFY_2 = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df0[0],4376)
#
#maxSweSY_1obs = 619.76 #28/3/2017 
#maxSweSY_1 = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df0[0],13042)
#maxSweSY_2obs = maxSwe2.copy() # [670.56] # 14/4/2017 6am
#maxSweSY_2 = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df0[0],13450)
#
#maxSwe_obs_df = pd.DataFrame(np.array([maxSweFY_1obs,maxSweSY_1obs,maxSweSY_2obs])).T; maxSwe_obs_df.columns=['maxFY_1','maxSY_1','maxSY_2']
#maxSwe_df = pd.DataFrame(np.array([maxSweFY_1,maxSweSY_1,maxSweSY_2])).T; maxSwe_df.columns=['maxFY_1','maxSY_1','maxSY_2']
#
#maxSwe_residual=[]
#for crit in maxSwe_df.columns:
#    maxSwe_residual.append(abs(maxSwe_obs_df[crit][0]-maxSwe_df[crit]))
#
#maxSwe_residual_df = pd.DataFrame(np.reshape(np.array(maxSwe_residual),(3,len(maxSwe_df))).T, columns=['maxSweFY_1','maxSweSY_1','maxSweSY_2'])

#%% find best params
#criteria_df = pd.concat([maxSwe_residual_df, dosd_residual_df], axis=1) 
#criteria_df.set_index(hru_names_df0[0],inplace=True)
#pareto_model_param = pd.DataFrame(criteria_df.index[((criteria_df['dosd2016'])>=1) & ((criteria_df['dosd2017'])>=-2) &
#                                                    ((criteria_df['dosd2016'])<=2) & ((criteria_df['dosd2017'])<=0) &
#                                                    ((criteria_df['maxSweFY_1']) <= 89) & 
#                                                    ((criteria_df['maxSweSY_1'])<=75) & ((criteria_df['maxSweSY_2'])<=48)
#                                                    ].tolist()) 
## 3  2  48  40  34  49                                                   
##%%
#DateSc = [i.strftime("%Y-%m") for i in tvalue]
#sax = np.arange(0,np.size(DateSc))
#sa_xticks = DateSc
#safig, saax = plt.subplots(1,1, figsize=(20,15))
#plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
#saax.xaxis.set_major_locator(ticker.AutoLocator())
#plt.yticks(fontsize=20)
#for hru2 in pareto_model_param[0]:
#    plt.plot(av_swe_df[hru2])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#
#plt.plot(swe_obs_df, 'k', markersize=10)
#
#plt.title('SWE (mm) SCT1', position=(0.04, 0.88), ha='left', fontsize=40)
#plt.xlabel('Time 2015-2017', fontsize=30)
#plt.ylabel('SWE(mm)', fontsize=30)
#plt.legend()
#plt.savefig('sweScT1_best1.png')


















