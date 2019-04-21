###       /bin/bash runTestCases_dockerScT4_calib.sh
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
with open("input_SWE_T4.csv") as scvd:
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

#snowPpt = []
#for indxswe in range (len(swe_obs_df)-1):
#    if (swe_obs_df['observed swe'][indxswe+1]-swe_obs_df['observed swe'][indxswe])<0:
#        snowPpt.append(0)
#    else: snowPpt.append(swe_obs_df['observed swe'][indxswe+1]-swe_obs_df['observed swe'][indxswe])
#snowPpt_df = pd.DataFrame(snowPpt, columns = ['snow ppt']) 
#snowPpt_df.set_index(sc_swe_obs_date[1:],inplace=True)  
#
#snowPpt_hr = snowPpt_df.resample('h').sum() 
#%% open scenario-Params test
p1 = [273.16] #,273.5,273.66 tempCritRain	
p2 = [2,3,4] #2, 3, 4] #mw_exp exponent for meltwater flow
p3 = [0.2,0.28,0.4] #windReductionParam        |       0.2800 |       0.0000 |       1.0000
p4 = [0.001,0.002] #zsnow

p5 = [0.85,0.9,0.94] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p6 = [0.6,0.7,0.75] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500
p7 = [0.6,0.7] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p8 = [0.3,0.4] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p9 = [1,3,5]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p10 = [350000,420000,500000]#,350000,200000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 

p11 = [0.3] #0.5albedoSootLoad
p12 = [0.86] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500
p13 = [0.01] #winterSAI
p14 = [0.1] #summerLAI
p15 = [0.01] #LAIMIN
p16 = [1] #LAIMAX
p17 = [0.3] #heightCanopyTop
p18 = [0.03] #heightCanopyBottom to calculate wind speed, wind speed reduction; coupute roughness length of the veg canopy; neutral ground resistance; canopy air temp;
p19 = [0.02] #0.1#z0Canopy  

def hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9):#, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21
    ix1 = np.arange(1,len(p1)+1)
    ix2 = np.arange(1,len(p2)+1)
    ix3 = np.arange(1,len(p3)+1)
    ix4 = np.arange(1,len(p4)+1)
    ix5 = np.arange(1,len(p5)+1)
    ix6 = np.arange(1,len(p6)+1)
    ix7 = np.arange(1,len(p7)+1)
    ix8 = np.arange(1,len(p8)+1)
    ix9 = np.arange(1,len(p9)+1)
    
    c = list(itertools.product(ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9))#,ix10,ix11,ix12,ix13,ix14,ix15,ix16,ix17,ix18,ix19,ix20,ix21
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p2, p3, p4, p5, p6, p7, p8, p9, p10)#, p11, p12, p13, p14, p15, p16, p17 p18, p19, p20, p21
#
hru_num = np.size(hruidxID)

out_names = ['llj']
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df0 = pd.DataFrame (hru_names1)

#%% reading output_swe files for open scenario
av_ncfiles = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\scT4_openCal_llj_2015-2016_senatorVariableDecayRate_1.nc",
                "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\scT4_openCal_llj_2016-2017_senatorVariableDecayRate_1.nc"
                ]

Date, tvalue = dateTime(av_ncfiles)

av_swe_df = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles,'scalarSWE',hru_names_df0,Date)
#%% plotting
sax = np.arange(0,np.size(Date))
sa_xticks = Date
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=20)
for hru in hru_names_df0[0][3000:5831]:
    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]

plt.plot(swe_obs_df, 'k', markersize=10)
#
plt.title('SWE (mm); SagehenT4', position=(0.04, 0.88), ha='left', fontsize=40)
plt.xlabel('Time 2015-2017', fontsize=30)
plt.ylabel('SWE(mm)', fontsize=30)
#plt.legend()
#plt.show()
plt.savefig('swe_ScT4_open2.png')

#%%   calculating day of snow disapperance open scenario
dayofsnowDiss0penT4_df=calculateDay0fSnowDissappearance(av_sweT4_df,hru_names_df0)

dayofsnowDiss_obsT4 = np.array([np.array([4656]),np.array([14424])]).T
dayofsnowDiss_obsT4_df = pd.DataFrame(dayofsnowDiss_obsT4,columns=['2016','2017'])

#%% find 2 max swe in each year
maxSwe2016_1obs = [276.86]
maxSwe2016_1 = readSpecificDatafromAllHRUs(av_swe_df2,hru_names_df[0],3417)
maxSwe2016_2obs = [294.64]
maxSwe2016_2 = readSpecificDatafromAllHRUs(av_swe_df2,hru_names_df[0],3990)

maxSwe2017_1obs = [637.54]
maxSwe2017_1 = readSpecificDatafromAllHRUs(av_swe_df2,hru_names_df[0],12558)
maxSwe2017_2obs = [670.56]
maxSwe2017_2 = readSpecificDatafromAllHRUs(av_swe_df2,hru_names_df[0],13488)

maxSwe_obs_df = pd.DataFrame(np.array([maxSwe2016_1obs,maxSwe2016_2obs,maxSwe2017_1obs,maxSwe2017_2obs]).T,columns=['2016_1','2016_2','2017_1','2017_2'])
maxSwe_df = pd.DataFrame(np.array([maxSwe2016_1,maxSwe2016_2,maxSwe2017_1,maxSwe2017_2]).T,columns=['2016_1','2016_2','2017_1','2017_2'])

maxSwe_residual=[]
for years2 in maxSwe_df.columns:
    maxSwe_residual.append(abs(maxSwe_obs_df[years2][0]-maxSwe_df[years2]))

maxSwe_residual_df = pd.DataFrame(np.reshape(np.array(maxSwe_residual),(4,3072)).T, columns=['maxSwe2016_1','maxSwe2016_2','maxSwe2017_1','maxSwe2017_2'])

plt.subplots(1,1, figsize=(20,15))
for hru in av_swe_df.columns:
    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
plt.plot(swe_obs_df, 'k', markersize=10)

plt.plot(maxSwe_df['2016_1'])
plt.plot(maxSwe_obs_df['2016_1'], 'ko')
#%% find best params
criteria_df = pd.concat([maxSwe_residual_df, dosd_residual_df], axis=1) 
criteria_df.set_index(hru_names_df[0],inplace=True)
pareto_model_param = pd.DataFrame(criteria_df.index[((criteria_df['maxSwe2016_1']) <= 50) & ((criteria_df['maxSwe2016_2'])<=75) & 
                                                    ((criteria_df['maxSwe2017_1'])<=150) & ((criteria_df['maxSwe2017_2'])<=150) &
                                                    ((criteria_df['dosd2016'])<=10) & ((criteria_df['dosd2017'])<=9)  ].tolist())#

DateSc_2 = [i.strftime("%Y-%m") for i in tvalueSc]
sax = np.arange(0,np.size(DateSc))
sa_xticks = DateSc_2
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=20)
for hru2 in pareto_model_param[0]:
    plt.plot(av_swe_df[hru2])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]

plt.plot(swe_obs_df, 'k', markersize=10)

plt.title('SWE (mm) SCT1', position=(0.04, 0.88), ha='left', fontsize=40)
plt.xlabel('Time 2015-2017', fontsize=30)
plt.ylabel('SWE(mm)', fontsize=30)
plt.legend()
plt.savefig('initialResults_sc\sweScT1_TPBP.png')


















#%% scalarLWNetGround
#
#av_nlwrG_df_open = readNcfdatasetF0rEachVariableT0dataframe(av_ncfilesT4,'scalarLWNetGround',hru_names_df0,DateSc)
#
#av_nlwrG_df_open2016 = av_nlwrG_df_open[0:8785]
#av_nlwrG_df_open2017 = av_nlwrG_df_open[8784:]
#
#av_nlwrG_df_open2016SS = av_nlwrG_df_open2016[0:6217]
#av_nlwrG_df_open2017SS = av_nlwrG_df_open2017[0:6193]
#
#sumNlwrG_df_open2016SS = (av_nlwrG_df_open2016SS.sum(axis=0)[1:])/1000
#sumNlwrG_df_open2017SS = (av_nlwrG_df_open2017SS.sum(axis=0)[1:])/1000
#
##%% scalarGroundAbsorbedSolar
#
#av_nswrG_df_open = readNcfdatasetF0rEachVariableT0dataframe(av_ncfilesT4,'scalarGroundAbsorbedSolar',hru_names_df0,DateSc)
#
#av_nswrG_df_open2016 = av_nswrG_df_open[0:8785]
#av_nswrG_df_open2017 = av_nswrG_df_open[8784:]
#
#av_nswrG_df_open2016SS = av_nswrG_df_open2016[0:6217]
#av_nswrG_df_open2017SS = av_nswrG_df_open2017[0:6193]
#
#sumNswrG_df_open2016SS = (av_nswrG_df_open2016SS.sum(axis=0)[1:])/1000
#sumNswrG_df_open2017SS = (av_nswrG_df_open2017SS.sum(axis=0)[1:])/1000
#
##%% reading output_swe files for veg scenario1
#av_ncfiles_vg1 = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1ssc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1ssc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dlb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dlb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dlc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dlc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dsb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dsb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dsc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1dsc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1slb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1slb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1slc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1slc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1ssb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg1ssb_2016-2017_senatorVariableDecayRate_1.nc",
#                  ]
#
#av_swe_df_vg1 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg1,'scalarSWE',hru_names_df1,DateSc)
#
##%%   calculating day of snow disapperance veg1 scenario
#av_swe_df4000 = av_swe_df_vg1[:][4000:8784]
#av_swe_df13000 = av_swe_df_vg1[:][13000:17137]
#    
#zerosnowdate = []
#for val in hru_names_df1[0]:
#    zerosnowdate.append(np.where(av_swe_df4000[val]==0))
#    zerosnowdate.append(np.where(av_swe_df13000[val]==0))
#
#zerosnowdate_omg = [item[0] for item in zerosnowdate] #change tuple to array
#
#dayofsnowDiss = []
#for wow in range (len(zerosnowdate_omg)):
#    if len(zerosnowdate_omg[wow])>1:
#        dayofsnowDiss.append(zerosnowdate_omg[wow][0])
#    else: dayofsnowDiss.append(4700)
#    
#dayofsnowDiss_resh = np.reshape(dayofsnowDiss,(len(dayofsnowDiss)/2,2))
#    
#dayofsnowDiss_vg1_df = pd.DataFrame(np.array(dayofsnowDiss_resh))
#dayofsnowDiss_vg1_df.columns = ['2016','2017']
#dayofsnowDiss_vg1_df.index = hru_names_df1[0]
##%% scalarLWNetGround
#
#av_nlwrG_df_vg1 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg1,'scalarLWNetGround',hru_names_df1,DateSc)
#
#av_nlwrG_df_vg12016 = av_nlwrG_df_vg1[0:8785]
#av_nlwrG_df_vg12017 = av_nlwrG_df_vg1[8784:]
#
#av_nlwrG_df_vg12016SS = av_nlwrG_df_vg12016[0:6217]
#av_nlwrG_df_vg12017SS = av_nlwrG_df_vg12017[0:6193]
#
#sumNlwrG_df_vg12016SS = (av_nlwrG_df_vg12016SS.sum(axis=0)[1:])/1000
#sumNlwrG_df_vg12017SS = (av_nlwrG_df_vg12017SS.sum(axis=0)[1:])/1000
#
##%% scalarGroundAbsorbedSolar
#
#av_nswrG_df_vg1 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg1,'scalarGroundAbsorbedSolar',hru_names_df1,DateSc)
#
#av_nswrG_df_vg12016 = av_nswrG_df_vg1[0:8785]
#av_nswrG_df_vg12017 = av_nswrG_df_vg1[8784:]
#
#av_nswrG_df_vg12016SS = av_nswrG_df_vg12016[0:6217]
#av_nswrG_df_vg12017SS = av_nswrG_df_vg12017[0:6193]
#
#sumNswrG_df_vg12016SS = (av_nswrG_df_vg12016SS.sum(axis=0)[1:])/1000
#sumNswrG_df_vg12017SS = (av_nswrG_df_vg12017SS.sum(axis=0)[1:])/1000
#
##%% reading output_swe files for veg scenario2
#av_ncfiles_vg2 = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2ssc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2ssc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dlb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dlb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dlc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dlc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dsb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dsb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dsc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2dsc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2slb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2slb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2slc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2slc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2ssb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg2ssb_2016-2017_senatorVariableDecayRate_1.nc"]
#
#av_swe_df_vg2 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg2,'scalarSWE',hru_names_df1,DateSc)
#
##%%   calculating day of snow disapperance veg2 scenario
#av_swe_df40002 = av_swe_df_vg2[:][4000:8784]
#av_swe_df130002 = av_swe_df_vg2[:][13000:17137]
#    
#zerosnowdate2 = []
#for val2 in hru_names_df1[0]:
#    zerosnowdate2.append(np.where(av_swe_df40002[val]==0))
#    zerosnowdate2.append(np.where(av_swe_df130002[val]==0))
#
#zerosnowdate_omg2 = [item[0] for item in zerosnowdate2] #change tuple to array
#
#dayofsnowDiss2 = []
#for wow2 in range (len(zerosnowdate_omg2)):
#    if len(zerosnowdate_omg2[wow2])>1:
#        dayofsnowDiss2.append(zerosnowdate_omg2[wow2][0])
#    else: dayofsnowDiss2.append(4700)
#    
#dayofsnowDiss_resh2 = np.reshape(dayofsnowDiss2,(len(dayofsnowDiss2)/2,2))
#    
#dayofsnowDiss_vg2_df = pd.DataFrame(np.array(dayofsnowDiss_resh2))
#dayofsnowDiss_vg2_df.columns = ['2016','2017']
#dayofsnowDiss_vg2_df.index = hru_names_df1[0]
#
##%% scalarGroundNetNrgFlux
#
#av_nlwrG_df_vg2 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg2,'scalarLWNetGround',hru_names_df1,DateSc)
#
#av_nlwrG_df_vg22016 = av_nlwrG_df_vg2[0:8785]
#av_nlwrG_df_vg22017 = av_nlwrG_df_vg2[8784:]
#
#av_nlwrG_df_vg22016SS = av_nlwrG_df_vg22016[0:6217]
#av_nlwrG_df_vg22017SS = av_nlwrG_df_vg22017[0:6193]
#
#sumNlwrG_df_vg22016SS = (av_nlwrG_df_vg22016SS.sum(axis=0)[1:])/1000 #kw/m2
#sumNlwrG_df_vg22017SS = (av_nlwrG_df_vg22017SS.sum(axis=0)[1:])/1000
#
##%% scalarGroundAbsorbedSolar
#av_nswrG_df_vg2 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg2,'scalarGroundAbsorbedSolar',hru_names_df1,DateSc)
#
#av_nswrG_df_vg22016 = av_nswrG_df_vg2[0:8785]
#av_nswrG_df_vg22017 = av_nswrG_df_vg2[8784:]
#
#av_nswrG_df_vg22016SS = av_nswrG_df_vg22016[0:6217]
#av_nswrG_df_vg22017SS = av_nswrG_df_vg22017[0:6193]
#
#sumNswrG_df_vg22016SS = (av_nswrG_df_vg22016SS.sum(axis=0)[1:])/1000 #kw/m2
#sumNswrG_df_vg22017SS = (av_nswrG_df_vg22017SS.sum(axis=0)[1:])/1000
#
##%% reading output_swe files for veg scenario3
#av_ncfiles_vg3 = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3ssc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3ssc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dlb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dlb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dlc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dlc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dsb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dsb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dsc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3dsc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3slb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3slb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3slc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3slc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3ssb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg3ssb_2016-2017_senatorVariableDecayRate_1.nc"]
#
#av_swe_df_vg3 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg3,'scalarSWE',hru_names_df1,DateSc)
#
##%%   calculating day of snow disapperance veg3 scenario
#av_swe_df40003 = av_swe_df_vg3[:][4000:8784]
#av_swe_df130003 = av_swe_df_vg3[:][13000:17137]
#    
#zerosnowdate3 = []
#for val3 in hru_names_df1[0]:
#    zerosnowdate3.append(np.where(av_swe_df40003[val3]==0))
#    zerosnowdate3.append(np.where(av_swe_df130003[val3]==0))
#
#zerosnowdate_omg3 = [item[0] for item in zerosnowdate3] #change tuple to array
#
#dayofsnowDiss3 = []
#for wow3 in range (len(zerosnowdate_omg3)):
#    if len(zerosnowdate_omg3[wow3])>1:
#        dayofsnowDiss3.append(zerosnowdate_omg3[wow3][0])
#    else: dayofsnowDiss3.append(4700)
#    
#dayofsnowDiss_resh3 = np.reshape(dayofsnowDiss3,(len(dayofsnowDiss3)/2,2))
#    
#dayofsnowDiss_vg3_df = pd.DataFrame(np.array(dayofsnowDiss_resh3))
#dayofsnowDiss_vg3_df.columns = ['2016','2017']
#dayofsnowDiss_vg3_df.index = hru_names_df1[0]
#
##%% scalarGroundNetNrgFlux
#
#av_nlwrG_df_vg3 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg3,'scalarLWNetGround',hru_names_df1,DateSc)
#
#av_nlwrG_df_vg32016 = av_nlwrG_df_vg3[0:8785]
#av_nlwrG_df_vg32017 = av_nlwrG_df_vg3[8784:]
#
#av_nlwrG_df_vg32016SS = av_nlwrG_df_vg32016[0:6217]
#av_nlwrG_df_vg32017SS = av_nlwrG_df_vg32017[0:6193]
#
#sumNlwrG_df_vg32016SS = (av_nlwrG_df_vg32016SS.sum(axis=0)[1:])/1000 #kw/m2
#sumNlwrG_df_vg32017SS = (av_nlwrG_df_vg32017SS.sum(axis=0)[1:])/1000
#
##%% scalarGroundAbsorbedSolar
#av_nswrG_df_vg3 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg3,'scalarGroundAbsorbedSolar',hru_names_df1,DateSc)
#
#av_nswrG_df_vg32016 = av_nswrG_df_vg3[0:8785]
#av_nswrG_df_vg32017 = av_nswrG_df_vg3[8784:]
#
#av_nswrG_df_vg32016SS = av_nswrG_df_vg32016[0:6217]
#av_nswrG_df_vg32017SS = av_nswrG_df_vg32017[0:6193]
#
#sumNswrG_df_vg32016SS = (av_nswrG_df_vg32016SS.sum(axis=0)[1:])/1000 #kw/m2
#sumNswrG_df_vg32017SS = (av_nswrG_df_vg32017SS.sum(axis=0)[1:])/1000
#
##%% reading output_swe files for veg scenario4
#av_ncfiles_vg4 = ["C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4ssc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4ssc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dlb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dlb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dlc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dlc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dsb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dsb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dsc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4dsc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4slb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4slb_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4slc_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4slc_2016-2017_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4ssb_2015-2016_senatorVariableDecayRate_1.nc",
#                  "C:\Users\HHS\summaTestCases_2.x\output\sagehencreek\sagehen_T4_veg4ssb_2016-2017_senatorVariableDecayRate_1.nc"]
#
#av_swe_df_vg4 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg4,'scalarSWE',hru_names_df1,DateSc)
#
##%%   calculating day of snow disapperance veg4 scenario
#av_swe_df40004 = av_swe_df_vg4[:][4000:8784]
#av_swe_df130004 = av_swe_df_vg4[:][13000:17137]
#    
#zerosnowdate4 = []
#for val4 in hru_names_df1[0]:
#    zerosnowdate4.append(np.where(av_swe_df40004[val4]==0))
#    zerosnowdate4.append(np.where(av_swe_df130004[val4]==0))
#
#zerosnowdate_omg4 = [item[0] for item in zerosnowdate4] #change tuple to array
#
#dayofsnowDiss4 = []
#for wow4 in range (len(zerosnowdate_omg4)):
#    if len(zerosnowdate_omg4[wow4])>1:
#        dayofsnowDiss4.append(zerosnowdate_omg4[wow4][0])
#    else: dayofsnowDiss4.append(4700)
#    
#dayofsnowDiss_resh4 = np.reshape(dayofsnowDiss4,(len(dayofsnowDiss4)/2,2))
#    
#dayofsnowDiss_vg4_df = pd.DataFrame(np.array(dayofsnowDiss_resh4))
#dayofsnowDiss_vg4_df.columns = ['2016','2017']
#dayofsnowDiss_vg4_df.index = hru_names_df1[0]
#
##%% scalarGroundNetNrgFlux
#
#av_nlwrG_df_vg4 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg4,'scalarLWNetGround',hru_names_df1,DateSc)
#
#av_nlwrG_df_vg42016 = av_nlwrG_df_vg4[0:8785]
#av_nlwrG_df_vg42017 = av_nlwrG_df_vg4[8784:]
#
#av_nlwrG_df_vg42016SS = av_nlwrG_df_vg42016[0:6217]
#av_nlwrG_df_vg42017SS = av_nlwrG_df_vg42017[0:6193]
#
#sumNlwrG_df_vg42016SS = (av_nlwrG_df_vg42016SS.sum(axis=0)[1:])/1000 #kw/m2
#sumNlwrG_df_vg42017SS = (av_nlwrG_df_vg42017SS.sum(axis=0)[1:])/1000
#
##%% scalarGroundAbsorbedSolar
#av_nswrG_df_vg4 = readNcfdatasetF0rEachVariableT0dataframe(av_ncfiles_vg4,'scalarGroundAbsorbedSolar',hru_names_df1,DateSc)
#
#av_nswrG_df_vg42016 = av_nswrG_df_vg4[0:8785]
#av_nswrG_df_vg42017 = av_nswrG_df_vg4[8784:]
#
#av_nswrG_df_vg42016SS = av_nswrG_df_vg42016[0:6217]
#av_nswrG_df_vg42017SS = av_nswrG_df_vg42017[0:6193]
#
#sumNswrG_df_vg42016SS = (av_nswrG_df_vg42016SS.sum(axis=0)[1:])/1000 #kw/m2
#sumNswrG_df_vg42017SS = (av_nswrG_df_vg42017SS.sum(axis=0)[1:])/1000
#
##%%  boxplot delta dayofsnowDiss_obs_df
#d1 = [(dayofsnowDiss_vg1_df['2016']-dayofsnowDiss0penT4_df['2016'][1])/24,(dayofsnowDiss_vg2_df['2016']-dayofsnowDiss0penT4_df['2016'][1])/24,
#      (dayofsnowDiss_vg3_df['2016']-dayofsnowDiss0penT4_df['2016'][1])/24,(dayofsnowDiss_vg4_df['2016']-dayofsnowDiss0penT4_df['2016'][1])/24,
#      (dayofsnowDiss_vg1_df['2017']-dayofsnowDiss0penT4_df['2017'][1])/24,(dayofsnowDiss_vg2_df['2017']-dayofsnowDiss0penT4_df['2017'][1])/24,
#      (dayofsnowDiss_vg3_df['2017']-dayofsnowDiss0penT4_df['2017'][1])/24,(dayofsnowDiss_vg4_df['2017']-dayofsnowDiss0penT4_df['2017'][1])/24]
#fig, ax = plt.subplots(1,1, figsize=(20,15))
#bp1 = ax.boxplot(d1, patch_artist=True)
#bp1['boxes'][0].set(color='orange', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp1['boxes'][1].set(color='pink', linewidth=2, facecolor = 'olive', hatch = '/')
#bp1['boxes'][2].set(color='red', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][3].set(color='darkred', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][4].set(color='navy', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][5].set(color='blue', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][6].set(color='skyblue', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][7].set(color='turquoise', linewidth=2, facecolor = 'pink', hatch = '/')
#
#plt.title('T4: dry (2016) ------------------------- T4: wet (2017)',fontsize=40)
#plt.xticks([1,2,3,4,5,6,7,8], ['tall_dense','tall_sperse','short_dense','short_sperse','tall_dense','tall_sperse','short_dense','short_sperse'], fontsize=30, rotation=25)#'open2016','open2017',
#plt.yticks(fontsize=30)
#plt.xlabel('Scenarios', fontsize=30)
#plt.ylabel('delta ssd', fontsize=30)
#plt.savefig('delta_dosd0penVeg1234T4.png')
#
##%%  boxplot  scalarGroundNetNrgFlux
#
##%%  boxplot  scalarGroundNetNrgFlux
#a=sumNswrG_df_vg12016SS-sumNswrG_df_open2016SS[1]
#b=sumNswrG_df_vg22016SS-sumNswrG_df_open2016SS[1]
#c=sumNswrG_df_vg32016SS-sumNswrG_df_open2016SS[1]
#d=sumNswrG_df_vg42016SS-sumNswrG_df_open2016SS[1]
#e=sumNlwrG_df_vg12016SS-sumNlwrG_df_open2016SS[1]
#f=sumNlwrG_df_vg22016SS-sumNlwrG_df_open2016SS[1]
#g=sumNlwrG_df_vg32016SS-sumNlwrG_df_open2016SS[1]
#h=sumNlwrG_df_vg42016SS-sumNlwrG_df_open2016SS[1]
#
#l=sumNswrG_df_vg12017SS-sumNswrG_df_open2017SS[1]
#m=sumNlwrG_df_vg22017SS-sumNswrG_df_open2017SS[1]
#n=sumNswrG_df_vg32017SS-sumNswrG_df_open2017SS[1]
#p=sumNswrG_df_vg42017SS-sumNswrG_df_open2017SS[1]
#q=sumNlwrG_df_vg12017SS-sumNlwrG_df_open2017SS[1]
#r=sumNswrG_df_vg22017SS-sumNlwrG_df_open2017SS[1]
#s=sumNlwrG_df_vg32017SS-sumNlwrG_df_open2017SS[1]
#t=sumNlwrG_df_vg42017SS-sumNlwrG_df_open2017SS[1]
#
#d2 = [a,b,c,d,e,f,g,h,l,m,n,p,q,r,s,t]
#fig2, ax2 = plt.subplots(1,1, figsize=(20,15))
#bp2 = ax2.boxplot(d2, patch_artist=True)
##bp2['boxes'][0].set(color='darkred', linewidth=2, facecolor = 'olive', hatch = '/')
##bp2['boxes'][1].set(color='indianred', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][2].set(color='red', linewidth=2, facecolor = 'olive', hatch = '/')
##bp2['boxes'][3].set(color='palevioletred', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][4].set(color='navy', linewidth=2, facecolor = 'skyblue', hatch = '/')
##bp2['boxes'][5].set(color='darkblue', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][6].set(color='blue', linewidth=2, facecolor = 'skyblue', hatch = '/')
##bp2['boxes'][7].set(color='lightblue', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][8].set(color='darkred', linewidth=2, facecolor = 'olive', hatch = '/')
##bp2['boxes'][9].set(color='indianred', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][10].set(color='red', linewidth=2, facecolor = 'olive', hatch = '/')
##bp2['boxes'][11].set(color='palevioletred', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][12].set(color='navy', linewidth=2, facecolor = 'skyblue', hatch = '/')
##bp2['boxes'][13].set(color='darkblue', linewidth=2, facecolor = 'pink', hatch = '/')
##bp2['boxes'][14].set(color='blue', linewidth=2, facecolor = 'skyblue', hatch = '/')
##bp2['boxes'][15].set(color='lightblue', linewidth=2, facecolor = 'pink', hatch = '/')
#
#
#plt.title('T4: dry (2016) ------------------------- T4: wet (2017)',fontsize=40)
#plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], ['swr_TD','swr_TS','swr_SD','swr_SS',
#                                                   'lwr_TD','lwr_TS','lwr_SD','lwr_SS',
#                                                   'swr_TD','swr_TS','swr_SD','swr_SS',
#                                                   'lwr_TD','lwr_TS','lwr_SD','lwr_SS'], 
#                                                    fontsize=25, rotation=35)#'2016veg1','2017veg1',
#plt.yticks(fontsize=25)
#plt.xlabel('Scenarios', fontsize=35)
#plt.ylabel('delta radiations on the ground (kw/m2)', fontsize=30)
#plt.savefig('deltaR0penVeg1234T42.png')
##%%
#DateSc_2 = [i.strftime("%Y-%m") for i in tvalueSc]
#
#sax = np.arange(0,np.size(DateSc_2))
#sa_xticks = DateSc_2
#safig, saax = plt.subplots(1,1, figsize=(20,15))
#plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
#saax.xaxis.set_major_locator(ticker.AutoLocator())
#plt.yticks(fontsize=20)
#for hru in hru_names_df0[0]:
#    plt.plot(av_sweT4_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#
#plt.plot(swe_obs_df, 'k', markersize=10)
##
#plt.title('SWE (mm)', position=(0.04, 0.88), ha='left', fontsize=40)
#plt.xlabel('Time 2015-2017', fontsize=30)
#plt.ylabel('SWE(mm)', fontsize=30)
#plt.legend()
##plt.show()
#plt.savefig('swe_ScT4_open.png')





















