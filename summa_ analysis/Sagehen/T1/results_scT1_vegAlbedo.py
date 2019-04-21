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

from resultsScT1NcFiles import scT1NcFilesVeg1
#from resultsScT4NcFiles import scT4NcFilesVeg2
#from resultsScT4NcFiles import scT4NcFilesVeg3
#from resultsScT4NcFiles import scT4NcFilesVeg4


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

def calculateDay0fSnowDissappearance(swe_df,hruname_df,firstYear,secondyear):
    av_swe_df2000 = swe_df[:][2000:8784]
    av_swe_df13000 = swe_df[:][13000:17137]
    
    zerosnowdate = []
    for val in hruname_df[0]:
        zerosnowdate.append(np.where(av_swe_df2000[val]==0))
        zerosnowdate.append(np.where(av_swe_df13000[val]==0))

    zerosnowdate_omg = [item[0] for item in zerosnowdate] #change tuple to array

    zerosnowdate_omg2 = [[]] * np.size(zerosnowdate_omg)
    for ijk in range(len(zerosnowdate_omg)/2):
        if np.size(zerosnowdate_omg[2*ijk]) == 0:
            zerosnowdate_omg2[2*ijk].append(np.array([8783]))
        else: zerosnowdate_omg2[2*ijk].append(zerosnowdate_omg[2*ijk]+2000)
    
        if np.size(zerosnowdate_omg[2*ijk+1]) == 0:
            zerosnowdate_omg2[2*ijk+1].append(np.array([17137]))
        else: zerosnowdate_omg2[2*ijk+1].append(zerosnowdate_omg[2*ijk+1]+13000)

    zerosnowdate_omg3 = zerosnowdate_omg2[0]

    dayofsnowDiss = []
    for dosd in range (len(zerosnowdate_omg3)/2):
        dayofsnowDiss.append([zerosnowdate_omg3[2*dosd][0],zerosnowdate_omg3[2*dosd+1][0]])
    
    dayofsnowDiss_df = pd.DataFrame(np.array(dayofsnowDiss))
    dayofsnowDiss_df.columns = [firstYear,secondyear]
    
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

def calculateResidualDOSD(FY,SY,dosdFY_obs,dosdSY_obs,dosd_df,hruNamesDF):
    dayofsnowDiss_obs = np.array([np.array([dosdFY_obs]),np.array([dosdSY_obs])]).T
    dayofsnowDiss_obs_df = pd.DataFrame(dayofsnowDiss_obs,columns=[FY,SY])
    dofs_residual=[]
    for year in dayofsnowDiss_obs_df.columns:
        dofs_residual.append((dayofsnowDiss_obs_df[year][0]-dosd_df[year])/24)
    dosd_residual_df = pd.DataFrame(np.reshape(np.array(dofs_residual),(2,len(dosd_df))).T, columns=['dosdFY','dosdSY'])
    dosd_residual_df.set_index(hruNamesDF[0],inplace=True)
    return dosd_residual_df

def generateDFNetRad(av_nR_df,FY,SY):
    av_nr_df_FY = av_nR_df[0:8785]
    av_nr_df_SY = av_nR_df[8784:]
    sumNR_df_FY = pd.DataFrame((av_nr_df_FY.sum(axis=0)[1:])/1000,columns=[FY])
    sumNR_df_SY = pd.DataFrame((av_nr_df_SY.sum(axis=0)[1:])/1000,columns=[SY])
    sumNR_df = pd.concat([sumNR_df_FY,sumNR_df_SY],axis=1)
    return sumNR_df
#%% SWE observed data T1
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


#%% albedo testing params
p1 = [0.5] #[0.5,0.45,0.5,0.45] #LAIMIN
p2 = [6] #[6,5,6,5] #LAIMAX
p3 = [1] #[1,0.5,1,0.5] #winterSAI
p4 = [5] #[5,2,5,2] #summerLAI
p5 = [25] #[25,25,10,10] #heightCanopyTop
p6 = [5] #[5,5,2,2] #heightCanopyBottom 20% of top
p7 = [40] #[40,27,27,20] #maxMassVegetation         |      25.0000 |       1.0000 |      50.0000

#p1 = [273.7] # tempCritRain
#p6 = [0.001] #zsnow
#p10 = [1100] #specificHeatVeg
#p12 = [0.5] #refInterceptCapRain

#sllj_rsu_veg1_real21112222211.0
#sllj_rsu_veg1_real21112222211.0

p8 = [2,3.5] # 4 mw_exp exponent for meltwater flow

p9 = [0.7,0.78] #0.80 albedoMax |       0.8500 |       0.7000 |       0.9500
p10 = [0.4,0.57] #albedoMinWinter     0.6000 |       0.6000 |       1.0000
p11 = [0.3,0.43] #albedoMinSpring    0.4500 |       0.3000 |       1.0000
p12 = [250000,500000] #700000 albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p13 = [0.7,0.78] #0.80 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p14 = [0.55,0.7] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500
p15 = [0.5,0.57] # 0.6 albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p16 = [0.2,0.38] # 0.4 albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p17 = [3.5,5] # 3 albedoRefresh |       1.0000 |       1.0000 |      10.0000
p18 = [0.3,0.5] # 0.5 albedoSootLoad |       0.3000 |       0.1000 |       0.5000

p19 = [9] #[6,13]refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)

p20 = [0.4] # [0.3,0.7]0.5 throughfallScaleSnow [0.3,0.45,0.6]
p21 = [0.5] #[0.4,0.8]throughfallScaleRain      |       0.6000 |       0.1000 |       0.9000

p22 = [0.35] #[0.1,0.6]ratioDrip2Unloading       |       0.4000 |       0.0000 |       1.0000
p23 = [0] #[0,0.0000014]snowUnloadingCoeff     [0,0.000001,0.0000014]   |       0.0000 |       0.0000 |       1.5d-6

p24 = [0.04] #leafDimension            0.0400 |       0.0100 |       0.1000
p25 = [0.01] #leafExchangeCoeff         |       0.0100 |       0.0010 |       0.1000
p26 = [0.3] #windReductionParam        |       0.2800 |       0.0000 |       1.0000
p27 = [3] #rootingDepth

def hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11):#
    ix1 = np.arange(1,len(p1)+1)
    ix2 = np.arange(1,len(p2)+1)
    ix3 = np.arange(1,len(p3)+1)
    ix4 = np.arange(1,len(p4)+1)
    ix5 = np.arange(1,len(p5)+1)
    ix6 = np.arange(1,len(p6)+1)
    ix7 = np.arange(1,len(p7)+1)
    ix8 = np.arange(1,len(p8)+1)
    ix9 = np.arange(1,len(p9)+1)
    ix10 = np.arange(1,len(p10)+1)
    ix11 = np.arange(1,len(p11)+1)
    
    c = list(itertools.product(ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9,ix10,ix11))#
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxIDI = hru_ix_ID(p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18)#
hruidxID = np.arange(100001,100002)
hru_num = np.size(hruidxID)

out_names = ['sllj_rsu_veg1_real']#'ellj_rsu','eslj_rsu','lllj_rsu','cdc','cdn','csc','csn','cdb','cdu','csb','csu','rdb','rdc','rdn','rdu','rsb','rsc','rsn','rsu'
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df = pd.DataFrame (hru_names1)

#%% reading output_swe files for veg1 scenario
Date, tvalue = dateTime(scT1NcFilesVeg1)
datePlot = [i.strftime("%Y-%m") for i in tvalue]

av_swe_dfVeg1 = readNcfdatasetF0rEachVariableT0dataframe(scT1NcFilesVeg1,'scalarSWE',hru_names_df,Date)
#av_swe_dfVeg2 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg2,'scalarSWE',hru_names_df,Date)
#av_swe_dfVeg3 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg3,'scalarSWE',hru_names_df,Date)
#av_swe_dfVeg4 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg4,'scalarSWE',hru_names_df,Date)

#%% plotting
sax = np.arange(0,np.size(Date))
sa_xticks = datePlot
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20) # 
#saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=20)
for hru in hru_names_df[0]:
    plt.plot(av_swe_dfVeg1[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#
plt.plot(swe_obs_df, 'k', markersize=10)
#
plt.title('SWE (mm); SagehenT1_veg1_albedoSame', position=(0.04, 0.88), ha='left', fontsize=30)
plt.xlabel('Time 2015-2017', fontsize=30)
plt.ylabel('SWE(mm)', fontsize=30)
#plt.legend()
#plt.show()
plt.savefig('swe_ScT1_veg4_albedo.png')

#%%   calculating day of snow disapperance open scenario
dosd_veg1_df=calculateDay0fSnowDissappearance(av_swe_dfVeg1,hru_names_df,'2016','2017')
#dosd_veg2_df=calculateDay0fSnowDissappearance(av_swe_dfVeg2,hru_names_df,'2016','2017')
#dosd_veg3_df=calculateDay0fSnowDissappearance(av_swe_dfVeg3,hru_names_df,'2016','2017')
#dosd_veg4_df=calculateDay0fSnowDissappearance(av_swe_dfVeg4,hru_names_df,'2016','2017')

dosd_residual_Veg1_df = calculateResidualDOSD('2016','2017',4656,14424,dosd_veg1_df,hru_names_df)
#dosd_residual_Veg2_df = calculateResidualDOSD('2016','2017',5870,15155,dosd_veg2_df)
#dosd_residual_Veg3_df = calculateResidualDOSD('2016','2017',5870,15155,dosd_veg3_df)
#dosd_residual_Veg4_df = calculateResidualDOSD('2016','2017',5870,15155,dosd_veg4_df)
dosd_residual_Veg1_min = dosd_residual_Veg1_df['dosdFY'][:].min()
pareto_dosd_veg1 = dosd_residual_Veg1_df['dosdFY'][:].idxmin()

dosd_resisual_veg1_min10p = dosd_residual_Veg1_df.nsmallest(200, 'dosdFY')
pareto_dosd_veg1_min10p = dosd_resisual_veg1_min10p.index.tolist()
#%%

#dosd_resisual_veg1_min15p = dosd_residual_Veg1_df.nlargest(300, 'dosdSY')
#pareto_dosd_veg1_min15p = dosd_resisual_veg1_min15p.index.tolist()

#maxSwe10p = []
#for hru in pareto_dosd_veg1_min10p:
#    maxSwe10p.append(av_swe_dfVeg1[hru][0:8500].max())
#    
sax = np.arange(0,np.size(Date))
sa_xticks = datePlot
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=20)
for hru in pareto_dosd_veg1_min10p:
    plt.plot(av_swe_dfVeg1[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
    
plt.plot(swe_obs_df, 'k', markersize=10)
#
plt.title('SWE (mm); SagehenT1_veg1_sameAlbedo', position=(0.04, 0.88), ha='left', fontsize=25)
plt.xlabel('Time 2015-2017', fontsize=30)
plt.ylabel('SWE(mm)', fontsize=30)
#plt.legend()
plt.show()
plt.savefig('swe_ScT1_veg1_samealbedo.png')

#%% scalarLWNetGround
#av_nlwrG_df_vg1 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg1,'scalarLWNetGround',hru_names_df,Date)
#av_nlwrG_df_vg2 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg2,'scalarLWNetGround',hru_names_df,Date)
#av_nlwrG_df_vg3 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg3,'scalarLWNetGround',hru_names_df,Date)
#av_nlwrG_df_vg4 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg4,'scalarLWNetGround',hru_names_df,Date)

#sumNlwrG_df_vg1 = generateDFNetRad(av_nlwrG_df_vg1,'2016','2017')
#sumNlwrG_df_vg2 = generateDFNetRad(av_nlwrG_df_vg2,'2016','2017')
#sumNlwrG_df_vg3 = generateDFNetRad(av_nlwrG_df_vg3,'2016','2017')
#sumNlwrG_df_vg4 = generateDFNetRad(av_nlwrG_df_vg4,'2016','2017')

#%% scalarGroundAbsorbedSolar
#av_nswrG_df_vg1 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg1,'scalarGroundAbsorbedSolar',hru_names_df,Date)
#av_nswrG_df_vg2 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg2,'scalarGroundAbsorbedSolar',hru_names_df,Date)
#av_nswrG_df_vg3 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg3,'scalarGroundAbsorbedSolar',hru_names_df,Date)
#av_nswrG_df_vg4 = readNcfdatasetF0rEachVariableT0dataframe(scT4NcFilesVeg4,'scalarGroundAbsorbedSolar',hru_names_df,Date)

#sumNswrG_df_vg1 = generateDFNetRad(av_nswrG_df_vg1,'2016','2017')
#sumNswrG_df_vg2 = generateDFNetRad(av_nswrG_df_vg2,'2016','2017')
#sumNswrG_df_vg3 = generateDFNetRad(av_nswrG_df_vg3,'2016','2017')
#sumNswrG_df_vg4 = generateDFNetRad(av_nswrG_df_vg4,'2016','2017')
#%%  boxplot delta dayofsnowDiss_obs_df
#d1 = [dosd_residual_Veg1_df['dosdFY'],dosd_residual_Veg1_df['dosdSY'],
#      dosd_residual_Veg2_df['dosdFY'],dosd_residual_Veg2_df['dosdSY'],
#      dosd_residual_Veg3_df['dosdFY'],dosd_residual_Veg3_df['dosdSY'],
#      dosd_residual_Veg4_df['dosdFY'],dosd_residual_Veg4_df['dosdSY']
#      ]

#fig, ax = plt.subplots(1,1, figsize=(20,15))
#bp1 = ax.boxplot(d1, patch_artist=True)
#bp1['boxes'][0].set(color='blue', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp1['boxes'][1].set(color='red', linewidth=2, facecolor = 'olive', hatch = '/')
#bp1['boxes'][2].set(color='navy', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][3].set(color='darkred', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][4].set(color='skyblue', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][5].set(color='orange', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][6].set(color='turquoise', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][7].set(color='pink', linewidth=2, facecolor = 'pink', hatch = '/')

#plt.title('T4: dry (2016) -----**************----- T4: wet (2017)',fontsize=40)
#['tall_dense','tall_sperse','short_dense','short_sperse','tall_dense','tall_sperse','short_dense','short_sperse']
#plt.xticks([1,2,3,4,5,6,7,8], ['2016veg1','2017veg1','2016veg2','2017veg2','2016veg3','2017veg3','2016veg4','2017veg4'], fontsize=30, rotation=25)#'open2016','open2017',
#plt.yticks(fontsize=30)
#plt.xlabel('vegSc', fontsize=30)
#plt.ylabel('residual day of snow dissappearance', fontsize=30)
#plt.savefig('rdosd_VegscT4.png')







