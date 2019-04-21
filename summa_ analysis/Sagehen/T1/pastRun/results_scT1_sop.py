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
p10 = [3]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p11 = [100000]#,200000,400000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 

p12 = [0.65] #albedoMinWinter           |       0.6500 |       0.6000 |       1.0000
p13 = [0.5] #albedoMinSpring           |       0.5000 |       0.3000 |       1.0000

p14 = [60] #newSnowDenMin 
p15= [75] #newSnowDenMult            |      75.0000 |      25.0000 |      75.0000

p16 = [0.01] #winterSAI
p17 = [0.1] #summerLAI
#p1 = [0.1] #LAIMIN
#p2 = [1] #LAIMAX
p18 = [0.3] #heightCanopyTop
p19 = [0.03] #heightCanopyBottom to calculate wind speed, wind speed reduction; coupute roughness length of the veg canopy; neutral ground resistance; canopy air temp;

p20 = [0.01] #0.1#z0Canopy                  |       0.1000

list_param = pd.DataFrame([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20])

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
#
hru_num = np.size(hruidxID)
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

#max_swe_obs = max(obs_swe['swe_mm'])
#max_swe_date_obs = obs_swe[obs_swe ['swe_mm']== max_swe_obs].index.tolist()    
#%%
out_names = ['pt11', 'pt12', 'pt13', 'pt21', 'pt22', 'pt23']
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df = pd.DataFrame (hru_names1)
#%% reading output_swe files
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

for varname in av_all[0].variables.keys():
    var = av_all[0].variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)

av_swe = []
for dfs in av_all:
    av_swe.append(pd.DataFrame(dfs['scalarSWE'][:]))

av_swe_2yearcons = []
for dfs2 in range (len(av_swe)/2):
    av_swe_2yearcons.append(pd.concat([av_swe[2*dfs2],av_swe[2*dfs2+1]], ignore_index=True))
av_swe_df = pd.concat (av_swe_2yearcons, axis=1)
av_swe_df.columns = hru_names_df[0]

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
#%% day of snow disappearance-final output
av_swe_df.set_index(pd.DatetimeIndex(DateSc),inplace=True)
counter = pd.DataFrame(np.arange(0,np.size(DateSc)),columns=['counter'])
counter.set_index(av_swe_df.index,inplace=True)
av_swe_df2 = pd.concat([counter, av_swe_df], axis=1)

#%% plotting
#DateSc_2 = [i.strftime("%Y-%m") for i in tvalueSc]
#sax = np.arange(0,np.size(DateSc))
#sa_xticks = DateSc_2
safig, saax = plt.subplots(1,1, figsize=(20,15))
#plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
#saax.xaxis.set_major_locator(ticker.AutoLocator())
#plt.yticks(fontsize=20)
for hru in av_swe_df.columns:
    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]

plt.plot(swe_obs_df, 'k', markersize=10)
#
#plt.title('SWE (mm)', position=(0.04, 0.88), ha='left', fontsize=40)
#plt.xlabel('Time 2015-2017', fontsize=30)
#plt.ylabel('SWE(mm)', fontsize=30)
##plt.legend()
##plt.show()
#plt.savefig('initialResults_sc\sweScT1_TP.png')
#%%   calculating day of snow disapperance
av_swe_df4000 = av_swe_df2[:][4000:8784]
av_swe_df13000 = av_swe_df2[:][13000:17137]

zerosnowdate = []
for val in hru_names_df[0]:
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

dayofsnowDiss_obs = np.array([np.array([4686]),np.array([14430])]).T
dayofsnowDiss_obs_df = pd.DataFrame(dayofsnowDiss_obs,columns=['2016','2017'])

dosd_residual=[]
for years in dayofsnowDiss_df.columns:
    dosd_residual.append(abs(dayofsnowDiss_obs_df[years][0]-dayofsnowDiss_df[years])/24)

dosd_residual_df = pd.DataFrame(np.reshape(np.array(dosd_residual),(2,3072)).T, columns=['dosd2016','dosd2017'])

#%%
#plt.xticks(x, hru[::3], rotation=25)
#for namefile in out_names:
#    x = list(np.arange(1,244))
#    fig = plt.figure(figsize=(20,15))
#    plt.bar(x,zerosnowdate_residual_df[namefile])
#    plt.title(namefile, fontsize=42)
#    plt.xlabel('hrus',fontsize=30)
#    plt.ylabel('residual dosd (day)', fontsize=30)
#    #vax.yaxis.set_label_coords(0.5, -0.1) 
#    plt.savefig('SA2/'+namefile)


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

#plt.subplots(1,1, figsize=(20,15))
#for hru in av_swe_df.columns:
#    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#plt.plot(swe_obs_df, 'k', markersize=10)
#
#plt.plot(maxSwe_df['2016_1'])
#plt.plot(maxSwe_obs_df['2016_1'], 'ko')
#%% find best params
criteria_df = pd.concat([maxSwe_residual_df, dosd_residual_df], axis=1) 
criteria_df.set_index(hru_names_df[0],inplace=True)
pareto_model_param = pd.DataFrame(criteria_df.index[((criteria_df['maxSwe2016_1']) <= 50) & ((criteria_df['maxSwe2016_2'])<=75) & 
                                                    ((criteria_df['maxSwe2017_1'])<=150) & ((criteria_df['maxSwe2017_2'])<=150) &
                                                    ((criteria_df['dosd2016'])<=10) & ((criteria_df['dosd2017'])<=9)  ].tolist())#

#DateSc_2 = [i.strftime("%Y-%m") for i in tvalueSc]
#sax = np.arange(0,np.size(DateSc))
#sa_xticks = DateSc_2
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
#plt.savefig('initialResults_sc\sweScT1_TPBP.png')
#%%    dayofsnowDiss_obs_df
d1 = [dayofsnowDiss_df['2016'],dayofsnowDiss_df['2017']-8784]
fig = plt.subplots(1,1, figsize=(20,15))
bp1 = plt.boxplot(d1, patch_artist=True)
bp1['boxes'][0].set(color='navy', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][1].set(color='blue', linewidth=2, facecolor = 'olive', hatch = '/')
#bp1['boxes'][2].set(color='skyblue', linewidth=2, facecolor = 'pink', hatch = '/')

plt.xticks([1,2], ['open2016','open2017'], fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Scenarios', fontsize=30)
plt.ylabel('day of snow dissappearance', fontsize=30)
plt.savefig('open.png')
#
#d0 = [desiredMaxSweT0,desiredMaxSweT2,desiredMaxSweT4]
#d1 = [desiredresDateT0,desiredresDateT2,desiredresDateT4]
#
#bp0 = plt.boxplot(d0, patch_artist=True)
#bp1 = plt.boxplot(d1, patch_artist=True)
#
#bp0['boxes'][0].set(color='red', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp0['boxes'][1].set(color='orange', linewidth=2, facecolor = 'olive', hatch = '/')
#bp0['boxes'][2].set(color='tan', linewidth=2, facecolor = 'pink', hatch = '/')
#
##plt.hold()
#
#bp1['boxes'][0].set(color='navy', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp1['boxes'][1].set(color='blue', linewidth=2, facecolor = 'olive', hatch = '/')
#bp1['boxes'][2].set(color='skyblue', linewidth=2, facecolor = 'pink', hatch = '/')
#
#plt.xticks([1, 2, 3], ['T0', 'T+2', 'T+4'])
#plt.savefig('resSwe2.png')
























