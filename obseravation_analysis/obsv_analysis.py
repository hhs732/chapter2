import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
import matplotlib.ticker as ticker
import csv
import scipy.io as sio
from datetime import datetime
import matplotlib.ticker as plticker

def calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(snwdpth_df):
    FDSS = []
    LDSG = []
    for col in snwdpth_df.columns[1:]:
        if (len (snwdpth_df.loc[snwdpth_df[col]==0]))>0:
            FDSS.append((snwdpth_df.loc[snwdpth_df[col]==0].index[0])/24)
        else:FDSS.append(365)
    
        if (len (snwdpth_df.loc[snwdpth_df[col]!=0]))>0:
            LDSG.append((snwdpth_df.loc[snwdpth_df[col]!=0].index[-1])/24)
        else:LDSG.append(365)
    return FDSS,LDSG

def calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(snwdpth_df):
    FDSS = []
    LDSG = []
    for col in snwdpth_df.columns:
        if (len (snwdpth_df.loc[snwdpth_df[col]==0]))>0:
            FDSS.append((snwdpth_df.loc[snwdpth_df[col]==0].index[0])/24)#
        else:FDSS.append(200)
    
        if (len (snwdpth_df.loc[snwdpth_df[col]!=0]))>0:
            LDSG.append((snwdpth_df.loc[snwdpth_df[col]!=0].index[-1])/24)#
        else:LDSG.append(290)
    return FDSS,LDSG

def calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundDayBased(snwdpth_df):
    FDSS = []
    LDSG = []
    for col in snwdpth_df.columns:
        if (len (snwdpth_df.loc[snwdpth_df[col]==0]))>0:
            FDSS.append((snwdpth_df.loc[snwdpth_df[col]==0].index[0]))#
        else:FDSS.append(200)
    
        if (len (snwdpth_df.loc[snwdpth_df[col]!=0]))>0:
            LDSG.append((snwdpth_df.loc[snwdpth_df[col]!=0].index[-1]))#
        else:LDSG.append(290)
    return FDSS,LDSG
#%% Sagehen Creek
#read snow is present when temperature is between 1 and -1 C and when the daily standard deviation in temperature is less than or equal to .353 C
with open('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/tree_sagehen.csv') as scvdT:
    readerT = csv.reader(scvdT)
    raw_BinDSS0T = [rT for rT in readerT]
raw_BinDSST = raw_BinDSS0T[1:]

raw_BinDSS_colT = []
for csv_counter1T in range (len (raw_BinDSST)):
    for csv_counter2T in range (len(rT)):
        raw_BinDSS_colT.append(float(raw_BinDSST[csv_counter1T][csv_counter2T]))
binDSST=np.reshape(raw_BinDSS_colT,(len (raw_BinDSST),len(rT))) 
binDSST_df = pd.DataFrame(binDSST)
binDSST_df.set_index(binDSST_df[0], inplace=True)
binDSST_df.drop(binDSST_df.columns[[0]], axis=1, inplace=True)
#%%
FDSS_lsT = []
LDSS_lsT = []
for col in range (1,len(binDSST_df.columns)):
    FDSS_lsT.append(binDSST_df.loc[binDSST_df[col]==0].index[0])
    if (len (binDSST_df.loc[binDSST_df[col]==1]))>0:
        LDSS_lsT.append((binDSST_df.loc[binDSST_df[col]==1].index[-1])+1)
    else:LDSS_lsT.append(70)
#%%
with open("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/open_sagehen.csv") as scvd0:
    reader0 = csv.reader(scvd0)
    raw_BinDSS00 = [r0 for r0 in reader0]
raw_BinDSS0 = raw_BinDSS00[1:]

raw_BinDSS_col0 = []
for csv_counter10 in range (len (raw_BinDSS0)):
    for csv_counter20 in range (len(r0)):
        raw_BinDSS_col0.append(float(raw_BinDSS0[csv_counter10][csv_counter20]))
binDSS0=np.reshape(raw_BinDSS_col0,(len (raw_BinDSS0),len(r0))) 
binDSS0_df = pd.DataFrame(binDSS0)
binDSS0_df.set_index(binDSS0_df[0], inplace=True)
binDSS0_df.drop(binDSS0_df.columns[[0]], axis=1, inplace=True)
#%%
FDSS_ls0 = []
LDSS_ls0 = []
for col in range (1,len(binDSS0_df.columns)):
    if (len (binDSS0_df.loc[binDSS0_df[col]==0]))>0:
        FDSS_ls0.append(binDSS0_df.loc[binDSS0_df[col]==0].index[0])
    else:FDSS_ls0.append(70)
    if (len (binDSS0_df.loc[binDSS0_df[col]==1]))>0:
        LDSS_ls0.append((binDSS0_df.loc[binDSS0_df[col]==1].index[-1])+1)
    else:LDSS_ls0.append(70)
#%%
d1 = [FDSS_ls0,FDSS_lsT]#LDSS_ls0,,LDSS_lsT

fig, ax = plt.subplots(1,1, figsize=(30,20))
bp1 = ax.boxplot(d1, patch_artist=True)
bp1['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][1].set(color='green', linewidth=2, facecolor = 'olive', hatch = '/')
#bp1['boxes'][2].set(color='khaki', linewidth=2, facecolor = 'pink', hatch = '/')
#bp1['boxes'][3].set(color='seagreen', linewidth=2, facecolor = 'pink', hatch = '/')

plt.title('Sagehen Creek Watershed (SCW), CA',fontsize=50)

plt.xticks([1,2], ['0pen site','under Tree'], fontsize=40, rotation=25)#3,4,'underTree_LSG','0pen_LSG','open2016','open2017',
plt.yticks(fontsize=40)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('Snow disappearance day - Julian days', fontsize=40)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/dtsDta_snowPresence_SC3.png')

#%% Jemez_VCM
snwDpth_vcm = sio.loadmat('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Jemez/valles_dep_0810.mat')
sd_vcm09_ut = pd.DataFrame(np.array([snwDpth_vcm['dep07'][323:,0],snwDpth_vcm['dep07'][323:,2],
                                     snwDpth_vcm['dep07'][323:,8]]).T,columns=['date','ut1','ut2'])
sd_vcm09_op = pd.DataFrame(np.array([snwDpth_vcm['dep07'][323:,0],snwDpth_vcm['dep07'][323:,1],
                                     snwDpth_vcm['dep07'][323:,4],snwDpth_vcm['dep07'][323:,6]]).T,
                                     columns=['date','op1','op2','op3'])
FDSS09_vcm_ut , LDSG09_vcm_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_vcm09_ut)
FDSS09_vcm_op , LDSG09_vcm_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_vcm09_op)

sd_vcm10_ut = pd.DataFrame(np.array([snwDpth_vcm['dep08'][914:,0],snwDpth_vcm['dep08'][914:,2],
                                     snwDpth_vcm['dep08'][914:,7],snwDpth_vcm['dep08'][914:,8]]).T,
                                     columns=['date','ut1','ut4','ut3'])
sd_vcm10_op = pd.DataFrame(np.array([snwDpth_vcm['dep08'][914:,0],snwDpth_vcm['dep08'][914:,1],
                                     snwDpth_vcm['dep08'][914:,4],snwDpth_vcm['dep08'][914:,6]]).T,
                                     columns=['date','op1','op2','op3'])
FDSS10_vcm_ut , LDSG10_vcm_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_vcm10_ut)
FDSS10_vcm_op , LDSG10_vcm_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_vcm10_op)

sd_vcm11_ut = pd.DataFrame(np.array([snwDpth_vcm['dep09'][1054:,0],snwDpth_vcm['dep09'][1054:,2],
                                     snwDpth_vcm['dep09'][1054:,7],snwDpth_vcm['dep09'][1054:,8]]).T,
                                     columns=['date','ut1','ut4','ut3'])
sd_vcm11_op = pd.DataFrame(np.array([snwDpth_vcm['dep09'][1054:,0],snwDpth_vcm['dep09'][1054:,1],
                                     snwDpth_vcm['dep09'][1054:,4],snwDpth_vcm['dep09'][1054:,6]]).T,
                                     columns=['date','op1','op2','op3'])
FDSS11_vcm_ut , LDSG11_vcm_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_vcm11_ut)
FDSS11_vcm_op , LDSG11_vcm_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_vcm11_op)
#%%
sd_jm10_date = pd.DatetimeIndex(sd_vcm10_op['date'])
tvalue_jm10 = num2date(sd_jm10_date)#, units=t_unit, calendar=t_cal)
date_jm10 = [i.strftime("%Y-%m-%d") for i in tvalue_jm10] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")  

sax = np.arange(0,np.size(date_jm10))
sa_xticks = date_jm10
safig, saax = plt.subplots(1,1, figsize=(30,20))
plt.xticks(sax, sa_xticks[::400], rotation=25, fontsize=40) # 
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=30)
plt.plot(sd_vcm10_op['op1'],linewidth=3,label='open_site1',color='maroon')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_vcm10_op['op2'],linewidth=3,label='open_site2',color='maroon')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_vcm10_op['op3'],linewidth=3,label='open_site3',color='maroon')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_vcm10_ut['ut1'],linewidth=3,label='undertree_site1',color='green')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_vcm10_ut['ut3'],linewidth=3,label='undertree_site2',color='darkgreen')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_vcm10_ut['ut4'],linewidth=3,label='undertree_site3',color='seagreen')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]

plt.xlabel('Time 2010', fontsize=40)
plt.ylabel('snow depth (cm)', fontsize=40)
plt.title('JRBN', fontsize=50, loc = 'right', x = 0.95, y=0.9)
plt.legend(fontsize=30, loc = 3)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/documents/jemez_obs20102.png')
#%%
with open('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Jemez/winter04-05_level0_level_1.csv') as vcmT1:
    readerv1 = csv.reader(vcmT1)
    raw_sd_vcm05 = [rv1 for rv1 in readerv1]
sd_vcm05_df = pd.DataFrame(raw_sd_vcm05,columns=['date','ut1','ut2','op1','op2']) 
sd_vcm05_df.index = np.arange(1582,len(sd_vcm05_df)+1582).astype(int)
sd_vcm05_ut = sd_vcm05_df[['ut1','ut2']].astype(float)
sd_vcm05_op = sd_vcm05_df[['op1','op2']].astype(float)

FDSS05_vcm_ut , LDSG05_vcm_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_vcm05_ut)
FDSS05_vcm_op , LDSG05_vcm_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_vcm05_op)
        
with open('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Jemez/winter05-06_level0_level_1.csv') as vcmT6:
    readerv6 = csv.reader(vcmT6)
    raw_sd_vcm06 = [rv6 for rv6 in readerv6]
sd_vcm06_df0 = pd.DataFrame(raw_sd_vcm06,columns=['date','ut1','ut2','ut3','op1','op2','op3'])
sd_vcm06_df = sd_vcm06_df0 [621:]
sd_vcm06_df.index = np.arange(625,len(sd_vcm06_df)+625).astype(int)
sd_vcm06_ut = sd_vcm06_df[['ut1','ut2','ut3']].astype(float)
sd_vcm06_op = sd_vcm06_df[['op1','op2','op3']].astype(float)

FDSS06_vcm_ut , LDSG06_vcm_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_vcm06_ut)
FDSS06_vcm_op , LDSG06_vcm_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_vcm06_op)

with open('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Jemez/winter06-07_level0_level_1.csv') as vcmT7:
    readerv7 = csv.reader(vcmT7)
    raw_sd_vcm07 = [rv7 for rv7 in readerv7]
sd_vcm07_df0 = pd.DataFrame(raw_sd_vcm07,columns=['date','ut1','ut2','ut3','op1','op2','op3'])
sd_vcm07_df = sd_vcm07_df0[1056:]
sd_vcm07_df.index = np.arange(0,len(sd_vcm07_df)+0).astype(int)
sd_vcm07_ut = sd_vcm07_df[['ut1','ut2','ut3']].astype(float)
sd_vcm07_op = sd_vcm07_df[['op1','op2','op3']].astype(float)

FDSS08_vcm_ut , LDSG08_vcm_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_vcm07_ut)
FDSS08_vcm_op , LDSG08_vcm_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_vcm07_op)

#%% ploting
d2=[FDSS05_vcm_op, FDSS05_vcm_ut,FDSS06_vcm_op, FDSS06_vcm_ut,  
    FDSS08_vcm_op, FDSS08_vcm_ut,FDSS09_vcm_op, FDSS09_vcm_ut,
    FDSS10_vcm_op, FDSS10_vcm_ut,FDSS11_vcm_op, FDSS11_vcm_ut]

fig2, ax2 = plt.subplots(1,1, figsize=(30,20))
bp2 = ax2.boxplot(d2, patch_artist=True)
bp2['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][6].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][7].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][8].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][9].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][10].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][11].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')

plt.title('Jemez (JMZ), NM',fontsize=50)

plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12], 
           ['0pen2005','UnTr2005','0pen2006','UnTr2006','0pen2008','UnTr2008',
            '0pen2009','UnTr2009','0pen2010','UnTr2010','0pen2011','UnTr2011'],
            fontsize=40, rotation=50)#'open2016','open2017',
plt.yticks(fontsize=40)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('day of snow disappearance - Julian days', fontsize=40)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/obsData_dsd_jemez2.png')
#%% ploting
d3=[LDSG05_vcm_op, LDSG05_vcm_ut,LDSG06_vcm_op, LDSG06_vcm_ut,  
    LDSG08_vcm_op, LDSG08_vcm_ut,LDSG09_vcm_op, LDSG09_vcm_ut,
    LDSG10_vcm_op, LDSG10_vcm_ut,LDSG11_vcm_op, LDSG11_vcm_ut]

fig3, ax3 = plt.subplots(1,1, figsize=(20,15))
bp3 = ax3.boxplot(d3, patch_artist=True)
bp3['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp3['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp3['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp3['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp3['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp3['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp3['boxes'][6].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp3['boxes'][7].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp3['boxes'][8].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp3['boxes'][9].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp3['boxes'][10].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp3['boxes'][11].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')

plt.title('JMZ, last day snow on the ground',fontsize=35)

plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12], 
           ['0p05_LDSG','uT05_LDSG','0p06_LDSG','uT06_LDSG',
            '0p07_LDSG','uT07_LDSG','0p08_LDSG','uT08_LDSG',
            '0p09_LDSG','uT09_LDSG','0p10_LDSG','uT10_LDSG'],
            fontsize=25, rotation=25)#'open2016','open2017',
plt.yticks(fontsize=30)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('Julian days', fontsize=30)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/documents/obsData_LDSG_jemez.png')
#%% Boulder_Niwot
snwDpth_nwt = sio.loadmat('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Niwot/niwot_dep_0810.mat')
sd_nwt09ut = pd.DataFrame(np.array([snwDpth_nwt['dep07'][1041:,0],snwDpth_nwt['dep07'][1041:,2],
                                    snwDpth_nwt['dep07'][1041:,7],snwDpth_nwt['dep07'][1041:,8]]).T,columns=['date','ut1','ut2','ut3'])
sd_nwt09op = pd.DataFrame(np.array([snwDpth_nwt['dep07'][1041:,0],snwDpth_nwt['dep07'][1041:,1],
                                    snwDpth_nwt['dep07'][1041:,4],snwDpth_nwt['dep07'][1041:,6]]).T,columns=['date','op1','op2','op3'])

FDSS09_nwt_ut , LDSG09_nwt_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_nwt09ut)
FDSS09_nwt_op , LDSG09_nwt_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_nwt09op)

sd_nwt10ut = pd.DataFrame(np.array([snwDpth_nwt['dep08'][1041:,0],snwDpth_nwt['dep08'][1041:,2],
                                    snwDpth_nwt['dep08'][1041:,7],snwDpth_nwt['dep08'][1041:,8]]).T,columns=['date','ut1','ut2','ut3'])
sd_nwt10op = pd.DataFrame(np.array([snwDpth_nwt['dep08'][1041:,0],snwDpth_nwt['dep08'][1041:,1],
                                    snwDpth_nwt['dep08'][1041:,4]]).T,columns=['date','op1','op2'])

FDSS10_nwt_ut , LDSG10_nwt_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_nwt10ut)
FDSS10_nwt_op , LDSG10_nwt_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_nwt10op)

#%%
tvalue_nw10 = num2date(sd_nwt10op['date'])#, units=t_unit, calendar=t_cal)
date_nw10 = [i.strftime("%Y-%m-%d") for i in tvalue_nw10] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")  

sax = np.arange(0,np.size(date_nw10))
sa_xticks = date_nw10
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::400], rotation=25, fontsize=30) # 
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=30)
plt.plot(sd_nwt10op['op1'],linewidth=3,label='open_site1',color='maroon')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_nwt10op['op2'],linewidth=3,label='open_site2',color='maroon')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_nwt10ut['ut1'],linewidth=3,label='undertree_site1',color='green')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_nwt10ut['ut2'],linewidth=3,label='undertree_site2',color='darkgreen')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]
plt.plot(sd_nwt10ut['ut3'],linewidth=3,label='undertree_site3',color='seagreen')#, sbx, swe_obs2006, 'k--', )#, ) param_nam_list[q] color_list[q]

plt.xlabel('Time 2010', fontsize=30)
plt.ylabel('snow depth (cm)', fontsize=30)
plt.title('NRC', fontsize=50, loc = 'right', x = 0.95, y=0.9)
plt.legend(fontsize=20, loc = 2)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/documents/boulder_obs20102.png')
#%%

sd_nwt11ut = pd.DataFrame(np.array([snwDpth_nwt['dep09'][1041:,0],snwDpth_nwt['dep09'][1041:,2],
                                    snwDpth_nwt['dep09'][1041:,7],snwDpth_nwt['dep09'][1041:,8]]).T,columns=['date','ut1','ut2','ut3'])
sd_nwt11op = pd.DataFrame(np.array([snwDpth_nwt['dep09'][1041:,0],snwDpth_nwt['dep09'][1041:,1],
                                    snwDpth_nwt['dep09'][1041:,4]]).T,columns=['date','op1','op2'])

FDSS11_nwt_ut , LDSG11_nwt_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_nwt11ut)
FDSS11_nwt_op , LDSG11_nwt_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased(sd_nwt11op)

with open('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Niwot/GLV_snow_WY05_06_level0_level_1.csv') as nwtT1:
    readern1 = csv.reader(nwtT1)
    raw_sd_nwt06 = [rn1 for rn1 in readern1]
sd_nwt06_df0 = pd.DataFrame(raw_sd_nwt06,columns=['date','ut1','ut2','op1','op2'])
sd_nwt06_df = sd_nwt06_df0[1642:]
sd_nwt06_df.index = np.arange(0,len(sd_nwt06_df)+0).astype(int)
sd_nwt07_ut = sd_nwt06_df[['ut1','ut2']].astype(float)
sd_nwt07_op = sd_nwt06_df[['op1','op2']].astype(float)

FDSS07_nwt_ut , LDSG07_nwt_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_nwt07_ut)
FDSS07_nwt_op , LDSG07_nwt_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundHourBased2(sd_nwt07_op)

with open('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/Niwot/GLV_snow_WY06_07_level0_level_1.csv') as nwtT2:
    readern2 = csv.reader(nwtT2)
    raw_sd_nwt07 = [rn2 for rn2 in readern2]
sd_nwt07_df0 = pd.DataFrame(raw_sd_nwt07,columns=['date','ut1','ut2','op1','op2'])
sd_nwt07_df = sd_nwt07_df0[77:]
sd_nwt07_df.index = np.arange(0,len(sd_nwt07_df)+0).astype(int)
sd_nwt08_ut = sd_nwt07_df[['ut1','ut2']].astype(float)
sd_nwt08_op = sd_nwt07_df[['op1','op2']].astype(float)

FDSS08_nwt_ut , LDSG08_nwt_ut = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundDayBased(sd_nwt08_ut)
FDSS08_nwt_op , LDSG08_nwt_op = calculatingDay0fSnowDisappearanceLastDay0fSnow0nGroundDayBased(sd_nwt08_op)

#%% ploting
d4=[FDSS07_nwt_op, FDSS07_nwt_ut, FDSS08_nwt_op, FDSS08_nwt_ut, FDSS09_nwt_op, FDSS09_nwt_ut, 
    FDSS10_nwt_op, FDSS10_nwt_ut, FDSS11_nwt_op, FDSS11_nwt_ut]

fig4, ax4 = plt.subplots(1,1, figsize=(30,20))
bp4 = ax4.boxplot(d4, patch_artist=True)
bp4['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp4['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp4['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp4['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp4['boxes'][6].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][7].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][8].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp4['boxes'][9].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')

plt.title('Niwot Ridge Creek (NRC), C0',fontsize=50)

plt.xticks([1,2,3,4,5,6,7,8,9,10],['0pen2007','UnTr2007','0pen08','UnTr2008','0pen2009','UnTr2009','0pen2010','UnTr2010','0pen2011','UnTr2011'],
            fontsize=40, rotation=50)#'open2016','open2017',
plt.yticks(fontsize=40)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('Day of snow disappearance - Julian days', fontsize=40)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/obsData_dsd_Niwot2.png')
#%%
#d5=[LDSG06_nwt_op, LDSG06_nwt_ut, LDSG07_nwt_op, LDSG07_nwt_ut, LDSG08_nwt_op, LDSG08_nwt_ut]
#
#fig5, ax5 = plt.subplots(1,1, figsize=(20,15))
#bp5 = ax5.boxplot(d5, patch_artist=True)
#bp5['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp5['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
#bp5['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
#bp5['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
#bp5['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp5['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
#
#plt.title('Niwot, last day snow on the ground (LDSG)',fontsize=35)
#
#plt.xticks([1,2,3,4,5,6],['0p06_LDSG','uT06_LDSG','0p07_LDSG','uT07_LDSG','0p08_LDSG','uT08_LDSG'],
#            fontsize=25, rotation=25)#'open2016','open2017',
#plt.yticks(fontsize=30)
##plt.xlabel('vegSc_year', fontsize=30)
#plt.ylabel('Julian days', fontsize=30)
#plt.savefig('obsData_LDSG_Niwot.png')
#
#%% KREW
ssd_krw_op10 = [149,153,155,157,158,158]
ssd_krw_op11 = [145,169,167,176,169,158]
ssd_krw_op12 = [119,123,124,130,130]

ssd_krw_ut10 = [138,151,158,138,157]
ssd_krw_ut11 = [131,155,173,132,151]
ssd_krw_ut12 = [114,113,124,130,130]

d5=[ssd_krw_op10,ssd_krw_ut10,ssd_krw_op11,ssd_krw_ut11,ssd_krw_op12,ssd_krw_ut12]

fig5, ax5 = plt.subplots(1,1, figsize=(30,20))
bp5 = ax5.boxplot(d5, patch_artist=True)
bp5['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp5['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp5['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp5['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp5['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp5['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')

plt.title('King River Experimental Watershed (KREW), CA',fontsize=50)

plt.xticks([1,2,3,4,5,6],['0pen2010','underTree2010','0pen2011','underTree2011','0pen2012','underTree2012'],
            fontsize=35, rotation=50)#'open2016','open2017'
plt.yticks(fontsize=40)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('day of snow disappearance - Julian days', fontsize=40)
plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/obsData_dsd_krew2.png')
#%%
d11 = [FDSS_ls0,FDSS_lsT]#LDSS_ls0,,LDSS_lsT
d12 = [ssd_krw_op10,ssd_krw_ut10,ssd_krw_op11,ssd_krw_ut11,ssd_krw_op12,ssd_krw_ut12]
d21 = [FDSS05_vcm_op, FDSS05_vcm_ut,FDSS06_vcm_op, FDSS06_vcm_ut,  
       FDSS08_vcm_op, FDSS08_vcm_ut,FDSS09_vcm_op, FDSS09_vcm_ut,
       FDSS10_vcm_op, FDSS10_vcm_ut,FDSS11_vcm_op, FDSS11_vcm_ut]
d22 = [FDSS07_nwt_op, FDSS07_nwt_ut, FDSS08_nwt_op, FDSS08_nwt_ut, FDSS09_nwt_op, FDSS09_nwt_ut, 
       FDSS10_nwt_op, FDSS10_nwt_ut, FDSS11_nwt_op, FDSS11_nwt_ut]
#%%
fig, ax = plt.subplots(2,2, figsize=(60,50))
bp1 = ax[0,0].boxplot(d11, patch_artist=True)
bp1['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][1].set(color='green', linewidth=2, facecolor = 'olive', hatch = '/')
ax[0,0].set_title('Sagehen Creek Watershed (SCW), CA',fontsize=50)
ax[0,0].set_xticklabels(['0pen2016','underTree2016'], fontsize=40)#, rotation=25[1,2]
ax[0,0].set_yticklabels(np.arange(60,150,10),fontsize=40)
ax[0,0].set_ylabel('Snow disappearance day (Julian days)', fontsize=50)
ax[0,0].legend([bp1["boxes"][0], bp1["boxes"][1]], ['Open Areas', 'Under Trees'],fontsize = 40,loc='upper left')#, loc='upper right'

bp5 = ax[0, 1].boxplot(d12, patch_artist=True)
bp5['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp5['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp5['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp5['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp5['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp5['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
ax[0, 1].set_title('King River Experimental Watershed (KREW), CA',fontsize=50)
ax[0, 1].set_xticks([1.5,3.5,5.5])
ax[0, 1].set_xticklabels(['2010','2011','2012'], fontsize=50)#, position =[1,2,3], ha = 'left''right', rotation=25[1,2]
ax[0, 1].set_yticklabels(np.arange(100,190,10),fontsize=50)
ax[0,1].legend([bp5["boxes"][0], bp5["boxes"][1]], ['Open Areas', 'Under Trees'],fontsize = 40,loc='upper left')#, loc='upper right'

bp2 = ax[1, 0].boxplot(d21, patch_artist=True)
bp2['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][6].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][7].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][8].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp2['boxes'][9].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp2['boxes'][10].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp2['boxes'][11].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
ax[1, 0].set_title('Jemez (JMZ), NM',fontsize=50)
ax[1, 0].set_xticks([1.5,3.5,5.5,7.5,9.5,11.5])
ax[1, 0].set_xticklabels(['2005','2006','2008','2009','2010','2011'],
                          fontsize=50)#, rotation=40'open2016','open2017',
ax[1, 0].set_yticklabels(np.arange(30,170,10),fontsize=50)
ax[1, 0].set_ylabel('Snow disappearance day (Julian days)', fontsize=50)
ax[1,0].legend([bp2["boxes"][0], bp2["boxes"][1]], ['Open Areas', 'Under Trees'],fontsize = 40,loc='upper left')#, loc='upper right'

bp4 = ax[1, 1].boxplot(d22, patch_artist=True)
bp4['boxes'][0].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp4['boxes'][1].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp4['boxes'][2].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][3].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][4].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp4['boxes'][5].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
bp4['boxes'][6].set(color='wheat', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][7].set(color='darkgreen', linewidth=2, facecolor = 'pink', hatch = '/')
bp4['boxes'][8].set(color='wheat', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp4['boxes'][9].set(color='darkgreen', linewidth=2, facecolor = 'olive', hatch = '/')
ax[1, 1].set_title('Niwot Ridge Creek (NRC), C0',fontsize=50)
ax[1, 1].set_xticks([1.5,3.5,5.5,7.5,9.5])
ax[1, 1].set_xticklabels(['2007','2008','2009','2010','2011'],
            fontsize=50)#, rotation=40'open2016','open2017',
ax[1, 1].set_yticklabels(np.arange(125,210,10),fontsize=50)
ax[1,1].legend([bp4["boxes"][0], bp4["boxes"][1]], ['Open Areas', 'Under Trees'],fontsize = 40,loc='upper left')#, loc='upper right'

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/obsData_dsd.png')





