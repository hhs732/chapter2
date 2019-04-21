import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
#%% read
with open("tree.csv") as scvdT:
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
with open("open.csv") as scvd0:
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
d1 = [FDSS_ls0,LDSS_ls0,FDSS_lsT,LDSS_lsT]

fig, ax = plt.subplots(1,1, figsize=(20,15))
bp1 = ax.boxplot(d1, patch_artist=True)
bp1['boxes'][0].set(color='orange', linewidth=2, facecolor = 'skyblue', hatch = '/')
bp1['boxes'][1].set(color='darkred', linewidth=2, facecolor = 'olive', hatch = '/')
bp1['boxes'][2].set(color='turquoise', linewidth=2, facecolor = 'pink', hatch = '/')
bp1['boxes'][3].set(color='green', linewidth=2, facecolor = 'pink', hatch = '/')

plt.title('Sagehen DTS snow presence data',fontsize=40)

plt.xticks([1,2,3,4], ['0pen_FDSD','0pen_LDSOG','underTree_FDSD','underTree_LDSOG'], fontsize=30, rotation=25)#'open2016','open2017',
plt.yticks(fontsize=30)
#plt.xlabel('vegSc_year', fontsize=30)
plt.ylabel('julian days', fontsize=30)
plt.savefig('dtsDta_snowPresence.png')























