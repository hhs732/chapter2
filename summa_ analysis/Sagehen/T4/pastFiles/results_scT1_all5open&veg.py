# -*- coding: utf-8 -*-
"""
Created on Sun Sep 09 12:15:00 2018

@author: HHS
"""

#%% SWE observed data T4
with open("input_SWE_T4.csv") as scvd:
    reader2 = csv.reader(scvd)
    raw_swe2 = [r2 for r2 in reader2]
sc_swe_column2 = []
for csv_counter12 in range (len (raw_swe2)):
    for csv_counter22 in range (2):
        sc_swe_column2.append(raw_swe2[csv_counter12][csv_counter22])
sc_swe2=np.reshape(sc_swe_column2,(len (raw_swe2),2))
sc_swe2 = sc_swe2[1:]
sc_swe_obs_date2 = pd.DatetimeIndex(sc_swe2[:,0])
sc_swe_obs2 = [float(value2) for value2 in sc_swe2[:,1]]
swe_obs_df2 = pd.DataFrame(sc_swe_obs2, columns = ['observed swe T4']) 
swe_obs_df2.set_index(sc_swe_obs_date2,inplace=True)