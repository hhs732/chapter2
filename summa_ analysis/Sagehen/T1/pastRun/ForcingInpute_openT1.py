#%matplotlib inline    /bin/bash runTestCases_dockerSC.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
#%% hru names
hruidxID = [101]
hru_num = np.size(hruidxID)

#%% ppt observed data T4
with open("ppt_T4.csv") as scvdp:
    readerp = csv.reader(scvdp)
    raw_ppt = [rp for rp in readerp]
sc_ppt_column = []
for csv_counterp1 in range (len (raw_ppt)):
    for csv_counterp2 in range (2):
        sc_ppt_column.append(raw_ppt[csv_counterp1][csv_counterp2])
sc_ppt=np.reshape(sc_ppt_column,(len (raw_ppt),2))
sc_ppt = sc_ppt[1:]
sc_ppt_obs_date = pd.DatetimeIndex(sc_ppt[:,0])
sc_ppt_obs = [float(value) for value in sc_ppt[:,1]]
ppt_df = pd.DataFrame(sc_ppt_obs, columns = ['ppt (mm)']) 
ppt_df.set_index(sc_ppt_obs_date,inplace=True)

ppt = []
for indxswe in range (len(ppt_df)-1):
    if (ppt_df['ppt (mm)'][indxswe+1]-ppt_df['ppt (mm)'][indxswe])<0:
        ppt.append(0)
    else: ppt.append(ppt_df['ppt (mm)'][indxswe+1]-ppt_df['ppt (mm)'][indxswe])
ppt_df = pd.DataFrame(ppt, columns = ['ppt (mm)']) 
ppt_df.set_index(sc_ppt_obs_date[1:],inplace=True)  

ppt_hr_df = ppt_df.resample('h').sum() 
#%% #Sagehen creek basin forcing data (tower 1) from 2014-1-7
with open("hhs_scT4_fd.csv") as scvd:
    reader = csv.reader(scvd)
    input_scT1 = [r for r in reader]
scT1_fd_column = []
for csv_counter1 in range (len (input_scT1)):
    for csv_counter2 in range (9):
        scT1_fd_column.append(input_scT1[csv_counter1][csv_counter2])
scT1_fd=np.reshape(scT1_fd_column,(len (input_scT1),9))
scT1_fd = scT1_fd[1:]
#scT1_fd_date = pd.DatetimeIndex(scT1_fd[:,1])
scT1_fd_time = scT1_fd[:,1]
scT1_lwr = np.array([[float(value)] for value in scT1_fd[:,2]])
scT1_swr = np.array([[float(value)] for value in scT1_fd[:,3]])
scT1_at = np.array([[float(value)] for value in scT1_fd[:,4]])
scT1_ppt = np.array([[float(value)] for value in scT1_fd[:,5]])
scT1_ws = np.array([[float(value)] for value in scT1_fd[:,6]])
scT1_ap = np.array([[float(value)] for value in scT1_fd[:,7]])
scT1_sh = np.array([[float(value)] for value in scT1_fd[:,8]])
#swe_obs_df = pd.DataFrame(sc_swe_obs, columns = ['observed swe']) 
#swe_obs_df.set_index(sc_swe_obs_date,inplace=True)
#%%T4 lat: 39.42222°  lon: 120.2989°; T1 [ 39.4321] [-120.2411]
scT1fd_in = Dataset('NewData_ava\T1_2016\shT1_osh_test.nc')

#%% sagehen timming
scT1time = scT1fd_in.variables['time'][:] # get values
t_unitST1 = scT1fd_in.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = scT1fd_in.variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueScT1 = num2date(scT1time, units=t_unitST1, calendar=t_cal)
scT1date_in = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueScT1]

#%% make new nc file 
new_fc_sc = Dataset("forcing_sagehenCreekT4_open.nc",'w',format='NETCDF3_CLASSIC')
# define dimensions 
hru = new_fc_sc.createDimension('hru', hru_num)
time = new_fc_sc.createDimension('time', None)
# define variables
hruid = new_fc_sc.createVariable('hruId', np.int32,('hru',))
lat = new_fc_sc.createVariable('latitude', np.float64,('hru',))
lon = new_fc_sc.createVariable('longitude', np.float64,('hru',))
ds = new_fc_sc.createVariable('data_step', np.float64)
times = new_fc_sc.createVariable('time', np.float64,('time',))
lwrad = new_fc_sc.createVariable('LWRadAtm', np.float64,('time','hru'), fill_value = -999.0)
swrad = new_fc_sc.createVariable('SWRadAtm', np.float64,('time','hru'), fill_value = -999.0)
airpres = new_fc_sc.createVariable('airpres', np.float64,('time','hru'), fill_value = -999.0)
airtemp = new_fc_sc.createVariable('airtemp', np.float64,('time','hru'), fill_value = -999.0)
pptrate = new_fc_sc.createVariable('pptrate', np.float64,('time','hru'), fill_value = -999.0)
spechum = new_fc_sc.createVariable('spechum', np.float64,('time','hru'), fill_value = -999.0)
windspd = new_fc_sc.createVariable('windspd', np.float64,('time','hru'), fill_value = -999.0)
# give variables units
times.units = 'days since 1990-01-01 00:00:00'
ds.units = 'seconds'
lwrad.units = 'W m-2'
swrad.units = 'W m-2'
airpres.units = 'Pa'
airtemp.units = 'K'
pptrate.units = 'kg m-2 s-1'
spechum.units = 'g g-1'
windspd.units = 'm s-1'
# give variables value type
lwrad.vtype = 'scalarv'
swrad.vtype = 'scalarv'
airpres.vtype = 'scalarv'
airtemp.vtype = 'scalarv'
pptrate.vtype = 'scalarv'
spechum.vtype = 'scalarv'
windspd.vtype = 'scalarv'

#%% T4 lat: 39.42222°  lon: 120.2989° elev:2370; T1 [ 39.4321] [-120.2411] elev = 1936m 
step = np.array([3600])

lat_sa = np.array([39.4222])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)

long_sa = np.array([-120.2989])
len_lon= np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)

#%% assign newly created variables with lists of values from NLDAS and Sagehen data
hruid[:] = hruidxID 
lat[:] = len_lat
lon[:] = len_lon
ds[:] = step

TimeScT1 = scT1fd_in.variables['time'][43152:] # get values
new_ix = np.array(TimeScT1)
times[:] = new_ix

lwr_sa = np.array(scT1_lwr)
lwr_sa_hru = np.repeat(lwr_sa[:,np.newaxis], hru_num, axis=1)
lwrad[:] = lwr_sa_hru

swr_sa = np.array(scT1_swr)
swr_sa_hru = np.repeat(swr_sa[:,np.newaxis], hru_num, axis=1)
swrad[:] = swr_sa_hru

ap_sa = np.array(scT1_ap)
ap_sa_hru = np.repeat(ap_sa[:,np.newaxis], hru_num, axis=1)
airpres[:] = ap_sa_hru

at_sa = np.array(scT1_at)
at_sa_hru = np.repeat(at_sa[:,np.newaxis], hru_num, axis=1) 
airtemp[:] = at_sa_hru

ws_sa = np.array(scT1_ws)
ws_sa_hru = np.repeat(ws_sa[:,np.newaxis], hru_num, axis=1) 
windspd[:] = ws_sa_hru

ppt_sa = np.array(scT1_ppt)
ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
pptrate[:] = ppt_sa_hru

#testfd1 = Dataset("sagehenCreek_forcing.nc")
#humidity1 = testfd1.variables['spechum'][:,0]

#sh_sa[np.isnan(sh_sa)] = 0
#for ix in range (len(sh_sa)):
#    if sh_sa[ix]==0:
#        sh_sa[ix] = humidity1[ix]
sh_sa = np.array(scT1_sh)
sh_sa_hru = np.repeat(sh_sa[:,np.newaxis], hru_num, axis=1) 
spechum[:] = sh_sa_hru

#%%******************************************************************************
test = new_fc_sc.variables['pptrate'][:]

# close the file to write it
new_fc_sc.close()
#%%
testfd = Dataset("forcing_sagehenCreekT4_open.nc")

print testfd.variables['latitude'][:]






