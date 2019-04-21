#%matplotlib inline    /bin/bash runTestCases_dockerSC.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
#%% hru names
hruidxID = list(np.arange(101,102))
hru_num = np.size(hruidxID)
#%% #Sagehen creek basin forcing data (tower 1)
with open("hhs_jemez_vcp.csv") as scvd:
    reader = csv.reader(scvd)
    inputfile = [r for r in reader]
inputfile_fd_column = []
for csv_counter1 in range (len (inputfile)):
    for csv_counter2 in range (9):
        inputfile_fd_column.append(inputfile[csv_counter1][csv_counter2])
inputfile_fd=np.reshape(inputfile_fd_column,(len (inputfile),9))
inputfile_fd = inputfile_fd[1:]
fd_date = pd.DatetimeIndex(inputfile_fd[:,0])
fd_time = pd.to_numeric(fd_date, errors='raise', downcast=None)
fd_at = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,1]]), columns = ['at'])
fd_ppt = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,2]]), columns = ['ppt'])
fd_ws = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,3]]), columns = ['ws'])
fd_sh = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,4]]), columns = ['sh'])
fd_ap = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,5]]), columns = ['ap'])
fd_swr = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,6]]), columns = ['swr'])
fd_lwr = pd.DataFrame(np.array([[float(value)] for value in inputfile_fd[:,7]]), columns = ['lwr'])
allFd_df05 = pd.concat([fd_at , fd_ppt, fd_ws, fd_sh, fd_ap, fd_swr, fd_lwr], axis=1)
allFd_df05.set_index(fd_date,inplace=True)

allFd_df_hour = allFd_df05.resample('H').mean()
ppt_hour = pd.DataFrame(2*allFd_df_hour['ppt'][:])
allFd_df_hour['ppt']= ppt_hour['ppt']
#%% sagehen timming
sa_fd = Dataset('SwampAngel_forcing_sa2.nc')
TimeJp = sa_fd.variables['time'][37272:] # get values
t_unit = sa_fd.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = sa_fd.variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueJp = num2date(TimeJp, units=t_unit, calendar=t_cal)
dateJp = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueJp]
#%% sagehen creek forcing data columns=['pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum']
#Temp and ppt average to test data
#temp_sc = scT1fd.variables['airtemp'][:]
#temp_data = pd.DataFrame(temp_sc,index=pd.DatetimeIndex(scDateT1))
##temp_data=pd.Series(pd.DataFrame(temp_sc),index=pd.DatetimeIndex(scDate))
#temp_meanyr=temp_data.resample("A").mean()
#
#ppt_sc = scT1fd.variables['pptrate'][:]
#ppt_data = pd.DataFrame(ppt_sc,index=pd.DatetimeIndex(scDateT1))
#ppt_meanyr=ppt_data.resample("A").sum()
#%% make new nc file
new_fc_sc = Dataset("sagehenCreekT1_forcing.nc",'w',format='NETCDF3_CLASSIC')
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
# read out to compare with original
#for varname in new_fc_sa.variables.keys():
#    var = new_fc_sa.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
#%% define hru id, time step (1hr), lat and lon
step = np.array([3600])

lat_sa = np.array(scT1fd.variables['latitude'][:])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)

long_sa = np.array(scT1fd.variables['longitude'][:])
len_lon= np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)

#%% assign newly created variables with lists of values from NLDAS and Sagehen data
hruid[:] = hruidxID 
lat[:] = len_lat
lon[:] = len_lon
ds[:] = step

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


#%%precipitation calibration
#ppt_sa0 = np.array(sa_df['pptrate'])
#ppt_sa1 = []
#for cpi in range (np.size(ppt_sa0)):
#    if at_sa [cpi] <= 274 :
#        ppt_sa1.append(ppt_sa0[cpi]+ppt_sa0[cpi]*(1-((np.exp(4.61-0.04*((ws_sa[cpi])**1.75)))/100)))
#    else: ppt_sa1.append(ppt_sa0[cpi])
#ppt_sa = np.array(ppt_sa1)
#ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
#pptrate[:] = ppt_sa_hru

#ppt_sa0 = np.array(sa_df['pptrate'])
#ppt_sa1 = []
#for cpi in range (np.size(ppt_sa0)):
#    if at_sa [cpi] <= 274 :
#        ppt_sa1.append(ppt_sa0[cpi]+ppt_sa0[cpi]*(1-((np.exp(4.61-0.16*(ws_sa[cpi]**1.28)))/100)))
#    else: ppt_sa1.append(ppt_sa0[cpi])
#ppt_sa = np.array(ppt_sa1)
#ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
#pptrate[:] = ppt_sa_hru

#%%specific humididty calculations***********************************************
#at0_sb = sbFD.variables['airtemp'][:]
#
#e_t = (ap_sb * sh_sb)/0.622
#p_da = ap_sb - e_t
#e_star_t = 611*(np.exp((17.27*(at0_sb-273.15))/(at0_sb-273.15+237.3)))
#rh = e_t/e_star_t
#
#e_star_t2 = 611*(np.exp((17.27*(at0_sb+2-273.15))/(at0_sb+2-273.15+237.3)))
#e_star_t4 = 611*(np.exp((17.27*(at0_sb+4-273.15))/(at0_sb+4-273.15+237.3)))
#
#e_t2 = rh * e_star_t2
#e_t4 = rh * e_star_t4
#
#p_t2 = p_da + e_t2
#p_t4 = p_da + e_t4
#
#sh_t2 = 0.622 * e_t2 / p_t2
#sh_t4 = 0.622 * e_t4 / p_t4
#%%******************************************************************************
test = new_fc_sc.variables['pptrate'][:]

# close the file to write it
new_fc_sc.close()
#%%
testfd = Dataset("sagehenCreekT1_forcing.nc")
testfd1 = Dataset("sagehenCreek_forcing.nc")
testfd2 = Dataset("shT1_force1_2.nc")
testfd3 = Dataset("shT1_osh_test.nc")

#print testfd.file_format
# read out variables, data types, and dimensions of original forcing netcdf
for varname in testfd.variables.keys():
    var = testfd.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)
#humidity = testfd.variables['spechum'][:,0]
#humidity1 = testfd1.variables['spechum'][:,0]
humidity2 = testfd2.variables['spechum'][:,0]
humidity3 = testfd3.variables['spechum'][:,0]

fsh = testfd.variables['spechum'][:,0]
#plt.figure(figsize=(20,15))
#plt.plot(humidity)
#plt.plot(humidity1)








