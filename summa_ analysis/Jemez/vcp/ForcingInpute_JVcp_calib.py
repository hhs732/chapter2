#%matplotlib inline    /bin/bash runTestCases_dockerSC.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
#%% hru names
p1 = [273.16,273.5] #,273.5,273.66 tempCritRain	
p2 = [1,1.05] # 1.045 frozenPrecipMultip	

p3 = [2,3,4] #2, 3, 4] #mw_exp exponent for meltwater flow
p4 = [0.28,0.4] #windReductionParam        |       0.2800 |       0.0000 |       1.0000
p5 = [0.001,0.002] #zsnow

p6 = [0.89,0.95] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p7 = [1,3,5]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p8 = [200000,500000,1000000]#,350000,200000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p9 = [0.3,0.5] #0.5albedoSootLoad

p10 = [0.75] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500
p11 = [0.65] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p12 = [0.3] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p13 = [0.85] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500

p14 = [0.01] #winterSAI
p15 = [0.1] #summerLAI
p16 = [0.01] #LAIMIN
p17 = [1] #LAIMAX
p18 = [0.3] #heightCanopyTop
p19 = [0.03] #heightCanopyBottom to calculate wind speed, wind speed reduction; coupute roughness length of the veg canopy; neutral ground resistance; canopy air temp;
p20 = [0.02] #0.1#z0Canopy  

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

hruidxID = hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9)#, p10, p11, p12, p13, p14, p15, p16, p17 p18, p19, p20, p21
#
hru_num = np.size(hruidxID)

#%% #Sagehen creek basin forcing data (tower 1) from 2014-1-7
with open("hhs_jemez_vcp.csv") as scvd:
    reader = csv.reader(scvd)
    input_scT1 = [r for r in reader]
scT1_fd_column = []
for csv_counter1 in range (len (input_scT1)):
    for csv_counter2 in range (8):
        scT1_fd_column.append(input_scT1[csv_counter1][csv_counter2])
scT1_fd=np.reshape(scT1_fd_column,(len (input_scT1),8))
scT1_fd = scT1_fd[1:]
scT1_fd_date = pd.DatetimeIndex(scT1_fd[:,0])
#scT1_fd_time = [float(x) for x in scT1_fd_date.values.astype(float)]
scT1_lwr = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,1]]), index = scT1_fd_date)
scT1_lwr_hr = scT1_lwr.resample('H').mean()
scT1_swr = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,2]]), index = scT1_fd_date)
scT1_swr_hr = scT1_swr.resample('H').mean()
scT1_at = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,3]]), index = scT1_fd_date)
scT1_at_hr = scT1_at.resample('H').mean()
scT1_ppt = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,4]]), index = scT1_fd_date)
scT1_ppt_hr = scT1_ppt.resample('H').sum()
scT1_ws = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,5]]), index = scT1_fd_date)
scT1_ws_hr = scT1_ws.resample('H').mean()
scT1_ap = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,6]]), index = scT1_fd_date)
scT1_ap_hr = scT1_ap.resample('H').mean()
scT1_sh = pd.DataFrame(np.array([[float(value)] for value in scT1_fd[:,7]]), index = scT1_fd_date)
scT1_sh_hr = scT1_sh.resample('H').mean()
#%%T4 lat: 39.42222°  lon: 120.2989°; T1 [ 39.4321] [-120.2411]
scT1fd_in = Dataset('shT1_osh_test.nc')
sbfd_in = Dataset('SenatorBeck_forcing.nc')

time_sb = sbfd_in.variables['time'][0:61079]
time_sc = scT1fd_in.variables['time'][:]
total_time = np.concatenate((time_sb, time_sc), axis=None)
#%% Jemez timming
scT1time = total_time[51647:69167] # get values
t_unitST1 = scT1fd_in.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = scT1fd_in.variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueScT1 = num2date(scT1time, units=t_unitST1, calendar=t_cal)
scT1date = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueScT1]

#%% make new nc file 
new_fc_sc = Dataset("forcing_JmzVcp_open_calib.nc",'w',format='NETCDF3_CLASSIC')
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

#%% http://www.fluxdata.org:8080/sitepages/siteinfo.aspx?us-vcm
step = np.array([3600])

lat_sa = np.array([35.8642])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)

long_sa = np.array([-106.5967])
len_lon= np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)

#%% assign newly created variables with lists of values from NLDAS and Sagehen data
hruid[:] = hruidxID 
lat[:] = len_lat
lon[:] = len_lon
ds[:] = step

TimeScT1 = scT1time # get values
new_ix = np.array(scT1time)
times[:] = new_ix

lwr_sa = np.array(scT1_lwr_hr)
lwr_sa_hru = np.repeat(lwr_sa[:,np.newaxis], hru_num, axis=1)
lwrad[:] = lwr_sa_hru

swr_sa = np.array(scT1_swr_hr)
swr_sa_hru = np.repeat(swr_sa[:,np.newaxis], hru_num, axis=1)
swrad[:] = swr_sa_hru

ap_sa = np.array(scT1_ap_hr)
ap_sa_hru = np.repeat(ap_sa[:,np.newaxis], hru_num, axis=1)
airpres[:] = ap_sa_hru

at_sa = np.array(scT1_at_hr)
at_sa_hru = np.repeat(at_sa[:,np.newaxis], hru_num, axis=1) 
airtemp[:] = at_sa_hru

ws_sa = np.array(scT1_ws_hr)
ws_sa_hru = np.repeat(ws_sa[:,np.newaxis], hru_num, axis=1) 
windspd[:] = ws_sa_hru

ppt_sa = np.array(scT1_ppt_hr)
ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
pptrate[:] = ppt_sa_hru

sh_sa = np.array(scT1_sh_hr)
sh_sa_hru = np.repeat(sh_sa[:,np.newaxis], hru_num, axis=1) 
spechum[:] = sh_sa_hru

#testfd1 = Dataset("sagehenCreek_forcing.nc")
#humidity1 = testfd1.variables['spechum'][:,0]

#sh_sa[np.isnan(sh_sa)] = 0
#for ix in range (len(sh_sa)):
#    if sh_sa[ix]==0:
#        sh_sa[ix] = humidity1[ix]
#%%******************************************************************************
test = new_fc_sc.variables['pptrate'][:]

# close the file to write it
new_fc_sc.close()
#%%
testfd = Dataset("forcing_JmzVcp_open_calib.nc")

print testfd.variables['LWRadAtm'][0:5]
print testfd.variables['SWRadAtm'][0:5]
print testfd.variables['airtemp'][0:5]
print testfd.variables['windspd'][0:5]


#%%time and date
#with open("hhs_time.csv") as scvd:
#    reader_time = csv.reader(scvd)
#    input_time = [t for t in reader_time]
#time_column = []
#for csv1 in range (len (input_time)):
#    for csv2 in range (2):
#        time_column.append(input_time[csv1][csv2])
#totalTime=np.reshape(time_column,(len (input_time),2))
#totalTimeDate = totalTime[1:]
#total_date = pd.DatetimeIndex(totalTimeDate[:,0])
#total_time = totalTimeDate[:,1]




