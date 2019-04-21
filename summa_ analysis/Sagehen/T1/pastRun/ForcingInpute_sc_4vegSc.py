#%matplotlib inline    /bin/bash runTestCases_dockerSC.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
#%% hru names
p1 = [0.5] #[0.5,0.45,0.5,0.45] #LAIMIN
p2 = [6] #[6,5,6,5] #LAIMAX
p3 = [1] #[1,0.5,1,0.5] #winterSAI
p4 = [5] #[5,2,5,2] #summerLAI
p5 = [25] #[25,25,10,10] #heightCanopyTop
p6 = [3] #[3,3,2.5,2.5] #heightCanopyBottom 20% of top
p7 = [35] #[35,20,30,20] #maxMassVegetation         |      25.0000 |       1.0000 |      50.0000

p8 = [4.5,6.6,8] #refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)
p9 = [0.3,0.45,0.6] #throughfallScaleSnow
p10 = [700,874,950] #specificHeatVeg   j/kg k         |     874.0000 |     500.0000 |    1500.0000
p11 = [0.02,0.04,0.06] #leafDimension             |       0.0400 |       0.0100 |       0.1000

p12 = [0.5] #0.1#z0Canopy                  |       0.1000

#best param so far : #pt12_112112121 # pt12_112212121 #pt22_112112121 #pt2
p13 = [273.66] #273.66 tempCritRain	
p14 = [1.1] # 1.045 frozenPrecipMultip	

p15 = [2] #2, 3, 4] #mw_exp exponent for meltwater flow

p16 = [0.94] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500
p17 = [0.89] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p18 = [0.75] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500 ?????????????????????
p19 = [0.8] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p20 = [0.35] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p21 = [0.3] #0.5albedoSootLoad
p22 = [3]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p23 = [100000]#,200000,400000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p24 = [0.65] #albedoMinWinter           |       0.6500 |       0.6000 |       1.0000
p25 = [0.5] #albedoMinSpring           |       0.5000 |       0.3000 |       1.0000

p26 = [60] #newSnowDenMin 
p27= [75] #newSnowDenMult            |      75.0000 |      25.0000 |      75.0000

def hru_ix_ID(p1, p2, p3, p4):
    ix1 = np.arange(1,len(p1)+1)
    ix2 = np.arange(1,len(p2)+1)
    ix3 = np.arange(1,len(p3)+1)
    ix4 = np.arange(1,len(p4)+1)

    c = list(itertools.product(ix1,ix2,ix3,ix4))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p8, p9, p10, p11)
#
hru_num = np.size(hruidxID)
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
#%%
scT1fd_in = Dataset('NewData_ava\T1_2016\shT1_osh_test.nc')
scT4fd = Dataset("NewData_ava\shT4_forcing_t2017.nc")
#for varname in scT1fd_in.variables.keys():
#    var = scT1fd_in.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
#%% sagehen timming
scT1time = scT1fd_in.variables['time'][:] # get values
t_unitST1 = scT1fd_in.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = scT1fd_in.variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueScT1 = num2date(scT1time, units=t_unitST1, calendar=t_cal)
scT1date_in = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueScT1]

#%% make new nc file T4 lat: 39.42222°  lon: 120.2989°
new_fc_sc = Dataset("forcing_sagehenCreekT4_4veg.nc",'w',format='NETCDF3_CLASSIC')
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
testfd = Dataset("forcing_sagehenCreekT4_4veg.nc")

#print testfd.file_format
# read out variables, data types, and dimensions of original forcing netcdf
for varname in testfd.variables.keys():
    var = testfd.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)
#humidity = testfd.variables['spechum'][:,0]
#humidity1 = testfd1.variables['spechum'][:,0]
testfd.variables['airtemp'][:]
#plt.figure(figsize=(20,15))
#plt.plot(humidity)
#plt.plot(humidity1)








