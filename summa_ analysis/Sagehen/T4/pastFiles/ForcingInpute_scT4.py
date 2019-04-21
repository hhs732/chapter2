#%matplotlib inline    /bin/bash runTestCases_dockerSC.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
#%% hru names
p1 = [273,16,273.66,273.75] #273.66 tempCritRain	
p2 = [1,1.05,1.1] # 1.045 frozenPrecipMultip	

p3 = [2,3] #2, 3, 4] #mw_exp exponent for meltwater flow

p4 = [0.89,0.94] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500
p5 = [0.85,0.9,0.94] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p6 = [0.7,0.75] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500
p7 = [0.6,0.7] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p8 = [0.3,0.4] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p9 = [0.3,0.5] #0.5albedoSootLoad
p10 = [3]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p11 = [600000]#,350000,200000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 

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
#p5 = [25] #maxMassVegetation 
#p8 = [2] #rootingDepth
#p21 = [6.6] #refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)
#p22 = [0.89] #throughfallScaleSnow
#p23 = [0.4] #ratioDrip2Unloading       |       0.4000 |       0.0000 |       1.0000
#p24 = [1] #rootDistExp |       1.0000 |       0.0100 |       1.0000
#p25 = [874] #specificHeatVeg   j/kg k         |     874.0000 |     500.0000 |    1500.0000
#p26 = [0.04] #leafDimension             |       0.0400 |       0.0100 |       0.1000

#p19 = [0.35] #0.2, 0.4 , 0.6] #fixedThermalCond_snow used in tcs_smnv model, but we are using Jordan model  
#p28 = [0.28] #windReductionParam        |       0.2800 |       0.0000 |       1.0000
#p29 = [0.2] #critRichNumber
#p30 = [1] #Mahrt87_eScale  
#p31 = [0.06] #Fcapil
#p32 = [0.015] #k_snow
#p33 = [9.4] #Louis79_bparam  
#p34 = [5.3] #Louis79_cStar 

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
    #ix10 = np.arange(1,len(p10)+1)
    #ix11 = np.arange(1,len(p11)+1)

    c = list(itertools.product(ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9))#,ix10,ix11,ix12,ix13,ix14,ix15,ix16,ix17,ix18,ix19,ix20,ix21))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p1, p2, p3, p4, p5, p6, p7, p8, p9)#, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21)
#
hru_num = np.size(hruidxID)
#%% #Sagehen creek basin forcing data (tower 1) from 2014-1-7
with open("hhs_scT4_fd.csv") as scvd:
    reader = csv.reader(scvd)
    input_scT1 = [r for r in reader]
scT1_fd_column = []
for csv_counter1 in range (len (input_scT1)):
    for csv_counter2 in range (8):
        scT1_fd_column.append(input_scT1[csv_counter1][csv_counter2])
scT4_fd=np.reshape(scT1_fd_column,(len (input_scT1),8))
scT4_fd = scT4_fd[1:]
scT4_time = np.array([float(value) for value in scT4_fd[:,0]])
scT4_lwr = np.array([[float(value)] for value in scT4_fd[:,1]])
scT4_swr = np.array([[float(value)] for value in scT4_fd[:,2]])
scT4_at = np.array([[float(value)] for value in scT4_fd[:,3]])
scT4_ppt = np.array([[float(value)] for value in scT4_fd[:,4]])
scT4_ws = np.array([[float(value)] for value in scT4_fd[:,5]])
scT4_ap = np.array([[float(value)] for value in scT4_fd[:,6]])
scT4_sh = np.array([[float(value)] for value in scT4_fd[:,7]])
#swe_obs_df = pd.DataFrame(sc_swe_obs, columns = ['observed swe']) 
#swe_obs_df.set_index(sc_swe_obs_date,inplace=True)
#%%
scT4fd_in16 = Dataset("NewData_ava/T4_2016/t4n_forcing_testnp.nc")
scT4fd_in17 = Dataset("NewData_ava/T4_2017/t4n_forcing_t2017.nc")
#%% sagehen timming
scT4time = scT4_time.copy()
#scT4time = scT4fd_in17.variables['time'][:]

t_unitST1 = scT4fd_in17.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = scT4fd_in17.variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueScT4 = num2date(scT4time, units=t_unitST1, calendar=t_cal)
scT4date_in = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueScT4]
#%% I calculated specific humidy. difference between data is less than 1%.
#time = pd.DataFrame(scT4fd_in16.variables['time'][:], columns = ['time'])
#scT1_lwr_in = pd.DataFrame(scT4fd_in16.variables['LWRadAtm'][:], columns = ['lwr'])
#scT1_swr_in = pd.DataFrame(scT4fd_in16.variables['SWRadAtm'][:], columns = ['swr'])
#scT1_at_in = pd.DataFrame(scT4fd_in16.variables['airtemp'][:], columns = ['at'])
#scT1_ap_in = pd.DataFrame(scT4fd_in16.variables['airpres'][:], columns = ['ap'])
#scT1_ppt_in = pd.DataFrame(scT4fd_in16.variables['pptrate'][:], columns = ['ppt'])
#scT1_sh_in = pd.DataFrame(scT4fd_in16.variables['spechum'][:], columns = ['sh'])
#scT1_ws_in = pd.DataFrame(scT4fd_in16.variables['windspd'][:], columns = ['ws'])
#
#scT4_allfd_in = pd.concat([time, scT1_lwr_in , scT1_swr_in, scT1_at_in, scT1_ap_in, scT1_ppt_in, scT1_sh_in, scT1_ws_in], axis=1)
##scT4_allfd_in.set_index(pd.DatetimeIndex(time),inplace=True)

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
new_fc_sc = Dataset("sagehenCreekT4_forcing.nc",'w',format='NETCDF3_CLASSIC')
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

lat_sa = np.array(scT4fd_in17.variables['latitude'][:])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)

long_sa = np.array(scT4fd_in17.variables['longitude'][:])
len_lon= np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)

#%% assign newly created variables with lists of values from NLDAS and Sagehen data
hruid[:] = hruidxID 
lat[:] = len_lat
lon[:] = len_lon
ds[:] = step

TimeSc = scT4time.copy()
new_ix = np.array(TimeSc)
times[:] = new_ix

lwr_sa = np.array(scT4_lwr)
lwr_sa_hru = np.repeat(lwr_sa[:,np.newaxis], hru_num, axis=1)
lwrad[:] = lwr_sa_hru

swr_sa = np.array(scT4_swr)
swr_sa_hru = np.repeat(swr_sa[:,np.newaxis], hru_num, axis=1)
swrad[:] = swr_sa_hru

ap_sa = np.array(scT4_ap)
ap_sa_hru = np.repeat(ap_sa[:,np.newaxis], hru_num, axis=1)
airpres[:] = ap_sa_hru

at_sa = np.array(scT4_at)
at_sa_hru = np.repeat(at_sa[:,np.newaxis], hru_num, axis=1) 
airtemp[:] = at_sa_hru

ws_sa = np.array(scT4_ws)
ws_sa_hru = np.repeat(ws_sa[:,np.newaxis], hru_num, axis=1) 
windspd[:] = ws_sa_hru

ppt_sa = np.array(scT4_ppt)
ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
pptrate[:] = ppt_sa_hru

#testfd1 = Dataset("sagehenCreek_forcing.nc")
#humidity1 = testfd1.variables['spechum'][:,0]

#sh_sa[np.isnan(sh_sa)] = 0
#for ix in range (len(sh_sa)):
#    if sh_sa[ix]==0:
#        sh_sa[ix] = humidity1[ix]
sh_sa = np.array(scT4_sh)
sh_sa_hru = np.repeat(sh_sa[:,np.newaxis], hru_num, axis=1) 
spechum[:] = sh_sa_hru

#%%******************************************************************************
test = new_fc_sc.variables['LWRadAtm'][:]

# close the file to write it
new_fc_sc.close()
#%%
testfd = Dataset("sagehenCreekT4_forcing.nc")
testfd1 = Dataset("sagehenCreek_forcing.nc")

#print testfd.file_format
# read out variables, data types, and dimensions of original forcing netcdf
for varname in testfd.variables.keys():
    var = testfd.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)
#humidity = testfd.variables['spechum'][:,0]
#humidity1 = testfd1.variables['spechum'][:,0]
print testfd.variables['longitude'][:]
#plt.figure(figsize=(20,15))
#plt.plot(humidity)
#plt.plot(humidity1)








