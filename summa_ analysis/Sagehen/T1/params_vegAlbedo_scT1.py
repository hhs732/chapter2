###       /bin/bash runTestCases_dockerScT1_calib.sh   snwDensity snwDensity
# 2007 - 2008 as wet year for sensirivity analysis 1st step
# I thinking having five scenarios, open, dense and tall, dense and short, sparse and tall, sparse and short. Assuming they are sensitive.
import numpy as np
from netCDF4 import Dataset
import itertools
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

p8 = [4] # 2,3.5 mw_exp exponent for meltwater flow

p9 = [0.80] # 0.7,0.78 albedoMax |       0.8500 |       0.7000 |       0.9500
p10 = [0.4,0.6] #albedoMinWinter     0.6000 |       0.6000 |       1.0000
p11 = [0.3,0.45] #albedoMinSpring    0.4500 |       0.3000 |       1.0000
p12 = [700000] # 250000,500000 albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p13 = [0.80] # 0.7,0.78 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p14 = [0.75] # 0.55,0.7 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500
p15 = [0.6] # 0.5,0.57 albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p16 = [0.4] # 0.2,0.38 albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p17 = [3] # 3.5,5 albedoRefresh |       1.0000 |       1.0000 |      10.0000
p18 = [0.5] # 0.3,0.5 albedoSootLoad |       0.3000 |       0.1000 |       0.5000

p19 = [9] #[6,13]refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)

p20 = [0.4] # [0.3,0.7]0.5 throughfallScaleSnow [0.3,0.45,0.6]
p21 = [0.5] #[0.4,0.8]throughfallScaleRain      |       0.6000 |       0.1000 |       0.9000

p22 = [0.35] #[0.1,0.6]ratioDrip2Unloading       |       0.4000 |       0.0000 |       1.0000
p23 = [0] #[0,0.0000014]snowUnloadingCoeff     [0,0.000001,0.0000014]   |       0.0000 |       0.0000 |       1.5d-6

p24 = [0.04] #leafDimension            0.0400 |       0.0100 |       0.1000
p25 = [0.01] #leafExchangeCoeff         |       0.0100 |       0.0010 |       0.1000
p26 = [0.3] #windReductionParam        |       0.2800 |       0.0000 |       1.0000
p27 = [3] #rootingDepth

hruidxID = np.arange(100001,100005)

hru_num = np.size(hruidxID)

#%% 
def param_fill(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27):#, p28
    b = list(itertools.product(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27))#, p28
    p1l =[]; p2l =[]; p3l =[]; p4l=[]; p5l =[]; p6l =[]; p7l =[]; p8l =[]; p9l =[]; p10l=[]; p11l =[]; p12l =[]; p13l =[]; p14l=[]; p15l = []; p16l =[]; 
    p17l=[]; p18l = []; p19l = []; p20l = []; p21l = []; p22l = []; p23l = []; p24l = []; p25l = []; p26l = []; p27l = []#; p28l = []
    for tup in b:
        p1l.append(tup[0]); p2l.append(tup[1]); p3l.append(tup[2]); p4l.append(tup[3]); p5l.append(tup[4]); p6l.append(tup[5]); p7l.append(tup[6]); 
        p8l.append(tup[7]); p9l.append(tup[8]); p10l.append(tup[9]); p11l.append(tup[10]); p12l.append(tup[11]); p13l.append(tup[12]); p14l.append(tup[13]); 
        p15l.append(tup[14]); p16l.append(tup[15]); p17l.append(tup[16]); p18l.append(tup[17]); p19l.append(tup[18]); p20l.append(tup[19]); p21l.append(tup[20]);
        p22l.append(tup[21]); p23l.append(tup[22]); p24l.append(tup[23]); p25l.append(tup[24]); p26l.append(tup[25]); p27l.append(tup[26])#; p28l.append(tup[27])
    return(p1l, p2l, p3l, p4l, p5l, p6l, p7l, p8l, p9l, p10l, p11l, p12l, p13l, p14l, p15l, p16l, p17l, p18l, p19l, p20l, p21l, p22l, p23l, p24l, p25l, p26l, p27l)#, p28l

valst1 = param_fill(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27) #, p28

#%% parameterTrial, Local attributes and initial conditions for senatore beck
pt = Dataset('summa_zParamTrial_variableDecayRate.nc')
la = Dataset('summa_zLocalAttributes_senatorSheltered.nc') #('settings/wrrPaperTestCases/figure07/summa_zLocalAttributes_riparianAspen.nc')
ic = Dataset('summa_zInitialCond.nc') #('settings/wrrPaperTestCases/figure07/summa_zInitialCond.nc')

for varname in pt.variables.keys():
    var = pt.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
paramfile = Dataset("summa_zParamTrial_scT1_vegReal.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['LAIMIN','LAIMAX','winterSAI','summerLAI','heightCanopyTop','heightCanopyBottom','maxMassVegetation',
                  'mw_exp','albedoMax','albedoMinWinter','albedoMinSpring','albedoDecayRate',
                  'albedoMaxVisible','albedoMinVisible','albedoMaxNearIR','albedoMinNearIR',
                  'albedoRefresh','albedoSootLoad','refInterceptCapSnow',
                  'throughfallScaleSnow','throughfallScaleRain','ratioDrip2Unloading','snowUnloadingCoeff',
                  'leafDimension','leafExchangeCoeff','windReductionParam','rootingDepth'] 
#'tempCritRain','z0Snow','specificHeatVeg','z0Canopy','refInterceptCapRain','newSnowDenMin'
for param in param_nam_list:
    paramfile.createVariable(param, np.float64,('hru',))

constant_params = ['rootDistExp','theta_sat','theta_res','vGn_alpha','vGn_n','k_soil','critSoilWilting','critSoilTranspire','frozenPrecipMultip']
for params in constant_params:
    paramfile.createVariable(params, np.float64,('hru',))

#%% # add values for the constant variables in HRUs for parameter Trail file
for varname in pt.variables.keys():
    var = pt.variables[varname][0]
    c = np.full((hru_num,),var)
    try :
        paramfile.variables[varname][:]=c
    except IndexError: # size of data array does not conform to slice
        pass
#%% creating changing variables and adding values for changing variables
j = 0 
for var in param_nam_list:
    paramfile.variables[var][:]=valst1[j]
    j=j+1

paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    
paramfile.close()

#%% 
varcheck = Dataset ('summa_zParamTrial_scT1_vegReal.nc')
print varcheck.variables['albedoMinVisible'][:]
print varcheck.variables['frozenPrecipMultip'][:]
#%% # local attributes file
local_atrbt = Dataset("summa_zLocalAttributes_scT1_veg.nc",'w',format='NETCDF3_CLASSIC')
# define dimensions 
hru = local_atrbt.createDimension('hru', hru_num) 
time = local_atrbt.createDimension('gru', 1)
# define variables
h2gid = local_atrbt.createVariable('hru2gruId', np.int32,('hru',))
dhruindx = local_atrbt.createVariable('downHRUindex', np.int32,('hru',))
slopeindx = local_atrbt.createVariable('slopeTypeIndex', np.int32,('hru',))
soilindx = local_atrbt.createVariable('soilTypeIndex', np.int32,('hru',))
vegindx = local_atrbt.createVariable('vegTypeIndex', np.int32,('hru',))
mh = local_atrbt.createVariable('mHeight', np.float64,('hru',))
cl = local_atrbt.createVariable('contourLength', np.float64,('hru',))
tanslope = local_atrbt.createVariable('tan_slope', np.float64,('hru',))
elev = local_atrbt.createVariable('elevation', np.float64,('hru',))
lon = local_atrbt.createVariable('longitude', np.float64,('hru',))
lat = local_atrbt.createVariable('latitude', np.float64,('hru',))
hruarea = local_atrbt.createVariable('HRUarea', np.float64,('hru',))
hruid = local_atrbt.createVariable('hruId', np.int32,('hru',))
gruid = local_atrbt.createVariable('gruId', np.int32,('gru',))
# give variables units
mh.units = 'm'
cl.units = 'm'
tanslope.units = 'm m-1'
elev.units = 'm'
lat.units = 'decimal degree north'
lon.units = 'decimal degree east'
hruarea.units = 'm^2'
#%% # add values for the constant variables in HRUs for local atribute file
for varname in la.variables.keys():
    var = la.variables[varname][0]
    #print var
    c2 = np.full((hru_num,),var)
    #print c2
    try :
        local_atrbt.variables[varname][:]=c2
    except IndexError: # size of data array does not conform to slice
        pass
    #local_atrbt.variables[varname][:]=c2
#%% T1 [ 39.4321] [-120.2411] elev : 1936m
lat_sa = np.array([39.4321])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)
lat[:] = len_lat

long_sa = np.array([-120.2411])
len_lon = np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
lon[:] = len_lon

elev_sa = np.array([1936]) 
elev_lon = np.repeat(elev_sa[:,np.newaxis], hru_num, axis=1); elev_lon=elev_lon.reshape(hru_num,)
elev[:] = elev_lon

vegindx_sc = np.array([14])
vegindx_lon = np.repeat(vegindx_sc[:,np.newaxis], hru_num, axis=1); vegindx_lon=vegindx_lon.reshape(hru_num,)
vegindx[:] = vegindx_lon

soilindx_sc = np.array([7])
soilindx_lon = np.repeat(soilindx_sc[:,np.newaxis], hru_num, axis=1); soilindx_lon=soilindx_lon.reshape(hru_num,)
soilindx[:] = soilindx_lon

mHeight_sc = np.array([7.62]) #mHeight
mHeight_len = np.repeat(mHeight_sc[:,np.newaxis], hru_num, axis=1); mHeight_len=mHeight_len.reshape(hru_num,)
mh[:] = mHeight_len

#%% # get the hru, gru, and hru2gru in local_atribute file
newgru = np.array([11111111])
local_atrbt.variables['gruId'][:] = newgru

c3 = np.repeat(newgru[:,np.newaxis], hru_num, axis=1); newlad = c3.reshape(hru_num,)
local_atrbt.variables['hru2gruId'][:] = c3

local_atrbt.variables['hruId'][:]=hruidxID

print local_atrbt.variables['hruId'][:]
local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_scT1_veg.nc')
print lacheck.variables['soilTypeIndex'][:]
print lacheck.variables['longitude'][:]
print lacheck.variables['mHeight'][:]
print lacheck.variables['hruId'][:]

#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
#lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_scT1_veg.nc",'w',format='NETCDF3_CLASSIC')
#print ic.variables.keys()

# define dimensions 
midtoto = in_condi.createDimension('midToto',8)
midsoil = in_condi.createDimension('midSoil',8)
idctoto = in_condi.createDimension('ifcToto',9)
scalarv = in_condi.createDimension('scalarv', 1)
# this is the number you will change to the number of HRU's from your param trial file
hrud = in_condi.createDimension('hru', hru_num)
# define variables
mlvfi = in_condi.createVariable('mLayerVolFracIce', np.float64, ('midToto', 'hru'))
scat = in_condi.createVariable('scalarCanairTemp', np.float64, ('scalarv', 'hru'))
nsnow = in_condi.createVariable('nSnow', np.int32, ('scalarv', 'hru'))
ilh = in_condi.createVariable('iLayerHeight', np.float64, ('ifcToto', 'hru'))
mlmh = in_condi.createVariable('mLayerMatricHead', np.float64, ('midSoil', 'hru'))
ssa = in_condi.createVariable('scalarSnowAlbedo', np.float64, ('scalarv', 'hru'))
dti = in_condi.createVariable('dt_init', np.float64, ('scalarv', 'hru'))
mlt = in_condi.createVariable('mLayerTemp', np.float64, ('midToto', 'hru'))
ssmp = in_condi.createVariable('scalarSfcMeltPond', np.float64, ('scalarv', 'hru'))
sct = in_condi.createVariable('scalarCanopyTemp', np.float64, ('scalarv', 'hru'))
ssd = in_condi.createVariable('scalarSnowDepth', np.float64, ('scalarv', 'hru'))
nsoil = in_condi.createVariable('nSoil', np.int32, ('scalarv', 'hru'))
sswe = in_condi.createVariable('scalarSWE', np.float64, ('scalarv', 'hru'))
scl = in_condi.createVariable('scalarCanopyLiq', np.float64, ('scalarv', 'hru'))
mlvf = in_condi.createVariable('mLayerVolFracLiq', np.float64, ('midToto', 'hru'))
mld = in_condi.createVariable('mLayerDepth', np.float64, ('midToto', 'hru'))
sci = in_condi.createVariable('scalarCanopyIce', np.float64, ('scalarv', 'hru'))
sas = in_condi.createVariable('scalarAquiferStorage', np.float64, ('scalarv', 'hru'))
#%% # add values for the intial condition variables in HRUs
for varname in ic.variables.keys():
    infovar = ic.variables[varname]
    var = ic.variables[varname][:]
    cic = np.repeat(var[:,np.newaxis], hru_num, axis=1); newic = cic.reshape(infovar.shape[0],hru_num)
    in_condi.variables[varname][:]=newic

print ic.variables['scalarCanairTemp'][:]

scat_sa = np.array([288.])
scat_lat = np.repeat(scat_sa[:,np.newaxis], hru_num, axis=1); scat_lat=scat_lat.reshape(hru_num,)
scat[:] = scat_lat

sct_sa = np.array([292.])
sct_lat = np.repeat(sct_sa[:,np.newaxis], hru_num, axis=1); sct_lat=sct_lat.reshape(hru_num,)
sct[:] = sct_lat

in_condi.close()
#%%
iccheck = Dataset("summa_zInitialCond_scT1_veg.nc")
for varname in iccheck.variables.keys():
    var = iccheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['scalarCanairTemp'][:]







