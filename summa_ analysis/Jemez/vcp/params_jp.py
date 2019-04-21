###       /bin/bash runTestCases_dockerSC.sh   snwDensity snwDensity
# 2007 - 2008 as wet year for sensirivity analysis 1st step
# I thinking having five scenarios, open, dense and tall, dense and short, sparse and tall, sparse and short. Assuming they are sensitive.
import numpy as np
from netCDF4 import Dataset
import itertools
import pandas as pd

#%% different scenarios
#20,1, 'SHDFAC NROOT   RS      RGL      HS      SNUP  MAXALB   LAIMIN  LAIMAX   EMISSMIN EMISSMAX ALBEDOMIN ALBEDOMAX   Z0MIN    Z0MAX'
#14,     .70,   4,    125.,    30.,   47.35,   0.08,    52.,    5.00,   6.40,   .950,    .950,     .12,      .12,      .50,     .50,     'Evergreen Needleleaf Forest'  

#p1 = [0.5,0.45,0.5,0.45] #LAIMIN
#p2 = [5,5,5,5] #LAIMAX
#p3 = [1,0.5,1,0.5] #winterSAI
#p4 = [5,2,5,2] #summerLAI
#p7 = [25,25,10,10] #heightCanopyTop
#p8 = [3,3,2.5,2.5] #heightCanopyBottom 20% of top
#p11 = [35,20,30,20] #maxMassVegetation         |      25.0000 |       1.0000 |      50.0000
#%% # scenario 1 (open space)
#p1 = [0.1] #LAIMIN
#p2 = [1] #LAIMAX
p3 = [0.01] #winterSAI
p4 = [0.01] #summerLAI
p5 = [0.5] #heightCanopyTop
p6 = [0.1] #heightCanopyBottom
#p7 = [25] #maxMassVegetation 

#p8 = [60] #newSnowDenMin 
p9 = [100000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p10 = [0.89] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p11 = [0.76] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p12 = [0.83] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p13 = [0.49] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500
p14 = [0.5] #albedoSootLoad
#p15 = [1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000

p16 = [273.66] #tempCritRain	
p17 = [1.045] #frozenPrecipMultip	

#p18 = [3] #2, 3, 4] #mw_exp exponent for meltwater flow
#p19 = [0.35] #0.2, 0.4 , 0.6] #fixedThermalCond_snow used in tcs_smnv model, but we are using Jordan model  

#p20 = [0.001] #z0Snow

#p8 = [2] #rootingDepth
#p21 = [6.6] #refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)
#p22 = [0.89] #throughfallScaleSnow
#p23 = [0.4] #ratioDrip2Unloading       |       0.4000 |       0.0000 |       1.0000
#p24 = [1] #rootDistExp |       1.0000 |       0.0100 |       1.0000
#p25 = [874] #specificHeatVeg   j/kg k         |     874.0000 |     500.0000 |    1500.0000
#p26 = [0.04] #leafDimension             |       0.0400 |       0.0100 |       0.1000
#p27 = [0.02] #z0Canopy                  |       0.0200 |       0.0010 |      10.0000

#p28 = [0.28] #windReductionParam        |       0.2800 |       0.0000 |       1.0000
#p29 = [0.2] #critRichNumber
#p30 = [1] #Mahrt87_eScale  
#p31 = [0.06] #Fcapil
#p32 = [0.015] #k_snow
#p33 = [9.4] #Louis79_bparam  
#p34 = [5.3] #Louis79_cStar 

list_param = pd.DataFrame([p3,p4,p5,p6,p9,p10,p11,p12,p13,p14,p16,p17])
hruidxID = list(np.arange(101,102))
hru_num = np.size(hruidxID)
#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
paramfile = Dataset("summa_zParamTrial_variableDecayRate_scT1.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file
 
hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['winterSAI','summerLAI','heightCanopyTop','heightCanopyBottom',#'maxMassVegetation','refInterceptCapSnow','ratioDrip2Unloading','critRichNumber',
                  #'rootingDepth','rootDistExp','throughfallScaleSnow','specificHeatVeg','leafDimension','newSnowDenMin',
                  'albedoDecayRate','albedoMaxVisible','albedoMinVisible','albedoMaxNearIR','albedoMinNearIR','albedoSootLoad',#'albedoRefresh',
                  #'z0Snow','z0Canopy','windReductionParam','mw_exp','fixedThermalCond_snow','Mahrt87_eScale','Fcapil','k_snow','Louis79_bparam','Louis79_cStar'
                  'tempCritRain','frozenPrecipMultip'] 
for param in param_nam_list:
    paramfile.createVariable(param, np.float64,('hru',))

constant_params = ['rootingDepth','rootDistExp','theta_sat','theta_res','vGn_alpha','vGn_n','k_soil','critSoilWilting','critSoilTranspire']
for params in constant_params:
    paramfile.createVariable(params, np.float64,('hru',))
#paramfile.close()
#%% parameterTrial, Local attributes and initial conditions for senatore beck
pt = Dataset('summa_zParamTrial_variableDecayRate.nc')
la = Dataset('summa_zLocalAttributes_senatorSheltered.nc') #('settings/wrrPaperTestCases/figure07/summa_zLocalAttributes_riparianAspen.nc')
ic = Dataset('summa_zInitialCond.nc') #('settings/wrrPaperTestCases/figure07/summa_zInitialCond.nc')
for j in pt.variables:
    print j
print pt.variables['winterSAI'][:] #0.45
print pt.variables['rootingDepth'][:] #1
print pt.variables['summerLAI'][:] #3
print pt.variables['heightCanopyTop'][:] #0.5
print pt.variables['heightCanopyBottom'][:] #0.05
print la.variables['vegTypeIndex'][:] #7 Ava:7
print la.variables['soilTypeIndex'][:] #8 Ava:9
print la.variables['elevation'][:] 

#%% # add values for the constant variables in HRUs for parameter Trail file
for varname in pt.variables.keys():
    var = pt.variables[varname][0]
    c = np.full((hru_num,),var)
    try :
        paramfile.variables[varname][:]=c
    except IndexError: # size of data array does not conform to slice
        pass
#%% creating changing variables and adding values for changing variables
#paramfile.variables['LAIMIN'][:]=p1
#paramfile.variables['LAIMAX'][:]=p2
paramfile.variables['winterSAI'][:]=p3
paramfile.variables['summerLAI'][:]=p4
paramfile.variables['heightCanopyTop'][:]=p5
paramfile.variables['heightCanopyBottom'][:]=p6

paramfile.variables['albedoDecayRate'][:]=p9
paramfile.variables['albedoMaxVisible'][:]=p10
paramfile.variables['albedoMinVisible'][:]=p11
paramfile.variables['albedoMaxNearIR'][:]=p12
paramfile.variables['albedoMinNearIR'][:]=p13
paramfile.variables['albedoSootLoad'][:]=p14

paramfile.variables['tempCritRain'][:]=p16
paramfile.variables['frozenPrecipMultip'][:]=p17

#paramfile.variables['maxMassVegetation'][:]=p5
#paramfile.variables['refInterceptCapSnow'][:]=p6
#paramfile.variables['ratioDrip2Unloading'][:]=p7
#paramfile.variables['critRichNumber'][:]=p8

#paramfile.variables['rootingDepth'][:]=p9
#paramfile.variables['rootDistExp'][:]=p10
#paramfile.variables['throughfallScaleSnow'][:]=p11
#paramfile.variables['specificHeatVeg'][:]=p12
#paramfile.variables['leafDimension'][:]=p13
#
#paramfile.variables['newSnowDenMin'][:]=p14
#paramfile.variables['z0Snow'][:]=p21
#paramfile.variables['z0Canopy'][:]=p22
#paramfile.variables['windReductionParam'][:]=p23
#paramfile.variables['mw_exp'][:]=p24
#paramfile.variables['fixedThermalCond_snow'][:]=p25
#
#paramfile.variables['Mahrt87_eScale'][:]=p26
#paramfile.variables['Fcapil'][:]=p27
#paramfile.variables['k_snow'][:]=p28
#paramfile.variables['Louis79_bparam'][:]=p29
#paramfile.variables['Louis79_cStar'][:]=p30

paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    print varname, var.dtype, var.dimensions, var.shape

#print paramfile.variables['hruIndex'][:]
paramfile.close()
#%% 
varcheck = Dataset ('summa_zParamTrial_variableDecayRate_scT1.nc')
print varcheck.variables['rootingDepth'][:]
print varcheck.variables['albedoSootLoad'][:]
#%% # local attributes file
local_atrbt = Dataset("summa_zLocalAttributes_sagehenCreekT1.nc",'w',format='NETCDF3_CLASSIC')
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
#%% add values for the changing variables in HRUs for local attribute file
scFD = Dataset('shT1_force1_2.nc')

lat_sc = np.array(scFD.variables['latitude'][:])
len_lat = np.repeat(lat_sc[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)
long_sc = np.array(scFD.variables['longitude'][:])
len_lon= np.repeat(long_sc[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
lat[:] = len_lat
lon[:] = len_lon

#vegindx_sc = np.array([14])
#vegindx_lon = np.repeat(vegindx_sc[:,np.newaxis], hru_num, axis=1); vegindx_lon=vegindx_lon.reshape(hru_num,)
#vegindx[:] = vegindx_lon

mHeight_sc = np.array([7.62]) #mHeight
mHeight_len = np.repeat(mHeight_sc[:,np.newaxis], hru_num, axis=1); mHeight_len=mHeight_len.reshape(hru_num,)
mh[:] = mHeight_len

elev_sc = np.array([1936]) #mHeight
elev_len = np.repeat(elev_sc[:,np.newaxis], hru_num, axis=1); elev_len=elev_len.reshape(hru_num,)
elev[:] = elev_len

#%% # get the hru, gru, and hru2gru in local_atribute file
newgru = np.array([1111])
local_atrbt.variables['gruId'][:] = newgru

c3 = np.repeat(newgru[:,np.newaxis], hru_num, axis=1); newlad = c3.reshape(hru_num,)
local_atrbt.variables['hru2gruId'][:] = c3

local_atrbt.variables['hruId'][:] = hruidxID

#print local_atrbt.variables['hruId'][:]
local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_sagehenCreekT1.nc')
print lacheck.variables['vegTypeIndex'][:]
print lacheck.variables['latitude'][:]
print lacheck.variables['mHeight'][:]
print lacheck.variables['elevation'][:]

#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
#lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_scT1.nc",'w',format='NETCDF3_CLASSIC')
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

print in_condi.variables['iLayerHeight'][:]

in_condi.close()
#%%
iccheck = Dataset("summa_zInitialCond_scT1.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['mLayerVolFracLiq'][:]







