###       /bin/bash runTestCases_docker.sh   snwDensity snwDensity
# 2007 - 2008 as wet year for sensirivity analysis 1st step
import numpy as np
from netCDF4 import Dataset
import itertools
#%%  all parameters
#20,1, 'SHDFAC NROOT   RS      RGL      HS      SNUP  MAXALB   LAIMIN  LAIMAX   EMISSMIN EMISSMAX ALBEDOMIN ALBEDOMAX   Z0MIN    Z0MAX'
#14,     .70,   4,    125.,    30.,   47.35,   0.08,    52.,    5.00,   6.40,   .950,    .950,     .12,      .12,      .50,     .50,     'Evergreen Needleleaf Forest'  
p1 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]#1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]# #LAIMIN
p2 = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]#5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]# #LAIMAX
p3 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]#0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #winterSAI
p4 = [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]#2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2] #summerLAI

p5 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2] #rootingDepth
p6 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20] #heightCanopyTop
p7 = [0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05] #heightCanopyBottom
p8 = [60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60] #newSnowDenMin 

p9 = [0.5,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #throughfallScaleSnow
p10 = [800,800,1000,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1000,800,800,800,800,800,800,800,800,800,800,800,800,800] #specificHeatVeg           |     874.0000 |     500.0000 |    1500.0000
p11 = [30,30,30,20,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,20,30,30,30,30,30,30,30,30,30,30,30,30] #maxMassVegetation         |      25.0000 |       1.0000 |      50.0000
p12 = [6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6] #refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000

p13 = [1000000,1000000,1000000,1000000,1000000,700000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,700000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p14 = [0.9,0.9,0.9,0.9,0.9,0.9,0.94,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.94,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p15 = [0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.7,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.7,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p16 = [0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.6,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.6,0.65,0.65,0.65,0.65,0.65,0.65,0.65] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p17 = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.35,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.35,0.3,0.3,0.3,0.3,0.3,0.3] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500
p18 = [3,3,3,3,3,3,3,3,3,3,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,3,3,3,3,3] # 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000

p19 = [0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.001,0.001] #[0.001, 0.002] #z0Snow
p20 = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01] #z0Canopy                  |       0.0200 |       0.0010 |      10.0000
p21 = [0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.4,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.4,0.28,0.28,0.28]	#windReductionParam        |       0.2800 |       0.0000 |       1.0000
p22 = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3] #2, 3, 4] #mw_exp exponent for meltwater flow
p23 = [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.5,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.5] #0.2, 0.4 , 0.6] #fixedThermalCond_snow

#p21 = [0.700, 1.000, 1.500] #Mahrt87_eScale  
#p14 = [0.040, 0.060, 0.080] #Fcapil
#p15 = [0.100, 0.015, 0.350] #k_snow
#p18 = [0.150, 0.200, 0.400] #critRichNumber  
#p19 = [9.300, 9.400, 9.500] #Louis79_bparam  
#p20 = [5.200, 5.300, 5.400] #Louis79_cStar 
#p15 = [105] #51, 70, 105] #constSnowDen 70.00, 100.0, 170.0    55 56 57 

hruidxID = list(np.arange(101,133))
hru_num = np.size(hruidxID)
#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sc.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['LAIMIN','LAIMAX','winterSAI','summerLAI','rootingDepth','heightCanopyTop','heightCanopyBottom','newSnowDenMin',
                  'throughfallScaleSnow','specificHeatVeg','maxMassVegetation','refInterceptCapSnow',
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR', 'albedoRefresh', 
                  'z0Snow', 'z0Canopy', 'windReductionParam', 'mw_exp', 'fixedThermalCond_snow'] 

for param in param_nam_list:
    paramfile.createVariable(param, np.float64,('hru',))

constant_params = ['frozenPrecipMultip','rootDistExp','theta_sat','theta_res','vGn_alpha','vGn_n','k_soil','critSoilWilting','critSoilTranspire']
for params in constant_params:
    paramfile.createVariable(params, np.float64,('hru',))
#paramfile.close()
#%% parameterTrial, Local attributes and initial conditions for senatore beck
pt = Dataset('summa_zParamTrial_variableDecayRate.nc')
la = Dataset('summa_zLocalAttributes_senatorSheltered.nc') #('settings/wrrPaperTestCases/figure07/summa_zLocalAttributes_riparianAspen.nc')
ic = Dataset('summa_zInitialCond.nc') #('settings/wrrPaperTestCases/figure07/summa_zInitialCond.nc')
for j in ic.variables:
    print j
#%% # add values for the constant variables in HRUs for parameter Trail file
for varname in pt.variables.keys():
    var = pt.variables[varname][0]
    c = np.full((hru_num,),var)
    try :
        paramfile.variables[varname][:]=c
    except IndexError: # size of data array does not conform to slice
        pass
#%% creating changing variables and adding values
paramfile.variables['LAIMIN'][:]=p1
paramfile.variables['LAIMAX'][:]=p2
paramfile.variables['winterSAI'][:]=p3
paramfile.variables['summerLAI'][:]=p4
paramfile.variables['rootingDepth'][:]=p5
paramfile.variables['heightCanopyTop'][:]=p6
paramfile.variables['heightCanopyBottom'][:]=p7
paramfile.variables['newSnowDenMin'][:]=p8

paramfile.variables['throughfallScaleSnow'][:]=p9
paramfile.variables['specificHeatVeg'][:]=p10
paramfile.variables['maxMassVegetation'][:]=p11
paramfile.variables['refInterceptCapSnow'][:]=p12

paramfile.variables['albedoDecayRate'][:]=p13
paramfile.variables['albedoMaxVisible'][:]=p14
paramfile.variables['albedoMinVisible'][:]=p15
paramfile.variables['albedoMaxNearIR'][:]=p16
paramfile.variables['albedoMinNearIR'][:]=p17
paramfile.variables['albedoRefresh'][:]=p18

paramfile.variables['z0Snow'][:]=p19
paramfile.variables['z0Canopy'][:]=p20
paramfile.variables['windReductionParam'][:]=p21
paramfile.variables['mw_exp'][:]=p22
paramfile.variables['fixedThermalCond_snow'][:]=p23

paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    print varname, var.dtype, var.dimensions, var.shape

#print paramfile.variables['hruIndex'][:]
paramfile.close()
#%% 
varcheck = Dataset ('summa_zParamTrial_variableDecayRate_sc.nc')
print varcheck.variables['theta_res'][:]
check2 =  varcheck.variables['albedoMaxNearIR'][:]
#%% # local attributes file
local_atrbt = Dataset("summa_zLocalAttributes_sagehenCreek.nc",'w',format='NETCDF3_CLASSIC')
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
#%% add values for the changing variables in HRUs for local atribute file
scFD = Dataset('shT1_force1_2.nc')

lat_sc = np.array(scFD.variables['latitude'][:])
len_lat = np.repeat(lat_sc[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)
long_sc = np.array(scFD.variables['longitude'][:])
len_lon= np.repeat(long_sc[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
lat[:] = len_lat
lon[:] = len_lon

vegindx_sc = np.array([14])
vegindx_lon= np.repeat(vegindx_sc[:,np.newaxis], hru_num, axis=1); vegindx_lon=vegindx_lon.reshape(hru_num,)
vegindx[:] = vegindx_lon
#%% # get the hru, gru, and hru2gru in local_atribute file
newgru = np.array([1111])
local_atrbt.variables['gruId'][:] = newgru

c3 = np.repeat(newgru[:,np.newaxis], hru_num, axis=1); newlad = c3.reshape(hru_num,)
local_atrbt.variables['hru2gruId'][:] = c3

local_atrbt.variables['hruId'][:] = hruidxID

#print local_atrbt.variables['hruId'][:]
local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_sagehenCreek.nc')
print lacheck.variables['vegTypeIndex'][:]
print lacheck.variables['latitude'][:]

#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
#lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_sc.nc",'w',format='NETCDF3_CLASSIC')
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
iccheck = Dataset("summa_zInitialCond_sc.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['mLayerVolFracLiq'][:]







