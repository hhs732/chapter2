###       /bin/bash runTestCases_dockerSC.sh   snwDensity snwDensity
# 2007 - 2008 as wet year for sensirivity analysis 1st step
# I thinking having five scenarios, open, dense and tall, dense and short, sparse and tall, sparse and short. Assuming they are sensitive.
import numpy as np
from netCDF4 import Dataset
import itertools
import pandas as pd

#%% different scenarios #best param so far : #pt12_112112121 # pt12_112212121 #pt22_112112121 #pt2

#p1 = [0.1] #LAIMIN
#p2 = [1] #LAIMAX
p3 = [0.01] #winterSAI
p4 = [0.1] #summerLAI
p5 = [0.3] #heightCanopyTop
p6 = [0.03] #heightCanopyBottom to calculate wind speed, wind speed reduction; coupute roughness length of the veg canopy; neutral ground resistance; canopy air temp;
p7 = [0.01] #0.1#z0Canopy                  |       0.1000
p8 = [6.6] #refInterceptCapSnow 

p9 = [273.66] #273.66 tempCritRain	
p10 = [1.1] # 1.045 frozenPrecipMultip	
p11 = [2] #2, 3, 4] #mw_exp exponent for meltwater flow

p12 = [0.89] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500 0.87
p13 = [0.94] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p14 = [0.68] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500 ?????????????????????
p15 = [0.75] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p16 = [0.45] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p17 = [0.3] #0.5albedoSootLoad
p18 = [3]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p19 = [100000]#,200000,400000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p20 = [0.65] #albedoMinWinter           |       0.6500 |       0.6000 |       1.0000
p21 = [0.5] #albedoMinSpring           |       0.5000 |       0.3000 |       1.0000

p22 = [60] #newSnowDenMin 
p23= [75] #newSnowDenMult            |      75.0000 |      25.0000 |      75.0000

hruidxID = [101]
hru_num = np.size(hruidxID)
#%% function to create lists of each parameter, this will iterate through to make sure all combinations are covered
def param_fill(p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23): 
    b = list(itertools.product(p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20 , p21, p22, p23))#
    p3l =[]; p4l=[]; p5l =[]; p6l =[]; p7l =[]; p8l =[]; p9l =[]; p10l=[]; p11l =[]; p12l =[]; p13l =[]; p14l=[]; p15l = []; 
    p16l =[]; p17l=[]; p18l = []; p19l = []; p20l = []; p21l=[]; p22l = []; p23l=[]
    for tup in b:
        p3l.append(tup[0]); p4l.append(tup[1]); p5l.append(tup[2]); p6l.append(tup[3]); p7l.append(tup[4]); 
        p8l.append(tup[5]); p9l.append(tup[6]); p10l.append(tup[7]); p11l.append(tup[8]); p12l.append(tup[9]); p13l.append(tup[10]);
        p14l.append(tup[11]); p15l.append(tup[12]); p16l.append(tup[13]); p17l.append(tup[14]); p18l.append(tup[15]); p19l.append(tup[16]); 
        p20l.append(tup[17]); p21l.append(tup[18]); p22l.append(tup[19]); p23l.append(tup[20])
    return(p3l, p4l, p5l, p6l, p7l, p8l, p9l, p10l, p11l, p12l, p13l, p14l, p15l, p16l, p17l, p18l, p19l, p20l, p21l, p22l, p23l)#  

# call the function on the parameters
valst1 = param_fill(p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23) 

#%% parameterTrial, Local attributes and initial conditions for senatore beck
pt = Dataset('summa_zParamTrial_variableDecayRate.nc')
la = Dataset("NewData_ava\T4_2016\summa_zParamTrial_T4_16.nc") #('settings/wrrPaperTestCases/figure07/summa_zLocalAttributes_riparianAspen.nc')
ic = Dataset('summa_zInitialCond.nc') #('settings/wrrPaperTestCases/figure07/summa_zInitialCond.nc')

paramfile_in =Dataset("NewData_ava\T1_2016\summa_zParamTrial_T1_16.nc")
la_in = Dataset("NewData_ava\T1_2016\SH_T1att.nc")

for j in la_in.variables:
    print j

la_in.variables['mHeight'][:]

#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
paramfile = Dataset("summa_zParamTrial_variableDecayRate_open.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['winterSAI','summerLAI','heightCanopyTop','heightCanopyBottom','z0Canopy','refInterceptCapSnow',
                  #'LAIMIN','LAIMAX','maxMassVegetation','throughfallScaleSnow','specificHeatVeg','leafDimension',
                  'tempCritRain','frozenPrecipMultip','mw_exp',
                  'albedoMax','albedoMaxVisible','albedoMinVisible','albedoMaxNearIR','albedoMinNearIR',
                  'albedoSootLoad','albedoRefresh','albedoDecayRate','albedoMinWinter','albedoMinSpring',
                  'newSnowDenMin','newSnowDenMult'] 
for param in param_nam_list:
    paramfile.createVariable(param, np.float64,('hru',))

constant_params = ['rootDistExp','critSoilWilting','critSoilTranspire','kAnisotropic','zScale_TOPMODEL','qSurfScale','theta_mp','theta_sat','theta_res',
                   'vGn_alpha','vGn_n','f_impede','k_soil','k_macropore']
for params in constant_params:
    paramfile.createVariable(params, np.float64,('hru',))
#%% # add values for the constant variables in HRUs for parameter Trail file
for varname in paramfile_in.variables.keys():
    var = paramfile_in.variables[varname][0]
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
# don't forget the HRU Index!!
paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    print varname, var.dtype, var.dimensions, var.shape

#print paramfile.variables['albedoSootLoad'][:]
paramfile.close()

#%% 
varcheck = Dataset ('summa_zParamTrial_variableDecayRate_open.nc')
print varcheck.variables['albedoMinSpring'][:]

#%% # local attributes file
local_atrbt = Dataset("summa_zLocalAttributes_scT4_open.nc",'w',format='NETCDF3_CLASSIC')
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
for varname in la_in.variables.keys():
    var = la_in.variables[varname][0]
    #print var
    c2 = np.full((hru_num,),var)
    #print c2
    try :
        local_atrbt.variables[varname][:]=c2
    except IndexError: # size of data array does not conform to slice
        pass
    #local_atrbt.variables[varname][:]=c2
#%% add values for the changing variables in HRUs for local attribute file
#scFD = Dataset('shT1_force1_2.nc')
#T4 lat: 39.42222°  lon: 120.2989° elev:2370; T1 [ 39.4321] [-120.2411] elev = 1936m 
vegindx_sc = np.array([7])
vegindx_lon = np.repeat(vegindx_sc[:,np.newaxis], hru_num, axis=1); vegindx_lon=vegindx_lon.reshape(hru_num,)
vegindx[:] = vegindx_lon

lat_sa = np.array([39.4222])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)
lat[:] = len_lat

long_sa = np.array([-120.2989])
len_lon = np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
lon[:] = len_lon

elev_sa = np.array([2370])
elev_lon = np.repeat(elev_sa[:,np.newaxis], hru_num, axis=1); elev_lon=elev_lon.reshape(hru_num,)
elev[:] = elev_lon
#%% # get the hru, gru, and hru2gru in local_atribute file
newgru = np.array([11111])
local_atrbt.variables['gruId'][:] = newgru

c3 = np.repeat(newgru[:,np.newaxis], hru_num, axis=1); newlad = c3.reshape(hru_num,)
local_atrbt.variables['hru2gruId'][:] = c3

local_atrbt.variables['hruId'][:] = hruidxID

#print local_atrbt.variables['hruId'][:]
local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_scT4_open.nc')
print lacheck.variables['vegTypeIndex'][:]
print lacheck.variables['longitude'][:]
print lacheck.variables['latitude'][:]
print lacheck.variables['mHeight'][:]
print lacheck.variables['elevation'][:]
print lacheck.variables['hruId'][:]


#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
#lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_open.nc",'w',format='NETCDF3_CLASSIC')
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
iccheck = Dataset("summa_zInitialCond_open.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['mLayerVolFracLiq'][:]

#%%
#ckeckparam = Dataset ("C:\Users\HHS\summaTestCases_2.x\settings\sagehencreek\summa_zParamTrial_variableDecayRate_scT1_4vegSc.nc")
ckeckla = Dataset ("C:\Users\HHS\summaTestCases_2.x\settings\sagehencreek\summa_zLocalAttributes_scT1_open.nc")
#
#for varname in ckeckparam.variables.keys():
#    var = ckeckparam.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
#print ckeckparam.variables['hruIndex'][:]
#
#for varname in ckeckla.variables.keys():
#    var = ckeckla.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
#print ckeckla.variables['hruId'][:]
#print ckeckla.variables['vegTypeIndex'][:]
#print ckeckla.variables['longitude'][:]
#print ckeckla.variables['latitude'][:]
print ckeckla.variables['mHeight'][:]
#print ckeckla.variables['elevation'][:]


























