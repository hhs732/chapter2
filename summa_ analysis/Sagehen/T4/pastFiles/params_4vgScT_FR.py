###       /bin/bash runTestCases_dockerScT1.sh   snwDensity snwDensity
# 2007 - 2008 as wet year for sensirivity analysis 1st step
# I thinking having five scenarios, open, dense and tall, dense and short, sparse and tall, sparse and short. Assuming they are sensitive.
#scalarSWE                 | 1 
#scalarSnowDepth           | 1
#scalarGroundNetNrgFlux    | 1
#scalarGroundAbsorbedSolar | 1
#scalarLWNetGround         | 1
#scalarLatHeatTotal        | 1 
#scalarSenHeatTotal        | 1 
import numpy as np
from netCDF4 import Dataset
import itertools
import pandas as pd

#%% different scenarios
#20,1, 'SHDFAC NROOT   RS      RGL      HS      SNUP  MAXALB   LAIMIN  LAIMAX   EMISSMIN EMISSMAX ALBEDOMIN ALBEDOMAX   Z0MIN    Z0MAX'
#14,     .70,   4,    125.,    30.,   47.35,   0.08,    52.,    5.00,   6.40,   .950,    .950,     .12,      .12,      .50,     .50,     'Evergreen Needleleaf Forest'  
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sc_4veg4.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

p1 = [0.45] #[0.5,0.45,0.5,0.45] #LAIMIN
p2 = [5] #[6,5,6,5] #LAIMAX
p3 = [0.5] #[1,0.5,1,0.5] #winterSAI
p4 = [2] #[5,2,5,2] #summerLAI
p5 = [10] #[25,25,10,10] #heightCanopyTop
p6 = [2] #[5,5,2,2] #heightCanopyBottom 20% of top
p7 = [20] #[40,27,27,20] #maxMassVegetation         |      25.0000 |       1.0000 |      50.0000
#%% # scenario 2 (open space)
p8 = [4.5,6.6,8] #refInterceptCapSnow       |       6.6000 |       1.0000 |      10.0000 #refInterceptCapSnow   =  reference canopy interception capacity per unit leaf area (snow) (kg m-2)
p9 = [0.3,0.45,0.6] #throughfallScaleSnow
p10 = [700,874,950] #specificHeatVeg   j/kg k         |     874.0000 |     500.0000 |    1500.0000
p11 = [0.02,0.04,0.06] #leafDimension             |       0.0400 |       0.0100 |       0.1000

p12 = [0.5] #0.1#z0Canopy                  |       0.1000

#best param so far : #pt12_112112121 # pt12_112212121 #pt22_112112121 #pt2
p13 = [273.66] #273.66 tempCritRain	
p14 = [1.1] # 1.045 frozenPrecipMultip	
p15 = [2] #2, 3, 4] #mw_exp exponent for meltwater flow

p16 = [0.89] #0.89albedoMax |       0.8500 |       0.7000 |       0.9500
p17 = [0.94] #0.89 albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p18 = [0.68] #0.75 albedoMinVisible 0.76|       0.7500 |       0.5000 |       0.7500 ?????????????????????
p19 = [0.75] #albedoMaxNearIR 0.83|       0.6500 |       0.5000 |       0.7500
p20 = [0.45] #albedoMinNearIR  0.49|       0.3000 |       0.1500 |       0.4500
p21 = [0.3] #0.5albedoSootLoad
p22 = [3]#,1] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p23 = [100000]#,200000,400000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p24 = [0.65] #albedoMinWinter           |       0.6500 |       0.6000 |       1.0000
p25 = [0.5] #albedoMinSpring           |       0.5000 |       0.3000 |       1.0000

p26 = [60] #newSnowDenMin 
p27= [75] #newSnowDenMult            |      75.0000 |      25.0000 |      75.0000

def hru_ix_ID(p1, p2, p3, p4):#, p5, p6, p7, p8, p9, p10, p11):#, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21):
    ix1 = np.arange(1,len(p1)+1)
    ix2 = np.arange(1,len(p2)+1)
    ix3 = np.arange(1,len(p3)+1)
    ix4 = np.arange(1,len(p4)+1)

    c = list(itertools.product(ix1,ix2,ix3,ix4))#,,ix5,ix6,ix7,ix8,ix9,ix10,ix11ix12,ix13,ix14,ix15,ix16,ix17,ix18,ix19,ix20,ix21))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p8, p9, p10, p11)#,, p5, p6, p7, p8, p9, p10, p11 p12, p13, p14, p15, p16, p17, p18, p19, p20, p21)
#
hru_num = np.size(hruidxID)
#%% function to create lists of each parameter, this will iterate through to make sure all combinations are covered
def param_fill(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27): 
    b = list(itertools.product(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20 , p21, p22, p23, p24, p25, p26, p27))#
    p1l =[]; p2l =[]; p3l =[]; p4l=[]; p5l =[]; p6l =[]; p7l =[]; p8l =[]; p9l =[]; p10l=[]; p11l =[]; p12l =[]; p13l =[]; p14l=[]; p15l = []; 
    p16l =[]; p17l=[]; p18l = []; p19l = []; p20l = []; p21l=[]; p22l=[]; p23l=[]; p24l=[]; p25l=[]; p26l=[]; p27l=[]#
    for tup in b:
        p1l.append(tup[0]); p2l.append(tup[1]); p3l.append(tup[2]); p4l.append(tup[3]); p5l.append(tup[4]); p6l.append(tup[5]); p7l.append(tup[6]); 
        p8l.append(tup[7]); p9l.append(tup[8]); p10l.append(tup[8]); p11l.append(tup[10]); p12l.append(tup[11]); p13l.append(tup[12]);
        p14l.append(tup[13]); p15l.append(tup[14]); p16l.append(tup[15]); p17l.append(tup[16]); p18l.append(tup[17]); p19l.append(tup[18]); 
        p20l.append(tup[19]); p21l.append(tup[20]); p22l.append(tup[21]); p23l.append(tup[22]); p24l.append(tup[23]); p25l.append(tup[24]);
        p26l.append(tup[25]); p27l.append(tup[26])
    return(p1l, p2l, p3l, p4l, p5l, p6l, p7l, p8l, p9l, p10l, p11l, p12l, p13l, p14l, p15l, p16l, p17l, p18l, p19l, p20l, p21l, p22l, p23l, p24l, p25l, p26l, p27l)#  

# call the function on the parameters
valst1 = param_fill(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27) 

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

hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['LAIMIN','LAIMAX','winterSAI','summerLAI','heightCanopyTop','heightCanopyBottom','maxMassVegetation',
                  'refInterceptCapSnow','throughfallScaleSnow','specificHeatVeg','leafDimension','z0Canopy',
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
    #print varname, var.dtype, var.dimensions, var.shape

#print paramfile.variables['albedoSootLoad'][:]
paramfile.close()

#%% 
varcheck = Dataset ('summa_zParamTrial_variableDecayRate_sc_4veg1.nc')
print varcheck.variables['hruIndex'][:]
print varcheck.variables['maxMassVegetation'][:]

#%% # local attributes file
local_atrbt = Dataset("summa_zLocalAttributes_sagehenCreekT4_4veg.nc",'w',format='NETCDF3_CLASSIC')
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
#%% T4 lat: 39.42222°  lon: 120.2989° elev:2370; T1 [ 39.4321] [-120.2411] elev = 1936m 
vegindx_sc = np.array([14])
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

local_atrbt.variables['hruId'][:] = varcheck.variables['hruIndex'][:]

#print local_atrbt.variables['hruId'][:]
local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_sagehenCreekT4_4veg.nc')
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

in_condi = Dataset("summa_zInitialCond_sc_4veg.nc",'w',format='NETCDF3_CLASSIC')
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
iccheck = Dataset("summa_zInitialCond_sc_4veg.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['mLayerVolFracLiq'][:]




























