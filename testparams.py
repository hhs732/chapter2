from netCDF4 import Dataset

pt1 = Dataset('C:/Users/HHS/summaTestCases_2.x/settings/wrrPaperTestCases/figure05/summa_zParamTrial_hedpomPP.nc')
pt2 = Dataset("C:/Users/HHS/summaTestCases_2.x/settings/wrrPaperTestCases/figure05/summa_zParamTrial_storckPP.nc") 
pt0 = Dataset("C:/Users/HHS/summaTestCases_2.x/settings/wrrPaperTestCases/figure06/summa_zParamTrial_baseSheltered.nc")
pt01 = Dataset("C:/Users/HHS/summaTestCases_2.x/settings/wrrPaperTestCases/figure06/summa_zParamTrial_variableDecayRate_test.nc") 

for varname in pt01.variables.keys():
    var = pt01.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    

pt1List = [pt1.variables['tempCritRain'][:],pt1.variables['tempRangeTimestep'][:],pt1.variables['heightCanopyTop'][:],
           pt1.variables['heightCanopyBottom'][:],pt1.variables['refInterceptCapSnow'][:],pt1.variables['snowUnloadingCoeff'][:],
           pt1.variables['ratioDrip2Unloading'][:]]

pt2List = [pt2.variables['tempCritRain'][:],pt2.variables['tempRangeTimestep'][:],pt2.variables['heightCanopyTop'][:],
           pt2.variables['heightCanopyBottom'][:],pt2.variables['refInterceptCapSnow'][:],pt2.variables['snowUnloadingCoeff'][:],
           pt2.variables['ratioDrip2Unloading'][:]]

pt0List = [pt0.variables['frozenPrecipMultip'][:],pt0.variables['rootingDepth'][:],pt0.variables['rootDistExp'][:],
           pt0.variables['theta_sat'][:],pt0.variables['theta_res'][:],pt0.variables['vGn_alpha'][:],
           pt0.variables['vGn_n'][:],pt0.variables['k_soil'][:],pt0.variables['critSoilWilting'][:],
           pt0.variables['critSoilTranspire'][:],pt0.variables['winterSAI'][:],pt0.variables['summerLAI'][:],
           pt0.variables['heightCanopyTop'][:],pt0.variables['heightCanopyBottom'][:],pt0.variables['albedoDecayRate'][:]]

pt01List = [pt01.variables['frozenPrecipMultip'][:],pt01.variables['rootingDepth'][:],pt01.variables['rootDistExp'][:],
           pt01.variables['theta_sat'][:],pt01.variables['theta_res'][:],pt01.variables['vGn_alpha'][:],
           pt01.variables['vGn_n'][:],pt01.variables['k_soil'][:],pt01.variables['critSoilWilting'][:],
           pt01.variables['critSoilTranspire'][:],pt01.variables['winterSAI'][:],pt01.variables['summerLAI'][:],
           pt01.variables['heightCanopyTop'][:],pt01.variables['heightCanopyBottom'][:],pt01.variables['albedoDecayRate'][:]]































