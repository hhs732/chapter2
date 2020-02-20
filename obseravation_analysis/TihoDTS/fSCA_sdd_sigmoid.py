import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy import stats
from scipy import optimize

with open("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/TihoDTS/fSCA-SDD2.csv") as safd:
    reader = csv.reader(safd)
    data_forcing = [r for r in reader]
data_forcing2 = data_forcing[1:]
sa_fd_column = []
for csv_counter1 in range (len (data_forcing2)):
    for csv_counter2 in range (30):
        sa_fd_column.append(float(data_forcing2[csv_counter1][csv_counter2]))
sdd_fsca=np.reshape(sa_fd_column,(len (data_forcing2),30))
sdd_obs = sdd_fsca[:,0]

#slope, intercept, r_value, p_value, std_err = stats.linregress(sdd_fsca[:,0], sdd_fsca[:,2])

def fit_func(x, a, b, c, d):
    return a/(b+c*np.exp(-d*x))
p0101213=[0.1,-10.21,100.21,-0.05]
p01527=[0.03,-20.21,10.21,-0.01]
params20, params_covariance20 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,16],
                                                   p0=[0.03,-20.21,10.21,-0.01])

plt.figure(figsize=(15, 10))
plt.scatter(sdd_obs, sdd_fsca[:,10], color='purple', s = 20**2, label='Data')
fsca_sigmoid = fit_func(np.array([70,82.,90,100,110,120,130]),params20[0]*1,
                        params20[1]*1,params20[2]*1,params20[3]*1)
plt.plot(np.array([70, 82.,90,100,110,120,130]), fsca_sigmoid,linewidth= 8,color = 'black', label='Fitted function')

params_sigmoid10 = []
covariance_sigmoid10 =[]
for fsca in range (9):#len(sdd_fsca[0,:])-1
    params, params_covariance = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,1+fsca],
                                                   p0=[0.03,-20.21,10.21,-0.01])
    params_sigmoid10.append(params)
    covariance_sigmoid10.append(params_covariance)
    

p01527=[0.03,-20.21,10.21,-0.01]
params_sigmoid102, params_covariance10 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,10],
                                                   p0=[0.1,-10.21,100.21,-0.05])
params_sigmoid12, params_covariance12 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,12],
                                                    p0=[0.1,-10.21,100.21,-0.05])
params_sigmoid13, params_covariance13 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,13],
                                                    p0=[0.1,-10.21,100.21,-0.05])
params_sigmoid15_27 = []
covariance_sigmoid15_27 =[]
for fsca in range (13):#len(sdd_fsca[0,:])-1
    params, params_covariance = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,15+fsca],
                                                   p0=[0.03,-20.21,10.21,-0.01])
    params_sigmoid15_27.append(params)
    covariance_sigmoid15_27.append(params_covariance)

x_sdd = np.array([60.,65.,70.,75.,80.,85.,90.,95.,100.,105.,110.,115.,120.,125.,130.])
fsca_sigmoid1 = fit_func(x_sdd,params_sigmoid10[0][0],params_sigmoid10[0][1],params_sigmoid10[0][2],
                         params_sigmoid10[0][3])
fsca_sigmoid2 = fit_func(x_sdd,params_sigmoid10[1][0],params_sigmoid10[1][1],params_sigmoid10[1][2],
                         params_sigmoid10[1][3])
fsca_sigmoid3 = fit_func(x_sdd,params_sigmoid10[2][0],params_sigmoid10[2][1],params_sigmoid10[2][2],
                         params_sigmoid10[2][3])
fsca_sigmoid4 = fit_func(x_sdd,params_sigmoid10[3][0],params_sigmoid10[3][1],params_sigmoid10[3][2],
                         params_sigmoid10[3][3])
fsca_sigmoid5 = fit_func(x_sdd,params_sigmoid10[4][0],params_sigmoid10[4][1],params_sigmoid10[4][2],
                         params_sigmoid10[4][3])
fsca_sigmoid6 = fit_func(x_sdd,params_sigmoid10[5][0],params_sigmoid10[5][1],params_sigmoid10[5][2],
                         params_sigmoid10[5][3])
fsca_sigmoid7 = fit_func(x_sdd,params_sigmoid10[6][0],params_sigmoid10[6][1],params_sigmoid10[6][2],
                         params_sigmoid10[6][3])
fsca_sigmoid8 = fit_func(x_sdd,params_sigmoid10[7][0],params_sigmoid10[7][1],params_sigmoid10[7][2],
                         params_sigmoid10[7][3])
fsca_sigmoid9 = fit_func(x_sdd,params_sigmoid10[8][0],params_sigmoid10[8][1],params_sigmoid10[8][2],
                         params_sigmoid10[8][3])
                        
fsca_sigmoid10 = fit_func(x_sdd,params_sigmoid102[0],params_sigmoid102[1],params_sigmoid102[2],
                          params_sigmoid102[3])
fsca_sigmoid12 = fit_func(x_sdd,params_sigmoid12[0],params_sigmoid12[1],params_sigmoid12[2],
                          params_sigmoid12[3])
fsca_sigmoid13 = fit_func(x_sdd,params_sigmoid13[0],params_sigmoid13[1],params_sigmoid13[2],
                          params_sigmoid13[3])

fsca_sigmoid15 = fit_func(x_sdd,params_sigmoid15_27[0][0],params_sigmoid15_27[0][1],params_sigmoid15_27[0][2],
                          params_sigmoid15_27[0][3])
fsca_sigmoid16 = fit_func(x_sdd,params_sigmoid15_27[1][0],params_sigmoid15_27[1][1],params_sigmoid15_27[1][2],
                          params_sigmoid15_27[1][3])
fsca_sigmoid17 = fit_func(x_sdd,params_sigmoid15_27[2][0],params_sigmoid15_27[2][1],params_sigmoid15_27[2][2],
                          params_sigmoid15_27[2][3])
fsca_sigmoid18 = fit_func(x_sdd,params_sigmoid15_27[3][0],params_sigmoid15_27[3][1],params_sigmoid15_27[3][2],
                          params_sigmoid15_27[3][3])
fsca_sigmoid19 = fit_func(x_sdd,params_sigmoid15_27[4][0],params_sigmoid15_27[4][1],params_sigmoid15_27[4][2],
                          params_sigmoid15_27[4][3])
fsca_sigmoid20 = fit_func(x_sdd,params_sigmoid15_27[5][0],params_sigmoid15_27[5][1],params_sigmoid15_27[5][2],
                          params_sigmoid15_27[5][3])
fsca_sigmoid21 = fit_func(x_sdd,params_sigmoid15_27[6][0],params_sigmoid15_27[6][1],params_sigmoid15_27[6][2],
                          params_sigmoid15_27[6][3])
fsca_sigmoid22 = fit_func(x_sdd,params_sigmoid15_27[7][0],params_sigmoid15_27[7][1],params_sigmoid15_27[7][2],
                          params_sigmoid15_27[7][3])
fsca_sigmoid23 = fit_func(x_sdd,params_sigmoid15_27[8][0],params_sigmoid15_27[8][1],params_sigmoid15_27[8][2],
                          params_sigmoid15_27[8][3])
fsca_sigmoid24 = fit_func(x_sdd,params_sigmoid15_27[9][0],params_sigmoid15_27[9][1],params_sigmoid15_27[9][2],
                          params_sigmoid15_27[9][3])
fsca_sigmoid25 = fit_func(x_sdd,params_sigmoid15_27[10][0],params_sigmoid15_27[10][1],params_sigmoid15_27[10][2],
                          params_sigmoid15_27[10][3])
fsca_sigmoid26 = fit_func(x_sdd,params_sigmoid15_27[11][0],params_sigmoid15_27[11][1],params_sigmoid15_27[11][2],
                          params_sigmoid15_27[11][3])
fsca_sigmoid27 = fit_func(x_sdd,params_sigmoid15_27[12][0],params_sigmoid15_27[12][1],params_sigmoid15_27[12][2],
                          params_sigmoid15_27[12][3])








 
titles = ['26Mar2016','17Apr2016','18May2016','p']
ylabel = ['fSCA','fSCA','fSCA','slope']
xlabel = ['SDD (days of year)','SDD (days of year)','SDD (days of year)','Daye']
sdd_fsca_fig = [sdd_fsca[:,17],sdd_fsca[:,39],sdd_fsca[:,69],sdd_fsca[:,17]]
parA = [10,1,1,1]
parB = [10,1,1,1]
parC = [0.6,1,1,1]
parD = [1,1,0.7,1]
x_fig = [[70, 82.,90,100,110,120,130],[82.,90,100,110,120,130,135],[70, 82.,90,100,110,120,130],[82.,90,100,110,120,130,135]]
fig,axes = plt.subplots(figsize=(50,35)) #,sharey='row', sharex=True, squeeze=True
for img in range (4):
    plt.subplot(221+img)
    plt.scatter(sdd_obs, sdd_fsca_fig[img], color='purple', s = 20**2, label='Data')
    fsca_sigmoid = fit_func(np.array(x_fig[img]),params39[0]*parA[img],
                            params39[1]*parB[img],params39[2]*parC[img],params39[3]*parD[img])
    plt.plot(np.array(x_fig[img]), fsca_sigmoid,linewidth= 8,color = 'black', label='Fitted function')
    
    legend = ['sigmoid curve fit','observation DTS data'] #.format(properties_rfm_pred_2080[sites][-1])
    plt.legend(legend,fontsize = 35)#, loc = 'upper center'

    plt.xticks(fontsize=45) #rotation='vertical', 
    plt.yticks(fontsize=40)
    plt.title(titles[img], fontsize=55)
    plt.ylabel(ylabel[img], fontsize=50)
    plt.xlabel(xlabel[img], fontsize=50)

plt.savefig('C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/TihoDTS/fSCA_sdd_sigmoid1.png')















