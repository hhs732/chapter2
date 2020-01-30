import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy import stats
from scipy import optimize

with open("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/TihoDTS/fSCA-SDD.csv") as safd:
    reader = csv.reader(safd)
    data_forcing = [r for r in reader]
data_forcing2 = data_forcing[1:]
sa_fd_column = []
for csv_counter1 in range (len (data_forcing2)):
    for csv_counter2 in range (71):
        sa_fd_column.append(float(data_forcing2[csv_counter1][csv_counter2]))
sdd_fsca=np.reshape(sa_fd_column,(len (data_forcing2),71))


#slope, intercept, r_value, p_value, std_err = stats.linregress(sdd_fsca[:,0], sdd_fsca[:,2])

def fit_func(x, a, b, c, d):
    return a/(b+c*np.exp(-d*x))

sdd_obs = sdd_fsca[:,0]

params39, params_covariance39 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,39],
                                                   p0=[0.03,-20.21,10.21,-0.01])
params17, params_covariance17 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,17],
                                                   p0=[1,1,0.5,0.5])
params65, params_covariance65 = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,52],
                                                   p0=[0.5,0.5,0.5,0.5])
params_sigmoid = []
covariance_sigmoid =[]
for fsca in range (33):#len(sdd_fsca[0,:])-1
    params, params_covariance = optimize.curve_fit(fit_func, sdd_obs, sdd_fsca[:,1+fsca],
                                                   p0=[0.03,-20.21,10.21,-0.01])
    params_sigmoid.append(params)
    covariance_sigmoid.append(params_covariance)
    
    
plt.figure(figsize=(15, 10))
plt.scatter(sdd_obs, sdd_fsca[:,39], color='purple', s = 20**2, label='Data')
fsca_sigmoid = fit_func(np.array([70,82.,90,100,110,120,130]),params39[0]*1,
                        params39[1]*1,params39[2]*1,params39[3]*1)
plt.plot(np.array([70, 82.,90,100,110,120,130]), fsca_sigmoid,linewidth= 8,color = 'black', label='Fitted function')

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















