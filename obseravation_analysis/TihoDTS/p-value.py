import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy import stats

with open("C:/1UNRuniversityFolder/Dissertation/Chapter 2-snow-forest/obseravation_analysis/TihoDTS/fSCA-SDD.csv") as safd:
    reader = csv.reader(safd)
    data_forcing = [r for r in reader]
data_forcing2 = data_forcing[1:]
sa_fd_column = []
for csv_counter1 in range (len (data_forcing2)):
    for csv_counter2 in range (71):
        sa_fd_column.append(float(data_forcing2[csv_counter1][csv_counter2]))
sdd_fsca=np.reshape(sa_fd_column,(len (data_forcing2),71))

slope, intercept, r_value, p_value, std_err = stats.linregress(sdd_fsca[:,0], sdd_fsca[:,2])
