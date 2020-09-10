from scipy.io import loadmat
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt

fileA = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\GroundwaterDataA.mat'
keyA = 'GroundwaterDataA'

gwts = loadmat(fileA)[keyA]     # import mat file using loadmat - and only save timeseries array for dictionary key of interest

dt = []
gwz = []

for rec in gwts:
    dt.append(rec[0])
    gwz.append(rec[1])

del(gwts)

ymin = max(0.0, min(gwz)*0.9)
ymax = max(gwz)*1.1

fig = plt.figure()
ax = fig.add_subplot(111,facecolor='whitesmoke')
ax.plot(dt,gwz,color='blue')

ax.set_ylim(ymax,ymin)
ax.set_xlabel('Day of year')
ax.set_ylabel('Depth to water table [units?]')
plt.title('%s' % keyA)

plt.show()

    
