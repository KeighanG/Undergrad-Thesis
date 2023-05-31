import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import csv
from numpy import ndarray
#------------------------------------------------------------------------------#
#Define grid spacial resolutions, x1 indicates Hartmann grid
lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)

#Load in Hartmann model data
ascii_grid = np.loadtxt("/home/kgemmell/weathering/CESM/CumFluxLand/glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("/home/kgemmell/weathering/CESM/CumFluxLand/glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Setting Ocean tiles to 0 which in the sum over all will lead to no contribution to the weathering flux
Land[Land<0]=0
Land[Land>0]=1

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(ascii_grid.shape)
for i in range(0,len(lat1)-1):
  for j in range(0,len(lon1)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat1[i]*np.pi/180.)-np.sin(lat1[i+1]*np.pi/180.))*np.abs(lon1[j]-lon1[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

#-----------------------Load in Weathering Arrays------------------------------#
Control=np.load('/home/kgemmell/weathering/CESM/CESMcontrol/ControlHomogeneousW.npy')
two=np.load('/home/kgemmell/weathering/CESM/2xCO2/2xHomogeneousW.npy')
four=np.load('/home/kgemmell/weathering/CESM/4xCO2/4xHomogeneousW.npy')

#Set NaN to 0 such that the sum calculations are not affected
Control[np.isnan(Control)]=0
#Calculate the total latitdudinal weathering flux by multiplying the gridbox average weathering
#flux with the land area in each grid box and then make a zonal mean
Control=np.sum(Control*A,axis=1)

#Same as above for two times CO_2
two[np.isnan(two)]=0
two=np.sum(two*A,axis=1)

#Same as above for 4 times CO_2
four[np.isnan(four)]=0
four=np.sum(four*A,axis=1)

#Load in CMIP-6 Data
lat=np.linspace(90,-90,360)
CMIP1=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIP1lat.npy')

#Plot the cumulative weathering flux by flipping the first half of the total zonal weathering flux
#such that the array goes from 0 to 90N and add the hemisphere from 0 to 90 S. Then take the cumulative
#sum and devide by the total global averaging flux to normalize
plt.plot(np.flip(lat[0:180]),CMIP1,'--k',label='piControl Lithology')

plt.plot(np.flip(lat1[0:180]),np.cumsum(np.flip(Control[0:180])+Control[180:])/np.sum(Control)*100,'k',label='piControl')
#Same as above for two times CO_2
plt.plot(np.flip(lat1[0:180]),np.cumsum(np.flip(two[0:180])+two[180:])/np.sum(two)*100,'r',label='2x $CO_2$')
#Same as above for four times CO_2
plt.plot(np.flip(lat1[0:180]),np.cumsum(np.flip(four[0:180])+four[180:])/np.sum(four)*100,'b',label='4x $CO_2$')
plt.legend(title='Case', loc='lower right',fontsize=14)

plt.xlabel('Equator to Pole Latitude (+ or - up to 90)',fontsize=14)
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.ylabel('Cumulative Flux (%)',fontsize=14)
#plt.title('Cumulative Flux by Latitude Starting at Equator for CMIP-6 Experiments',fontsize=20)
plt.show()

