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
Control=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/LithW.npy')
two=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/2xLithW.npy')
four=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/4xLithW.npy')

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


CMIP115=np.sum(np.flip(Control[150:180])+Control[180:210])/np.sum(Control)*100
CMIP130=np.sum(np.flip(Control[120:150])+Control[210:240])/np.sum(Control)*100
CMIP160=np.sum(np.flip(Control[0:120])+Control[240:360])/np.sum(Control)*100
print(CMIP115)
print(CMIP130)
print(CMIP160)

CMIP215=np.sum(np.flip(two[150:180])+two[180:210])/np.sum(two)*100
CMIP230=np.sum(np.flip(two[120:150])+two[210:240])/np.sum(two)*100
CMIP260=np.sum(np.flip(two[0:120])+two[240:360])/np.sum(two)*100
print(CMIP215)
print(CMIP230)
print(CMIP260)

CMIP415=np.sum(np.flip(four[150:180])+four[180:210])/np.sum(four)*100
CMIP430=np.sum(np.flip(four[120:150])+four[210:240])/np.sum(four)*100
CMIP460=np.sum(np.flip(four[0:120])+four[240:360])/np.sum(four)*100
print(CMIP415)
print(CMIP430)
print(CMIP460)
#These printed numbers are copied into the latitudinal bands script

