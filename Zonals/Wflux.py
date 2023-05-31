from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
lat=np.linspace(90,-90,360)
lon=np.linspace(0,360,720)

control=np.load('LithW.npy')
twox=np.load('2xLithW.npy')
fourx=np.load('4xLithW.npy')

zonalcontrol=np.nanmean(control,axis=1)
zonal2x=np.nanmean(twox,axis=1)
zonal4x=np.nanmean(fourx,axis=1)

#------------------------------------------------------------------------------#
#Define grid spacial resolutions, x1 indicates Hartmann grid
X, Y =np.meshgrid(lon,lat)
lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)
XI,YI=np.meshgrid(lon1,lat1)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Ea= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bsil= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bcarb= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Removing Oceans from Hartmann data
Land[Land<0]=np.nan
Land[Land>0]=1

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(ascii_grid.shape)
for i in range(0,len(lat1)-1):
  for j in range(0,len(lon1)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat1[i]*np.pi/180.)-np.sin(lat1[i+1]*np.pi/180.))*np.abs(lon1[j]-lon1[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

A_weighted_global=A/np.sum(A)

A=A*Land #Gives area of Land in each grid box
A2=np.nansum(A) #Gives total area of land
A_weighted=A/A2 #Gives weighting of land in each grid box


ZonalLand=np.nansum(Land,axis=1)

zonalcontrol=zonalcontrol*ZonalLand
zonal2x=zonal2x*ZonalLand
zonal4x=zonal4x*ZonalLand
#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(lat,zonalcontrol,'k',label='PI Control')
plt.plot(lat,zonal2x,label='2x $CO_2$')
plt.plot(lat,zonal4x,label='4x $CO_2$')

plt.legend(title='Model Run', loc='upper right')

plt.ylabel('Weathering Flux (mol km$^{-1}$ a$^{-1}$)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)

plt.show()
plt.savefig('ZonalWFlux.png')

