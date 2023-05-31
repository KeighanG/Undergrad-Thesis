from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
lat=np.linspace(90,-90,360)
lon=np.linspace(0,360,720)

eight=np.load('Flux8.npy')
nine=np.load('Flux9.npy')
one=np.load('Flux10.npy')
eleven=np.load('Flux11.npy')
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

tot8=eight*A/44.011
tot9=nine*A/44.011
tot1=one*A/44.011
tot11=eleven*A/44.011

#Multiply by area, and then sum. Gives mol/yr/latitude
zonal8=np.nansum(tot8,axis=1)
zonal9=np.nansum(tot9,axis=1)
zonal1=np.nansum(tot1,axis=1)
zonal11=np.nansum(tot11,axis=1)

zonal8[zonal8<=0]=np.nan
zonal9[zonal9<=0]=np.nan
zonal1[zonal1<=0]=np.nan
zonal11[zonal11<=0]=np.nan
#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(lat,zonal8,label='$S/S_0=0.8$')
plt.plot(lat,zonal9,label='$S/S_0=0.9$')
plt.plot(lat,zonal1,label='$S/S_0=1.0$')
plt.plot(lat,zonal11,label='$S/S_0=1.1$')


plt.legend(title='Solar Constant Ratio to Present Day', loc='lower center',fontsize=14)

plt.ylabel('Weathering Flux (mol yr$^{-1}$)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
plt.yscale('log')
#plt.title('Zonal Contribution to Total Weathering for CESM Model Runs: Lithology',fontsize=16)
plt.show()

