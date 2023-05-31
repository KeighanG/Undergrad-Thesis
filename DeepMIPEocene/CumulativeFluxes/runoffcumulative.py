from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import csv
from numpy import ndarray
#------------------------------------------------------------------------------#
#Read in Data
path='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
Landmask=d.variables['sftlf'][:,:,:]
Land=np.mean(Landmask,axis=0)
Land=Land.data
Landfrac=Land
Land[Land<=0]=np.nan
#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(Land.shape)
for i in range(0,len(lat)-1):
  for j in range(0,len(lon)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat[i]*np.pi/180.)-np.sin(lat[i+1]*np.pi/180.))*np.abs(lon[j]-lon[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

A=A*Land #Gives area of Land in each grid box
LandTotal=np.nansum(A) #Gives total area of land
A_weighted=A/LandTotal #Gives weighting of land in each grid box
LandFlat=np.ravel(A) #Flattened land array

#-----------------------Load in Weathering Arrays------------------------------#
#1x CO2
path1='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2tsDeepMIP.nc'
d1=Dataset(path1,'r')
TS1=d1.variables['ts'][:,:,:]
TS1=np.mean(TS1,axis=0)
TS1[np.isnan(Land)]=np.nan
TS1=TS1.data

path1a='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2PrecDeepMIP.nc'
d1a=Dataset(path1a,'r')
Pr1=d1a.variables['pr'][:,:,:]
Prec1=np.mean(Pr1,axis=0)
Prec1[np.isnan(Land)]=np.nan
Prec1=Prec1.data

#3x CO2
path3='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2tsDeepMIP.nc'
d3=Dataset(path3,'r')
TS3=d3.variables['ts'][:,:,:]
TS3=np.mean(TS3,axis=0)
TS3[np.isnan(Land)]=np.nan
TS3=TS3.data

path3a='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2PrecDeepMIP.nc'
d3a=Dataset(path3a,'r')
Pr3=d3a.variables['pr'][:,:,:]
Prec3=np.mean(Pr3,axis=0)
Prec3[np.isnan(Land)]=np.nan
Prec3=Prec3.data

#6x CO2
path6='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2tsDeepMIP.nc'
d6=Dataset(path6,'r')
TS6=d6.variables['ts'][:,:,:]
TS6=np.mean(TS6,axis=0)
TS6[np.isnan(Land)]=np.nan
TS6=TS6.data

path6a='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2PrecDeepMIP.nc'
d6a=Dataset(path6a,'r')
Pr6=d6a.variables['pr'][:,:,:]
Prec6=np.mean(Pr6,axis=0)
Prec6[np.isnan(Land)]=np.nan
Prec6=Prec6.data

#9x CO2
path9='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2tsDeepMIP.nc'
d9=Dataset(path9,'r')
TS9=d1.variables['ts'][:,:,:]
TS9=np.mean(TS9,axis=0)
TS9[np.isnan(Land)]=np.nan
TS9=TS9.data

path9a='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2PrecDeepMIP.nc'
d9a=Dataset(path9a,'r')
Pr9=d9a.variables['pr'][:,:,:]
Prec9=np.mean(Pr9,axis=0)
Prec9[np.isnan(Land)]=np.nan
Prec9=Prec9.data

#Define Constants
rho=1000 #kg/m3, density of freshwater
Rc=8.314 #J/mol/K
B=0.65 #From Colbourn
Lvap=2.26*10**6 #J/kg
C0=280 #ppm
T0=287.1968543493361 #K, from CESM piControl
Eah=55 #Homogenous Ea, KJ/mol
bsilh=0.02 #Homogeneous bsil
co21=280 #ppm
co23=840
co26=1680
co29=2520 #ppm

#Equations
#1x CO2
R1=0.018*Prec1**1.441
F1=R1

#3x CO2
R3=0.018*Prec3**1.441
F3=R3

#6x CO2
R6=0.018*Prec6**1.441
F6=R6

#9x CO2
R9=0.018*Prec9**1.441
F9=R9
#----------------------Make Cumulative Flux Curves-----------------------------#
#1x
Flat1=np.ravel(F1*A_weighted)
SortedF1=np.argsort(-Flat1)
Fluxes1=Flat1[SortedF1]

As1=LandFlat[SortedF1]

Total1=np.nansum(F1*A_weighted)

Cont1=((np.nancumsum(Fluxes1))/Total1)*100
Lands1=((np.nancumsum(As1))/LandTotal)*100

#3x
Flat3=np.ravel(F3*A_weighted)
SortedF3=np.argsort(-Flat3)
Fluxes3=Flat3[SortedF3]

As3=LandFlat[SortedF3]

Total3=np.nansum(F3*A_weighted)

Cont3=((np.nancumsum(Fluxes3))/Total3)*100
Lands3=((np.nancumsum(As3))/LandTotal)*100

#6x
Flat6=np.ravel(F6*A_weighted)
SortedF6=np.argsort(-Flat6)
Fluxes6=Flat6[SortedF6]

As6=LandFlat[SortedF6]

Total6=np.nansum(F6*A_weighted)

Cont6=((np.nancumsum(Fluxes6))/Total6)*100
Lands6=((np.nancumsum(As6))/LandTotal)*100

#9x
Flat9=np.ravel(F9*A_weighted)
SortedF9=np.argsort(-Flat9)
Fluxes9=Flat9[SortedF9]

As9=LandFlat[SortedF9]

Total9=np.nansum(F9*A_weighted)

Cont9=((np.nancumsum(Fluxes9))/Total9)*100
Lands9=((np.nancumsum(As9))/LandTotal)*100

#------------------------------------------------------------------------------#
#Control
Lands=np.load('ControlLands.npy')
Cont=np.load('ControlPrecFlux.npy')

#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(Lands,Cont,'k',label='piControl')
plt.plot(Lands1,Cont1,label='1x $CO_2$')
plt.plot(Lands3,Cont3,label='3x $CO_2$')
plt.plot(Lands6,Cont6,label='6x $CO_2$')
plt.plot(Lands9,Cont9,label='9x $CO_2$')

plt.ylabel('Cumulative Runoff (%)',fontsize=14)
plt.xlabel('Accumulated Land Area (%)',fontsize=14)
#plt.title('Cumulative Flux-Area Rating Curves for DeepMIP Eocene Data',fontsize=20)
plt.legend(title='Case', loc='upper left',fontsize=14)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.show()

