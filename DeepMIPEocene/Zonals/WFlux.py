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
LandFrac=Land
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

path1a='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2PrecDeepMIP.nc'
d1a=Dataset(path1a,'r')
Pr1=d1a.variables['pr'][:,:,:]
Prec1=np.mean(Pr1,axis=0)
Prec1[np.isnan(Land)]=np.nan

#3x CO2
path3='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2tsDeepMIP.nc'
d3=Dataset(path3,'r')
TS3=d3.variables['ts'][:,:,:]
TS3=np.mean(TS3,axis=0)
TS3[np.isnan(Land)]=np.nan

path3a='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2PrecDeepMIP.nc'
d3a=Dataset(path3a,'r')
Pr3=d3a.variables['pr'][:,:,:]
Prec3=np.mean(Pr3,axis=0)
Prec3[np.isnan(Land)]=np.nan

#6x CO2
path6='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2tsDeepMIP.nc'
d6=Dataset(path6,'r')
TS6=d6.variables['ts'][:,:,:]
TS6=np.mean(TS6,axis=0)
TS6[np.isnan(Land)]=np.nan

path6a='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2PrecDeepMIP.nc'
d6a=Dataset(path6a,'r')
Pr6=d6a.variables['pr'][:,:,:]
Prec6=np.mean(Pr6,axis=0)
Prec6[np.isnan(Land)]=np.nan


#9x CO2
path9='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2tsDeepMIP.nc'
d9=Dataset(path9,'r')
TS9=d9.variables['ts'][:,:,:]
TS9=np.mean(TS9,axis=0)
TS9[np.isnan(Land)]=np.nan

path9a='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2PrecDeepMIP.nc'
d9a=Dataset(path9a,'r')
Pr9=d9a.variables['pr'][:,:,:]
Prec9=np.mean(Pr9,axis=0)
Prec9[np.isnan(Land)]=np.nan

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
co23=840 #ppm
co26=1680 #ppm
co29=2520 #ppm

#Equations
#1x CO2
Prec1=(Prec1/rho)*60*60*24*365.25*1000
R1=0.018*Prec1**1.441

Pterm1=(co21/C0)**0.5
Tterm1=((np.exp((((-Eah*1000)/Rc)*(1/TS1-1/T0)))))
Rterm1=R1*bsilh
F1=Pterm1*Tterm1*Rterm1
F1[F1<0]=0
Fg1=F1*LandFrac


#3x CO2
Prec3=(Prec3/rho)*60*60*24*365.25*1000
R3=0.018*Prec3**1.441

Pterm3=(co23/C0)**0.5
Tterm3=((np.exp((((-Eah*1000)/Rc)*(1/TS3-1/T0)))))
Rterm3=R3*bsilh
F3=Pterm3*Tterm3*Rterm3
F3[F3<0]=0
Fg3=F3*LandFrac

#6x CO2
Prec6=(Prec6/rho)*60*60*24*365.25*1000
R6=0.018*Prec6**1.441

Pterm6=(co26/C0)**0.5
Tterm6=((np.exp((((-Eah*1000)/Rc)*(1/TS6-1/T0)))))
Rterm6=R6*bsilh
F6=Pterm6*Tterm6*Rterm6
F6[F6<0]=0
Fg6=F6*LandFrac

#9x CO2
Prec9=(Prec9/rho)*60*60*24*365.25*1000
R9=0.018*Prec9**1.441

Pterm9=(co29/C0)**0.5
Tterm9=((np.exp((((-Eah*1000)/Rc)*(1/TS9-1/T0)))))
Rterm9=R9*bsilh
F9=Pterm9*Tterm9*Rterm9
F9[F9<0]=0
Fg9=F9*LandFrac

#Calculate zonal contribution to total weathering
#Conversion to mol C
Mol1=Fg1*A/44.011 #mol C
Mol3=Fg3*A/44.011 # mol C
Mol6=Fg6*A/44.011 # mol C
Mol9=Fg9*A/44.011 # mol C

zonal1=np.nansum(Mol1,axis=1)
zonal3=np.nansum(Mol3,axis=1)
zonal6=np.nansum(Mol6,axis=1)
zonal9=np.nansum(Mol9,axis=1)


#Plotting the Data
plt.plot(lat,zonal1,label='1x $CO_2$')
plt.plot(lat,zonal3,label='3x $CO_2$')
plt.plot(lat,zonal6,label='6x $CO_2$')
plt.plot(lat,zonal9,label='9x $CO_2$')

plt.legend(title='Case', loc='upper right',fontsize=14)
plt.ylabel('Weathering Flux (mol yr$^{-1}$)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
plt.yscale('log')
#plt.title('Zonal Contribution to Total Weathering for DeepMIP Model Runs',fontsize=16)
plt.show()

