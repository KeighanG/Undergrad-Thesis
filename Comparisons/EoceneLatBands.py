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
lat=lat.data
lat=np.flip(lat)
lon=d.variables['lon'][:]
lon=lon.data
Landmask=d.variables['sftlf'][:,:,:]
Land=np.mean(Landmask,axis=0)
Land=Land.data
LandFrac=Land
Land[Land<=0]=0

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(Land.shape)
for i in range(0,len(lat)-1):
  for j in range(0,len(lon)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat[i]*np.pi/180.)-np.sin(lat[i+1]*np.pi/180.))*np.abs(lon[j]-lon[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

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

#----------------------Make Cumulative Flux Curves-----------------------------#
Fg1[np.isnan(Fg1)]=0
Fg1=np.sum(Fg1*A,axis=1)

Fg3[np.isnan(Fg3)]=0
Fg3=np.sum(Fg3*A,axis=1)

Fg6[np.isnan(Fg6)]=0
Fg6=np.sum(Fg6*A,axis=1)

Fg9[np.isnan(Fg9)]=0
Fg9=np.sum(Fg9*A,axis=1)

#CONTROL RUN
#------------------------------------------------------------------------------#
#Define grid spacial resolutions, x1 indicates Hartmann grid
lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Landc = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Removing Oceans from Hartmann data
Landc[Landc<0]=np.nan
Landc[Landc>0]=1

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
Ac=np.zeros(ascii_grid.shape)
for i in range(0,len(lat1)-1):
  for j in range(0,len(lon1)-1):
    Ac[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat1[i]*np.pi/180.)-np.sin(lat1[i+1]*np.pi/180.))*np.abs(lon1[j]-lon1[j+1])
Ac[-1,:]=Ac[0,:]
Ac[:,-1]=Ac[:,0]


Control=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/LithW.npy')

#Control
#Set NaN to 0 such that the sum calculations are not affected
Control[np.isnan(Control)]=0
#Calculate the total latitdudinal weathering flux by multiplying the gridbox average weathering
#flux with the land area in each grid box and then make a zonal mean
Control=np.sum(Control*Ac,axis=1)

eocene115=np.sum(np.flip(Fg1[40:48])+Fg1[48:56])/np.sum(Fg1)*100
eocene130=np.sum(np.flip(Fg1[32:40])+Fg1[56:64])/np.sum(Fg1)*100
eocene190=np.sum(np.flip(Fg1[0:32])+Fg1[64:96])/np.sum(Fg1)*100
print(eocene115)
print(eocene130)
print(eocene190)


eocene315=np.sum(np.flip(Fg3[40:48])+Fg3[48:56])/np.sum(Fg3)*100
eocene330=np.sum(np.flip(Fg3[32:40])+Fg3[56:64])/np.sum(Fg3)*100
eocene390=np.sum(np.flip(Fg3[0:32])+Fg3[64:96])/np.sum(Fg3)*100
print(eocene315)
print(eocene330)
print(eocene390)

eocene615=np.sum(np.flip(Fg6[40:48])+Fg6[48:56])/np.sum(Fg6)*100
eocene630=np.sum(np.flip(Fg6[32:40])+Fg6[56:64])/np.sum(Fg6)*100
eocene690=np.sum(np.flip(Fg3[0:32])+Fg6[64:96])/np.sum(Fg6)*100
print(eocene615)
print(eocene630)
print(eocene690)

eocene915=np.sum(np.flip(Fg9[40:48])+Fg9[48:56])/np.sum(Fg9)*100
eocene930=np.sum(np.flip(Fg9[32:40])+Fg9[56:64])/np.sum(Fg9)*100
eocene990=np.sum(np.flip(Fg9[0:32])+Fg9[64:96])/np.sum(Fg9)*100
print(eocene915)
print(eocene930)
print(eocene990)

