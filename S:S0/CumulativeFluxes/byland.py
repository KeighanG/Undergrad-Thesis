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
#Define grid spacial resolutions, x1 indicates Hartmann grid
lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

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
LandTotal=np.nansum(A) #Gives total area of land
A_weighted=A/LandTotal #Gives weighting of land in each grid box
LandFlat=np.ravel(A) #Flattened land array
#-----------------------Load in Weathering Arrays------------------------------#
eight=np.load('Flux8.npy')
nine=np.load('Flux9.npy')
one=np.load('Flux10.npy')
eleven=np.load('Flux11.npy')
#----------------------Make Cumulative Flux Curves-----------------------------#
#S/S0=0.8
Flat8=np.ravel(eight*A_weighted)
SortedF8=np.argsort(-Flat8)
Fluxes8=Flat8[SortedF8]

As8=LandFlat[SortedF8]

Total8=np.nansum(eight*A_weighted)

Cont8=((np.nancumsum(Fluxes8))/Total8)*100
Lands8=((np.nancumsum(As8))/LandTotal)*100

#S/S0=0.9
Flat9=np.ravel(nine*A_weighted)
SortedF9=np.argsort(-Flat9)
Fluxes9=Flat9[SortedF9]

As9=LandFlat[SortedF9]

Total9=np.nansum(nine*A_weighted)

Cont9=((np.nancumsum(Fluxes9))/Total9)*100
Lands9=((np.nancumsum(As9))/LandTotal)*100

#S/S0=1.0
Flat1=np.ravel(one*A_weighted)
SortedF1=np.argsort(-Flat1)
Fluxes1=Flat1[SortedF1]

As1=LandFlat[SortedF1]

Total1=np.nansum(one*A_weighted)

Cont1=((np.nancumsum(Fluxes1))/Total1)*100
Lands1=((np.nancumsum(As1))/LandTotal)*100

#S/S0=1.1
Flat11=np.ravel(eleven*A_weighted)
SortedF11=np.argsort(-Flat11)
Fluxes11=Flat11[SortedF11]

As11=LandFlat[SortedF11]

Total11=np.nansum(eleven*A_weighted)

Cont11=((np.nancumsum(Fluxes11))/Total11)*100
Lands11=((np.nancumsum(As11))/LandTotal)*100

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

A_weighted_globalc=Ac/np.sum(Ac)

Ac=Ac*Landc #Gives area of Land in each grid box
LandTotalc=np.nansum(Ac) #Gives total area of land
A_weightedc=Ac/LandTotalc #Gives weighting of land in each grid box
LandFlatc=np.ravel(Ac) #Flattened land array

Control=np.load('LithW.npy')

#Control
Flat=np.ravel(Control*A_weightedc)
SortedF=np.argsort(-Flat)
Fluxes=Flat[SortedF]

As=LandFlatc[SortedF]

Total=np.nansum(Control*A_weightedc)

Cont=((np.nancumsum(Fluxes))/Total)*100
Lands=((np.nancumsum(As))/LandTotalc)*100

#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(Lands,Cont,'k',label='piControl')
plt.plot(Lands8,Cont8,label='$S/S_0=0.8$')
plt.plot(Lands9,Cont9,label='$S/S_0=0.9$')
plt.plot(Lands1,Cont1,label='$S/S_0=1.0$')
plt.plot(Lands11,Cont11,'-b',label='$S/S_0=1.1$')

plt.ylabel('Accumulated Weathering Flux (%)',fontsize=14)
plt.xlabel('Accumulated Land Area (%)',fontsize=14)
#plt.title('Cumulative Flux-Area Rating Curves for CAM4 Model Runs',fontsize=20)
plt.legend(title='Case', loc='upper left',fontsize=14)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.show()


