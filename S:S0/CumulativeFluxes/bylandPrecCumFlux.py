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
Control=np.load('CESMPrec.npy')
one=np.load('SS0Prec.npy')

#Make Cumulative Curves
#Control
Flat=np.ravel(Control*A_weighted)
SortedF=np.argsort(-Flat)
Fluxes=Flat[SortedF]

As=LandFlat[SortedF]

Total=np.nansum(Control*A_weighted)

Cont=((np.nancumsum(Fluxes))/Total)*100
Lands=((np.nancumsum(As))/LandTotal)*100


#S/S0=1.0
Flat1=np.ravel(one*A_weighted)
SortedF1=np.argsort(-Flat1)
Fluxes1=Flat1[SortedF1]

As1=LandFlat[SortedF1]

Total1=np.nansum(one*A_weighted)

Cont1=((np.nancumsum(Fluxes1))/Total1)*100
Lands1=((np.nancumsum(As1))/LandTotal)*100

#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(Lands,Cont,'k',label='piControl')
plt.plot(Lands1,Cont1,label='$S/S_0=1.0$')

plt.ylabel('Accumulated Precipitation (%)',fontsize=14)
plt.xlabel('Accumulated Land Area (%)',fontsize=14)
#plt.title('Cumulative Precipitation-Area Rating Curves Comparison',fontsize=20)
plt.legend(title='Case', loc='upper left',fontsize=14)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.show()

