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
lon1=np.linspace(-180,180,720)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Removing Oceans from Hartmann data
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
eight=np.load('Flux8.npy')
nine=np.load('Flux9.npy')
one=np.load('Flux10.npy')
eleven=np.load('Flux11.npy')
#----------------------Make Cumulative Flux Curves-----------------------------#
eight[np.isnan(eight)]=0
eight=np.sum(eight*A,axis=1)

nine[np.isnan(nine)]=0
nine=np.sum(nine*A,axis=1)

one[np.isnan(one)]=0
one=np.sum(one*A,axis=1)

eleven[np.isnan(eleven)]=0
eleven=np.sum(eleven*A,axis=1)

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

SS0815=np.sum(np.flip(eight[150:180])+eight[180:210])/np.sum(eight)*100
SS0830=np.sum(np.flip(eight[120:150])+eight[210:240])/np.sum(eight)*100
SS0890=np.sum(np.flip(eight[0:120])+eight[240:360])/np.sum(eight)*100
print(SS0815)
print(SS0830)
print(SS0890)

SS0915=np.sum(np.flip(nine[150:180])+nine[180:210])/np.sum(nine)*100
SS0930=np.sum(np.flip(nine[120:150])+nine[210:240])/np.sum(nine)*100
SS0990=np.sum(np.flip(nine[0:120])+nine[240:360])/np.sum(nine)*100
print(SS0915)
print(SS0930)
print(SS0990)


SS0115=np.sum(np.flip(one[150:180])+one[180:210])/np.sum(one)*100
SS0130=np.sum(np.flip(one[120:150])+one[210:240])/np.sum(one)*100
SS0190=np.sum(np.flip(one[0:120])+one[240:360])/np.sum(one)*100
print(SS0115)
print(SS0130)
print(SS0190)


SS01115=np.sum(np.flip(eleven[150:180])+eleven[180:210])/np.sum(eleven)*100
SS01130=np.sum(np.flip(eleven[120:150])+eleven[210:240])/np.sum(eleven)*100
SS01190=np.sum(np.flip(eleven[0:120])+eleven[240:360])/np.sum(eleven)*100
print(SS01115)
print(SS01130)
print(SS01190)

