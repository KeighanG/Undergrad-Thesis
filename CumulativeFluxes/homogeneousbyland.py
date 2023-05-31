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

A_weighted_global=A/np.sum(A)

A=A*Land #Gives area of Land in each grid box
LandTotal=np.nansum(A) #Gives total area of land
A_weighted=A/LandTotal #Gives weighting of land in each grid box
LandFlat=np.ravel(A) #Flattened land array

#-----------------------Load in Weathering Arrays------------------------------#
Control=np.load('/home/kgemmell/weathering/CESM/CESMcontrol/ControlHomogeneousW.npy')
two=np.load('/home/kgemmell/weathering/CESM/2xCO2/2xHomogeneousW.npy')
four=np.load('/home/kgemmell/weathering/CESM/4xCO2/4xHomogeneousW.npy')


piControlCont=np.load('CMIPcont1.npy')
piControlLand=np.load('CMIPLand1.npy')
#----------------------Make Cumulative Flux Curves-----------------------------#
#Control
Flat=np.ravel(Control*A_weighted)
SortedF=np.argsort(-Flat)
Fluxes=Flat[SortedF]

As=LandFlat[SortedF]

Total=np.nansum(Control*A_weighted)

Cont=((np.nancumsum(Fluxes))/Total)*100
Lands=((np.nancumsum(As))/LandTotal)*100

#2x
Flat2=np.ravel(two*A_weighted)
SortedF2=np.argsort(-Flat2)
Fluxes2=Flat2[SortedF2]

As2=LandFlat[SortedF2]

Total2=np.nansum(two*A_weighted)

Cont2=((np.nancumsum(Fluxes2))/Total2)*100
Lands2=((np.nancumsum(As2))/LandTotal)*100

#4x
Flat4=np.ravel(four*A_weighted)
SortedF4=np.argsort(-Flat4)
Fluxes4=Flat4[SortedF4]

As4=LandFlat[SortedF4]

Total4=np.nansum(four*A_weighted)

Cont4=((np.nancumsum(Fluxes4))/Total4)*100
Lands4=((np.nancumsum(As4))/LandTotal)*100

#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(piControlLand,piControlCont,'--k',label='piControl Lithology')
plt.plot(Lands,Cont,'-k',label='piControl')
plt.plot(Lands2,Cont2,'-r',label='2x $CO_2$')
plt.plot(Lands4,Cont4,'-b',label='4x $CO_2$')

plt.ylabel('Accumulated Weathering Flux (%)',fontsize=14)
plt.xlabel('Accumulated Land Area (%)',fontsize=14)
#plt.title('Cumulative Flux-Area Rating Curves for CMIP-6 Experiments',fontsize=20)
plt.legend(title='Case', loc='lower right',fontsize=14)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.show()

