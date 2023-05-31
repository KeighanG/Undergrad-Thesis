import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
#Read in LandMask Data
path='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
Landmask=d.variables['sftlf'][:,:,:]
Land=np.mean(Landmask,axis=0)
LandMask=Land
Land[Land<=0]=np.nan
Land=np.roll(Land,72,axis=1)

#-------------------------Size of Each Gridbox---------------------------------#
#Eocene Configuration
r=6371000 #m
A=np.zeros(Land.shape)
for i in range(0,len(lat)-1):
  for j in range(0,len(lon)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat[i]*np.pi/180.)-np.sin(lat[i+1]*np.pi/180.))*np.abs(lon[j]-lon[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

L=A
A=A*Land #Gives area of Land in each grid box
LandTotal=np.nansum(A) #Gives total area of land
A_weighted=A/LandTotal #Gives weighting of land in each grid box
LandFlat=np.ravel(A) #Flattened land array

zonalL=np.nansum(L,axis=1)


Ahalf1pre=A_weighted[0:180,:]
Ahalfb=A_weighted[180:360,:]
Ahalfa=np.flip(Ahalf1pre)
zonalLandEocene=np.nansum(A,axis=1) #Gives the area of land per latitude

zonalFracEocene=zonalLandEocene/zonalL


#------------------------------------------------------------------------------#
#Present Day Configuration
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

Lc=Ac
zonalLc=np.nansum(Lc,axis=1)


A_weighted_globalc=Ac/np.sum(Ac)

Ac=Ac*Landc #Gives area of Land in each grid box
landfrac=Ac/Lc


zonalLandToday=np.nansum(Ac,axis=1)

zonalFracToday=zonalLandToday/zonalLc



#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(lat,zonalFracEocene,label='Eocene')
plt.plot(lat1,zonalFracToday,label='Present Day')


plt.legend(title='Continetal Configuration', loc='upper center',fontsize=14)
plt.ylabel('Fraction of Zonal Land Area',fontsize=14)
plt.xlabel('Latitude',fontsize=14)

plt.show()

