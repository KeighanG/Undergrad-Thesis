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
#Read in LandMask Data for Eocene
path='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]

#Land fraction
LANDFRAC=d.variables['sftlf'][:,:,:]
Land=np.mean(LANDFRAC,axis=0)
Landmask=Land
Land[Land<=0]=np.nan

path1='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d1=Dataset(path1,'r')
lat1=d.variables['lat'][:]
lon1=d.variables['lon'][:]

#Land fraction
LANDFRAC1=d1.variables['LANDFRAC'][:,:,:][240:,:,:]
Land1=np.mean(LANDFRAC1,axis=0)
Landmask1=Land1
Land1[Land1<=0]=np.nan

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(Land.shape)
for i in range(0,len(lat)-1):
  for j in range(0,len(lon)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat[i]*np.pi/180.)-np.sin(lat[i+1]*np.pi/180.))*np.abs(lon[j]-lon[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

Aeocene=Land*A*Landmask
Apresent=Land1*A*Landmask1

zonalLandEocene=np.nansum(Aeocene,axis=1)
zonalLandToday=np.nansum(Apresent,axis=1)

total1=np.nansum(zonalLandEocene)
total2=np.nansum(zonalLandToday)
#------------------------------------------------------------------------------#
#Plotting the Data
plt.plot(lat,zonalLandEocene,label='Eocene')
plt.plot(lat1,zonalLandToday,label='Present Day')

plt.legend(title='Continetal Configuration', loc='upper left',fontsize=14)
#plt.title('Zonal Land Distribution for Present Day and Eocene Continental Configurations',fontsize=20)
plt.ylabel('Land Area ($m^2$)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)

plt.show()

print(total1)
print(total2)
print(total1/total2)

