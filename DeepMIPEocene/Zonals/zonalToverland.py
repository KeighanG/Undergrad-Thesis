from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
#Read in Data
path='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
Landmask=d.variables['sftlf'][:,:,:]
Land=np.mean(Landmask,axis=0)
Land[Land<=0]=np.nan

path1='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2tsDeepMIP.nc'
d1=Dataset(path1,'r')
TS1=d1.variables['ts'][:,:,:]
TS1=np.mean(TS1,axis=0)
TS1[np.isnan(Land)]=np.nan
zonalT1=np.nanmean(TS1,axis=1)

path3='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2tsDeepMIP.nc'
d3=Dataset(path3,'r')
TS3=d3.variables['ts'][:,:,:]
TS3=np.mean(TS3,axis=0)
TS3[np.isnan(Land)]=np.nan
zonalT3=np.nanmean(TS3,axis=1)

path6='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2tsDeepMIP.nc'
d6=Dataset(path6,'r')
TS6=d6.variables['ts'][:,:,:]
TS6=np.mean(TS6,axis=0)
TS6[np.isnan(Land)]=np.nan
zonalT6=np.nanmean(TS6,axis=1)

path9='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2tsDeepMIP.nc'
d9=Dataset(path9,'r')
TS9=d9.variables['ts'][:,:,:]
TS9=np.mean(TS9,axis=0)
TS9[np.isnan(Land)]=np.nan
zonalT9=np.nanmean(TS9,axis=1)

#Plotting the Data
plt.plot(lat,zonalT1,label='1x $CO_2$')
plt.plot(lat,zonalT3,label='3x $CO_2$')
plt.plot(lat,zonalT6,label='6x $CO_2$')
plt.plot(lat,zonalT9,label='9x $CO_2$')

plt.legend(title='Case', loc='upper right',fontsize=14)

plt.ylabel('Temperature (K)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
#plt.title('Zonally Averaged Temperature Over Land for DeepMIP Eocene Data',fontsize=16)

plt.show()

