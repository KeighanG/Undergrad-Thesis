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
Land[Land<1]=np.nan

path1='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2PrecDeepMIP.nc'
d1=Dataset(path1,'r')
P1=d1.variables['pr'][:,:,:]
P1=np.mean(P1,axis=0)
P1[np.isnan(Land)]=np.nan
zonalP1=np.nanmean(P1,axis=1)

path3='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2PrecDeepMIP.nc'
d3=Dataset(path3,'r')
P3=d3.variables['pr'][:,:,:]
P3=np.mean(P3,axis=0)
P3[np.isnan(Land)]=np.nan
zonalP3=np.nanmean(P3,axis=1)

path6='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2PrecDeepMIP.nc'
d6=Dataset(path6,'r')
P6=d6.variables['pr'][:,:,:]
P6=np.mean(P6,axis=0)
P6[np.isnan(Land)]=np.nan
zonalP6=np.nanmean(P6,axis=1)

path9='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2PrecDeepMIP.nc'
d9=Dataset(path9,'r')
P9=d9.variables['pr'][:,:,:]
P9=np.mean(P9,axis=0)
P9[np.isnan(Land)]=np.nan
zonalP9=np.nanmean(P9,axis=1)

zonalP1=zonalP1*60*60*24*365.25
zonalP3=zonalP3*60*60*24*365.25
zonalP6=zonalP6*60*60*24*365.25
zonalP9=zonalP9*60*60*24*365.25

#Plotting the Data
plt.plot(lat,zonalP1,'k',label='1x $CO_2$')
plt.plot(lat,zonalP3,label='3x $CO_2$')
plt.plot(lat,zonalP6,label='6x $CO_2$')
plt.plot(lat,zonalP9,label='9x $CO_2$')

plt.legend(title='Atmospheric $CO_2$ Concentration', loc='upper right')

plt.ylabel('Precipitation (mm/yr)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
plt.title('Zonally Averaged Precipitation Over Land for DeepMIP Eocene Data',fontsize=16)

plt.show()

