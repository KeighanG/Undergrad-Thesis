from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
#Read in Data
path='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
#Land fraction
LANDFRAC=d.variables['LANDFRAC'][:,:,:][240:,:,:]
Land=np.mean(LANDFRAC,axis=0)
Land[Land<=0]=np.nan

path8='/home/czg/FYSP_clouds_archive/CAM4/cam4_08.nc'
d8=Dataset(path8,'r')
TS8=d8.variables['TS'][:,:,:]
TS8=np.mean(TS8,axis=0)
TS8[np.isnan(Land)]=np.nan
zonalT8=np.nanmean(TS8,axis=1)

path9='/home/czg/FYSP_clouds_archive/CAM4/cam4_09.nc'
d9=Dataset(path9,'r')
TS9=d9.variables['TS'][:,:,:]
TS9=np.mean(TS9,axis=0)
TS9[np.isnan(Land)]=np.nan
zonalT9=np.nanmean(TS9,axis=1)

path1='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d1=Dataset(path1,'r')
TS1=d1.variables['TS'][:,:,:]
TS1=np.mean(TS1,axis=0)
TS1[np.isnan(Land)]=np.nan
zonalT1=np.nanmean(TS1,axis=1)

path11='/home/czg/FYSP_clouds_archive/CAM4/cam4_11.nc'
d11=Dataset(path11,'r')
TS11=d11.variables['TS'][:,:,:]
TS11=np.mean(TS11,axis=0)
TS11[np.isnan(Land)]=np.nan
zonalT11=np.nanmean(TS11,axis=1)

#Plotting the Data
plt.plot(lat,zonalT8,label='$S/S_o$=0.8')
plt.plot(lat,zonalT9,label='$S/S_o$=0.9')
plt.plot(lat,zonalT1,'k',label='$S/S_o$=1.0')
plt.plot(lat,zonalT11,label='$S/S_o$=1.1')

plt.legend(title='Solar Constant Ratio to Present Day', loc='upper right',fontsize=14)

plt.ylabel('Temperature (K)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
#plt.title('Zonally Averaged Temperature Over Land for CAM4 Experiments',fontsize=20)

plt.show()

