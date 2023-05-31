import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
#Read in Data for CMIP6
path='/home/kgemmell/weathering/CESM/CESMcontrol/TS.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
TS1=d.variables['ts'][:,:,:]
TS1=np.nanmean(TS1,axis=0)

path2='/home/kgemmell/weathering/CESM/2xCO2/TS.nc'
d2=Dataset(path2,'r')
TS2=d2.variables['ts'][:,:,:]
TS2=np.nanmean(TS2,axis=0)

path4='/home/kgemmell/weathering/CESM/4xCO2/TS.nc'
d4=Dataset(path4,'r')
TS4=d4.variables['ts'][:,:,:]
TS4=np.nanmean(TS4,axis=0)

zonalT1=np.nanmean(TS1,axis=1) 
zonalT2=np.nanmean(TS2,axis=1)
zonalT4=np.nanmean(TS4,axis=1)

#------------------------------------------------------------------------------#
#Read in Data for CAM4
path1='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d1=Dataset(path1,'r')
lat1=d1.variables['lat'][:]
lon1=d1.variables['lon'][:]

path8='/home/czg/FYSP_clouds_archive/CAM4/cam4_08.nc'
d8=Dataset(path8,'r')
TS8=d8.variables['TS'][:,:,:]
TS8=np.mean(TS8,axis=0)
zonalT8=np.nanmean(TS8,axis=1)

path9='/home/czg/FYSP_clouds_archive/CAM4/cam4_09.nc'
d9=Dataset(path9,'r')
TS9=d9.variables['TS'][:,:,:]
TS9=np.mean(TS9,axis=0)
zonalT9=np.nanmean(TS9,axis=1)

path10='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d10=Dataset(path10,'r')
TS10=d10.variables['TS'][:,:,:]
TS10=np.mean(TS10,axis=0)
zonalT10=np.nanmean(TS10,axis=1)

path11='/home/czg/FYSP_clouds_archive/CAM4/cam4_11.nc'
d11=Dataset(path11,'r')
TS11=d11.variables['TS'][:,:,:]
TS11=np.mean(TS11,axis=0)
zonalT11=np.nanmean(TS11,axis=1)

#------------------------------------------------------------------------------#
#Read in Data for DeepMIP
path2='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d2=Dataset(path2,'r')
lat2=d2.variables['lat'][:]
lon2=d2.variables['lon'][:]

path1e='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2tsDeepMIP.nc'
d1e=Dataset(path1e,'r')
TS1e=d1e.variables['ts'][:,:,:]
TS1e=np.mean(TS1e,axis=0)
zonalT1e=np.nanmean(TS1e,axis=1)

path3='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2tsDeepMIP.nc'
d3=Dataset(path3,'r')
TS3=d3.variables['ts'][:,:,:]
TS3=np.mean(TS3,axis=0)
zonalT3=np.nanmean(TS3,axis=1)

path6='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2tsDeepMIP.nc'
d6=Dataset(path6,'r')
TS6=d6.variables['ts'][:,:,:]
TS6=np.mean(TS6,axis=0)
zonalT6=np.nanmean(TS6,axis=1)

path9e='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2tsDeepMIP.nc'
d9e=Dataset(path9e,'r')
TS9e=d9e.variables['ts'][:,:,:]
TS9e=np.mean(TS9e,axis=0)
zonalT9e=np.nanmean(TS9e,axis=1)

#-------------------------Plotting the Data------------------------------------#
plt.plot(lat,zonalT1,'-k',label='CESM piControl')
plt.plot(lat,zonalT2,'-r',label='CESM 2x $CO_2$')
plt.plot(lat,zonalT4,'-b',label='CESM 4x $CO_2$')

plt.plot(lat1,zonalT8,'-.c',label='$S/S_0=0.8$')
plt.plot(lat1,zonalT9,'-.y',label='$S/S_0=0.9$')
plt.plot(lat1,zonalT10,'-.m',label='$S/S_0=1.0$')
plt.plot(lat1,zonalT11,'-.g',label='$S/S_0=1.1$')

plt.plot(lat2,zonalT1e,'--',label='Eocene 1x $CO_2$')
plt.plot(lat2,zonalT3,'--',label='Eocene 3x $CO_2$')
plt.plot(lat2,zonalT6,'--',label='Eocene 6x $CO_2$')
plt.plot(lat2,zonalT9e,'--',label='Eocene 9x $CO_2$')

plt.legend(title='Case', loc='lower center',fontsize=14)
plt.ylabel('Temperature (K)',fontsize=14)
plt.title('Zonally Averaged Temperature for All Model Runs',fontsize=20)
plt.show()

