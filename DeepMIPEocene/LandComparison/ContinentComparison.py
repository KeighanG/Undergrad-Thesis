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
Land=np.roll(Land,72,axis=1)

path1='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d1=Dataset(path1,'r')
lat1=d.variables['lat'][:]
lon1=d.variables['lon'][:]

#Land fraction
LANDFRAC1=d1.variables['LANDFRAC'][:,:,:][240:,:,:]
Land1=np.mean(LANDFRAC1,axis=0)
Landmask1=Land1
Land1[Land1<=0]=np.nan
Land1=np.roll(Land1,72,axis=1)

plt.contourf(lon,lat,Land1,cmap='Greys')
plt.contourf(lon,lat,Land,cmap='Reds')
#plt.title('Present Day (Black) and Eocene (Red) Continetal Configuration Comparison',fontsize=20)
plt.show()

