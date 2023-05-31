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
#Load in Soil data
path='/home/kgemmell/weathering/Archive/Lithology/Soil/MU_GLOBAL.nc4'
d=Dataset(path,'r')
lat1=d.variables['lat'][:]
lon1=d.variables['lon'][:]
data1=d.variables['MU_GLOBAL'][:,:]
data=ndarray.tolist(data1)
Land=np.array(data)

#Remove Oceans
Land[Land<0]=np.nan
Land[Land>0]=1

#Load in soil shielding array
soil=np.load('SoilArray.npy')

#Plot the data
fig=plt.pcolormesh(lon1,lat1,soil,shading='auto')
plt.clim(0,1)
cbar=plt.colorbar(fig)
cbar.set_label('Soil Shielding Factor',fontsize=20)
cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.show()

