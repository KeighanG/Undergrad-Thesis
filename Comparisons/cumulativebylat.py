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
#Load in CMIP-6 Data
lat1=np.linspace(90,-90,360)
CMIP1=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIP1lat.npy')
CMIP2=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIP2lat.npy')
CMIP4=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIP4lat.npy')
#------------------------------------------------------------------------------#
#Load in S/S0 Data
SS08=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS08lat.npy')
SS09=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS09lat.npy')
SS01=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS01lat.npy')
SS011=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS011lat.npy')

#------------------------------------------------------------------------------#
#Load in Eocene Data
#Read in Data
path2='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d2=Dataset(path2,'r')
lat2=d2.variables['lat'][:]
lat2=lat2.data
lat2=np.flip(lat2)

eocene1x=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocene1xlat.npy')
eocene3x=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocene3xlat.npy')
eocene6x=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocene6xlat.npy')
eocene9x=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocene9xlat.npy')

#-------------------------Plotting the Data------------------------------------#
plt.plot(np.flip(lat1[0:180]),CMIP1,'-k',label='CMIP-6 piControl')
plt.plot(np.flip(lat1[0:180]),CMIP2,'-r',label='CMIP-6 2x $CO_2$')
plt.plot(np.flip(lat1[0:180]),CMIP4,'-b',label='CMIP-6 4x $CO_2$')

plt.plot(np.flip(lat1[0:180]),SS08,'-.c',label='$S/S_0=0.8$')
plt.plot(np.flip(lat1[0:180]),SS09,'-.y',label='$S/S_0=0.9$')
plt.plot(np.flip(lat1[0:180]),SS01,'-.m',label='$S/S_0=1.0$')
plt.plot(np.flip(lat1[0:180]),SS011,'-.g',label='$S/S_0=1.1$')

plt.plot(np.flip(lat2[0:48]),eocene1x,'--',label='Eocene 1x $CO_2$')
plt.plot(np.flip(lat2[0:48]),eocene3x,'--',label='Eocene 3x $CO_2$')
plt.plot(np.flip(lat2[0:48]),eocene6x,'--',label='Eocene 6x $CO_2$')
plt.plot(np.flip(lat2[0:48]),eocene9x,'--',label='Eocene 9x $CO_2$')

plt.legend(title='Experiment', loc='lower right',fontsize=16)
plt.xlabel('Equator to Pole Latitude (+ or - up to 90)',fontsize=16)
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.ylabel('Cumulative Flux (%)',fontsize=16)
#plt.title('Cumulative Flux by Latitude Starting at Equator',fontsize=20)
plt.show()

