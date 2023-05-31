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
CMIPLand1=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIPLand1.npy')
CMIPLand2=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIPLand2.npy')
CMIPLand4=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIPLand4.npy')

CMIPcont1=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIPcont1.npy')
CMIPcont2=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIPcont2.npy')
CMIPcont4=np.load('/home/kgemmell/weathering/CESM/CumFluxLand/CMIPcont4.npy')

#------------------------------------------------------------------------------#
#Load in S/S0 Data
SS0Land8=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0Land8.npy')
SS0Land9=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0Land9.npy')
SS0Land1=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0Land1.npy')
SS0Land11=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0Land11.npy')

SS0cont8=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0cont8.npy')
SS0cont9=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0cont9.npy')
SS0cont1=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0cont1.npy')
SS0cont11=np.load('/home/kgemmell/weathering/S:S0/CumulativeFluxes/SS0cont11.npy')

#------------------------------------------------------------------------------#
#Load in Eocene Data
eoceneLand1=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eoceneLand1.npy')
eoceneLand3=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eoceneLand3.npy')
eoceneLand6=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eoceneLand6.npy')
eoceneLand9=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eoceneLand9.npy')

eocenecont1=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocenecont1.npy')
eocenecont3=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocenecont3.npy')
eocenecont6=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocenecont6.npy')
eocenecont9=np.load('/home/kgemmell/weathering/DeepMIP/CumulativeFluxes/eocenecont9.npy')

#------------------------------------------------------------------------------#
#Plotting the Data

plt.plot(CMIPLand1,CMIPcont1,'-k',label='CMIP-6 piControl')
plt.plot(CMIPLand2,CMIPcont2,'-r',label='CMIP-6 2x $CO_2$')
plt.plot(CMIPLand4,CMIPcont4,'-b',label='CMIP-6 4x $CO_2$')

plt.plot(SS0Land8,SS0cont8,'-.c',label='$S/S_0=0.8$')
plt.plot(SS0Land9,SS0cont9,'-.y',label='$S/S_0=0.9$')
plt.plot(SS0Land1,SS0cont1,'-.m',label='$S/S_0=1.0$')
plt.plot(SS0Land11,SS0cont11,'-.g',label='$S/S_0=1.1$')

plt.plot(eoceneLand1,eocenecont1,'--',label='Eocene 1x $CO_2$')
plt.plot(eoceneLand3,eocenecont3,'--',label='Eocene 3x $CO_2$')
plt.plot(eoceneLand6,eocenecont6,'--',label='Eocene 6x $CO_2$')
plt.plot(eoceneLand9,eocenecont9,'--',label='Eocene 9x $CO_2$')

plt.ylabel('Accumulated Weathering Flux (%)',fontsize=14)
plt.xlabel('Accumulated Land Area (%)',fontsize=14)
#plt.title('Cumulative Flux-Area Rating Curves for CMIP-6 Experiments',fontsize=20)
plt.legend(title='Case', loc='lower right',fontsize=14)
plt.grid(color='lightgrey', linestyle='dotted', linewidth=0.5)
plt.show()

