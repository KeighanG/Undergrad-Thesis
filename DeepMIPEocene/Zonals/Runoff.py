from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
rho=1000
#Read in Data
path='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
Landmask=d.variables['sftlf'][:,:,:]
Land=np.mean(Landmask,axis=0)
LandFrac=Land
Land[Land<=0]=np.nan
#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(Land.shape)
for i in range(0,len(lat)-1):
  for j in range(0,len(lon)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat[i]*np.pi/180.)-np.sin(lat[i+1]*np.pi/180.))*np.abs(lon[j]-lon[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

A_weighted_global=A/np.sum(A)
A=A*Land #Gives area of Land in each grid box
A2=np.nansum(A) #Gives total area of land
A_weighted=A/A2 #Gives weighting of land in each grid box

#Runoffs
path1='/home/kgemmell/weathering/DeepMIP/1xCO2/1xCO2PrecDeepMIP.nc'
d1=Dataset(path1,'r')
P1=d1.variables['pr'][:,:,:]
P1=np.mean(P1,axis=0)
P1[np.isnan(Land)]=np.nan
Prec1=(P1/rho)*60*60*24*365.25*1000
R1=0.018*Prec1**1.441

path3='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2PrecDeepMIP.nc'
d3=Dataset(path3,'r')
P3=d3.variables['pr'][:,:,:]
P3=np.mean(P3,axis=0)
P3[np.isnan(Land)]=np.nan
Prec3=(P3/rho)*60*60*24*365.25*1000
R3=0.018*Prec3**1.441


path6='/home/kgemmell/weathering/DeepMIP/6xCO2/6xCO2PrecDeepMIP.nc'
d6=Dataset(path6,'r')
P6=d6.variables['pr'][:,:,:]
P6=np.mean(P6,axis=0)
P6[np.isnan(Land)]=np.nan
Prec6=(P6/rho)*60*60*24*365.25*1000
R6=0.018*Prec6**1.441


path9='/home/kgemmell/weathering/DeepMIP/9xCO2/9xCO2PrecDeepMIP.nc'
d9=Dataset(path9,'r')
P9=d9.variables['pr'][:,:,:]
P9=np.mean(P9,axis=0)
P9[np.isnan(Land)]=np.nan
Prec9=(P9/rho)*60*60*24*365.25*1000
R9=0.018*Prec9**1.441


tot1x=R1*A*LandFrac/1000 #m^3/a
tot3x=R3*A*LandFrac/1000 #m^3/a
tot6x=R6*A*LandFrac/1000 #m^3/a
tot9x=R9*A*LandFrac/1000 #m^3/a

#Multiply by area, and then sum
zonalR1=np.nansum(tot1x,axis=1)
zonalR3=np.nansum(tot3x,axis=1)
zonalR6=np.nansum(tot6x,axis=1)
zonalR9=np.nansum(tot9x,axis=1)

#Plotting the Data
plt.plot(lat,zonalR1,label='1x $CO_2$')
plt.plot(lat,zonalR3,label='3x $CO_2$')
plt.plot(lat,zonalR6,label='6x $CO_2$')
plt.plot(lat,zonalR9,label='9x $CO_2$')

plt.legend(title='Case', loc='upper right',fontsize=14)

plt.ylabel('Runoff ($m^3$/yr)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
#plt.title('Total Zonal Runoff for DeepMIP Eocene Data',fontsize=16)

plt.show()

