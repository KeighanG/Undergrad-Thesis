import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
rho=1000
#Read in Data
path='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
#Land fraction
LANDFRAC=d.variables['LANDFRAC'][:,:,:][240:,:,:]
Land=np.mean(LANDFRAC,axis=0)
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


path8='/home/czg/FYSP_clouds_archive/CAM4/cam4_08.nc'
d8=Dataset(path8,'r')
P8=d8.variables['PRECT'][:,:,:]
P8=np.mean(P8,axis=0)
P8[np.isnan(Land)]=np.nan
P8[np.isnan(P8)]=0
P8[P8<0]=0
Prec8=P8*60*60*24*365.25*1000
R8=0.018*Prec8**1.441

path9='/home/czg/FYSP_clouds_archive/CAM4/cam4_09.nc'
d9=Dataset(path9,'r')
P9=d9.variables['PRECT'][:,:,:]
P9=np.mean(P9,axis=0)
P9[np.isnan(P9)]=0
P9[P9<0]=0
P9[np.isnan(Land)]=np.nan
Prec9=P9*60*60*24*365.25*1000
R9=0.018*Prec9**1.441

path1='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d1=Dataset(path1,'r')
P1=d1.variables['PRECT'][:,:,:]
P1=np.mean(P1,axis=0)
P1[np.isnan(P1)]=0
P1[P1<0]=0
P1[np.isnan(Land)]=np.nan
Prec1=P1*60*60*24*365.25*1000
R1=0.018*Prec1**1.441

path11='/home/czg/FYSP_clouds_archive/CAM4/cam4_11.nc'
d11=Dataset(path11,'r')
P11=d11.variables['PRECT'][:,:,:]
P11=np.mean(P11,axis=0)
P11[np.isnan(P11)]=0
P11[P11<0]=0 
P11[np.isnan(Land)]=np.nan
Prec11=P11*60*60*24*365.25*1000
R11=0.018*Prec11**1.441

tot8=R8*A/1000 #m^3/a
tot9=R9*A/1000 #m^3/a
tot1=R1*A/1000 #m^3/a
tot11=R11*A/1000 #m^3/a

#Multiply by area, and then sum
zonalR8=np.nansum(tot8,axis=1)
zonalR9=np.nansum(tot9,axis=1)
zonalR1=np.nansum(tot1,axis=1)
zonalR11=np.nansum(tot11,axis=1)

#Plotting the Data
plt.plot(lat,zonalR8,label='$S/S_o$=0.8')
plt.plot(lat,zonalR9,label='$S/S_o$=0.9')
plt.plot(lat,zonalR1,'k',label='$S/S_o$=1.0')
plt.plot(lat,zonalR11,label='$S/S_o$=1.1')

plt.legend(title='Solar Constant Ratio to Present Day', loc='upper right')

plt.ylabel('Runoff ($m^3$/yr)',fontsize=14)
plt.xlabel('Latitude',fontsize=14)
plt.title('Total Zonal Runoff for CAM4 Experiments',fontsize=16)

plt.show()

print(np.nansum(zonalR1)/10**9)

