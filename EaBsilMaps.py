from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
#Read in Data
path='/home/kgemmell/weathering/CESM/CESMcontrol/TS.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
TS=d.variables['ts'][:,:,:]
lon=np.roll(lon,144) #To match longitude in Hartmann map
TS=np.mean(TS,axis=0)

path1='/home/kgemmell/weathering/CESM/CESMcontrol/Runoff.nc'
d1=Dataset(path1,'r')
Runoff=d1.variables['mrro'][:,:,:]
Run=np.mean(Runoff,axis=0)
Run[np.isnan(Run)]=0

path2='/home/kgemmell/weathering/CESM/CESMcontrol/Precip.nc'
d2=Dataset(path2,'r')
Precip=d2.variables['pr'][:,:,:]
Precip=np.mean(Precip,axis=0)
Precip[np.isnan(Precip)]=0

#Define Constants
rho=1000 #kg/m3, density of freshwater
Rc=8.314 #J/mol/K
B=0.65 #From Colbourn
Lvap=2.26*10**6 #J/kg

#Define grid spacial resolutions, x1 indicates Hartmann grid
X, Y =np.meshgrid(lon,lat)
lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)
XI,YI=np.meshgrid(lon1,lat1)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Ea= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Eax= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bsil= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Removing Oceans from Hartmann data
Land[Land<0]=np.nan
Land[Land>0]=1

#Making Activation Energy Array, Values from Beusen et al. (2005)
Eax[np.isnan(Land)]=np.nan
Eax[Eax==1]=60
Eax[Eax==2]=50
Eax[Eax==3]=60
Eax[Eax==4]=50
Eax[Eax==5]=60
Eax[Eax==6]=0
Eax[Eax==7]=60
Eax[Eax==8]=60
Eax[Eax==9]=60
Eax[Eax==10]=55 #Volcanic intermediate not listed, put in between acid and basic
Eax[Eax==11]=9999 #Water bodies, no weathering
Eax[Eax==12]=46 #From Hartmann et al. (2014)
Eax[Eax==13]=55
Eax[Eax==14]=0
Eax[Eax==15]=9999 #No data
Eax[Eax==16]=9999 #Ice and Glaciers, no weathering

F_copy=np.copy(Eax)
F_copy[Eax!=9999]=np.nan
#plotting the 0 value weathering fluxes
fig=plt.contourf(lon1,lat1,F_copy,colors='silver')

#Making Activation Energy Array, Values from Beusen et al. (2005)
Ea[np.isnan(Land)]=np.nan
Ea[Ea==1]=60
Ea[Ea==2]=50
Ea[Ea==3]=60
Ea[Ea==4]=50
Ea[Ea==5]=60
Ea[Ea==6]=0
Ea[Ea==7]=60
Ea[Ea==8]=60
Ea[Ea==9]=60
Ea[Ea==10]=55 #Volcanic intermediate not listed, put in between acid and basic
Ea[Ea==11]=np.nan #Water bodies, no weathering
Ea[Ea==12]=46 #From Hartmann et al. (2014)
Ea[Ea==13]=55
Ea[Ea==14]=0
Ea[Ea==15]=np.nan #No data
Ea[Ea==16]=np.nan #Ice and Glaciers, no weathering

fig=plt.contourf(lon1,lat1,Ea,vmax=60)
cbar=plt.colorbar(fig)
cbar.set_label('Activation Energy ($kJ/mol$)',fontsize=20)
plt.show()


#Making bsil array from Hartmann et al. (2014)
bsil[np.isnan(Land)]=np.nan
bsil[bsil==1]=0.021333
bsil[bsil==2]=0.04054
bsil[bsil==3]=0.019699
bsil[bsil==4]=0.04054
bsil[bsil==5]=0.021806
bsil[bsil==6]=0
bsil[bsil==7]=0.020762
bsil[bsil==8]=0.011918
bsil[bsil==9]=0.019307
bsil[bsil==10]=0.04054
bsil[bsil==11]=0
bsil[bsil==12]=0.076876
bsil[bsil==13]=0.019307
bsil[bsil==14]=0
bsil[bsil==15]=0
bsil[bsil==16]=0

fig=plt.contourf(lon1,lat1,bsil)
cbar=plt.colorbar(fig)
cbar.set_label('$b_{sil}$ value',fontsize=20)
plt.show()

