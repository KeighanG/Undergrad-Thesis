from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#------------------------------------------------------------------------------#
#Read in our model data
path='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
lon=np.roll(lon,72) #To match longitude in Hartmann map

#Define grid spacial resolutions, x1 indicates Hartmann grid
X, Y =np.meshgrid(lon,lat)
lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)
XI,YI=np.meshgrid(lon1,lat1)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Ea= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bsil= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bcarb= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Removing Oceans from Hartmann data
Land[Land<0]=np.nan
Land[Land>0]=1

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
Ea[Ea==11]=0 #Water bodies, no weathering
Ea[Ea==12]=46 #From Hartmann et al. (2014)
Ea[Ea==13]=55
Ea[Ea==14]=0
Ea[Ea==15]=0 #No data
Ea[Ea==16]=9999 #Ice and Glaciers, no weathering

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

#Making bcarb array from Hartmann et al. (2014)
bcarb[np.isnan(Land)]=np.nan
bcarb[bcarb==1]=0
bcarb[bcarb==2]=0
bcarb[bcarb==3]=0.0081
bcarb[bcarb==4]=0
bcarb[bcarb==5]=0.035915
bcarb[bcarb==6]=0.151243
bcarb[bcarb==7]=0
bcarb[bcarb==8]=0.020118
bcarb[bcarb==9]=0.007603
bcarb[bcarb==10]=0
bcarb[bcarb==11]=0
bcarb[bcarb==12]=0
bcarb[bcarb==13]=0.007603
bcarb[bcarb==14]=0.151243
bcarb[bcarb==15]=0
bcarb[bcarb==16]=0

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(ascii_grid.shape)
for i in range(0,len(lat1)-1):
  for j in range(0,len(lon1)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat1[i]*np.pi/180.)-np.sin(lat1[i+1]*np.pi/180.))*np.abs(lon1[j]-lon1[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

A_weighted_global=A/np.sum(A)

A=A*Land #Gives area of Land in each grid box
A2=np.nansum(A) #Gives total area of land
A_weighted=A/A2 #Gives weighting of land in each grid box
#-------------------------INITIALIZE VARIABLES---------------------------------#
#Read in data for initializing variables
path='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d=Dataset(path,'r')

#CO20
co2vmr=d.variables['co2vmr'][:]
co20=co2vmr[0]*1E+6
#DONE INITIALIZING VARIABLES#

#Define Constants
rho=1000 #kg/m3, density of freshwater
C0=280
Ea=Ea#J/mol
Rc=8.314 #J/mol/K
Lvap=2.26*10**6 #J/kg
T0=287.1968543493361
#------------------------------------------------------------------------------#
#0.8, for setting Colorbar
path8='/home/czg/FYSP_clouds_archive/CAM4/cam4_08.nc'
d8=Dataset(path8,'r')

#CO2
co2vmr8=d8.variables['co2vmr'][:]
co28=co2vmr8[0]*1000000

#Precipitation
PRECT8=d8.variables['PRECT'][:,:,:][240:,:,:]
Precip8=np.mean(PRECT8,axis=0)
regridPrecip8=griddata((X.flatten(),Y.flatten()), Precip8.flatten(), (XI,YI), method='cubic')
regridPrecip8[np.isnan(Land)]=np.nan


#Surface Temperature
Temp8=d8.variables['TS'][:,:,:][240:,:,:]
TS8=np.nanmean(Temp8,axis=0)
regridT8=griddata((X.flatten(),Y.flatten()), TS8.flatten(), (XI,YI), method='cubic')
regridT8[np.isnan(Land)]=np.nan

#Define variables
Prec8=regridPrecip8*60*60*24*365.25*1000
R8=0.018*Prec8**1.441
C8=co28
T8=regridT8

#Equations
Pterm8=(co28/C0)**0.5
Tterm8=((np.exp((((-Ea*1000)/Rc)*(1/regridT8-1/T0)))))
Rterm8=R8*(bsil)
F8=Pterm8*Tterm8*Rterm8
F8[F8<0]=0

#------------------------------------------------------------------------------#
#Load in data for variables
path='/home/czg/FYSP_clouds_archive/CAM4/cam4_09.nc'
d=Dataset(path,'r')

#CO2
co2vmr=d.variables['co2vmr'][:]
co2=co2vmr[0]*1000000

#Precip
PRECT=d.variables['PRECT'][:,:,:][240:,:,:]
Precip=np.mean(PRECT,axis=0)


regridPrecip=griddata((X.flatten(),Y.flatten()), Precip.flatten(), (XI,YI), method='cubic')
regridPrecip[np.isnan(Land)]=np.nan


#Surface Temperature
Temp=d.variables['TS'][:,:,:][240:,:,:]
TS=np.nanmean(Temp,axis=0)
regridT=griddata((X.flatten(),Y.flatten()), TS.flatten(), (XI,YI), method='cubic')

regridT[np.isnan(Land)]=np.nan

#Define variables
Prec=regridPrecip*60*60*24*365.25*1000
R=0.018*Prec**1.441
C=co2
T=regridT

#Equations
#Pterm=((2*(co2/C0))/(1+(co2/C0)))**0.4
Pterm=(co2/C0)**0.5
Tterm=((np.exp((((-Ea*1000)/Rc)*(1/regridT-1/T0)))))
Tterm[Tterm<0.01]=0
Rterm=R*(bsil)
F=Pterm*Tterm*Rterm
F[F<0]=0

Fg=F

#------------------------------------------------------------------------------#
#Conversion to mol C
MMconv=12.011/44.011 #MM of C / MM of CO2
ttomol=(10**6)/12.01 #(g/t)/(g/mol)of C = mol C

F=F*MMconv*ttomol*10**-6
F8=F8*MMconv*ttomol*10**-6

#------------------------------------------------------------------------------#
#Plotting
#Is used to set all non zero values nan which is used to plot a 0 weathering flux field
F_copy=np.copy(F)
F_copy[F!=0]=np.nan

#plotting the 0 value weathering fluxes
fig=plt.contourf(lon1,lat1,F_copy,colors='silver')
#Plotting the main data with a logarithmic colorbar
fig=plt.contourf(lon1,lat1,F,norm=colors.LogNorm(),levels=np.logspace(-3.5,np.log10(np.nanmax(F8)+1),200),cmap='viridis')
#plt.title('$S/S_o$=0.9 Weathering Flux',fontsize=20)
plt.show()

Total=Fg*A
final=np.nansum(Total,axis=1)
final=np.nansum(final)/10**12
Flux=final/44.011
print(Flux)

np.save('Flux9',Fg)

