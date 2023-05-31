from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy.ma as ma
#------------------------------------------------------------------------------#
#Read in Data
path='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
Landmask=d.variables['sftlf'][:,:,:]
Land=np.mean(Landmask,axis=0)
LandMask=Land.data
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

path1='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2tsDeepMIP.nc'
d1=Dataset(path1,'r')
TS=d1.variables['ts'][:,:,:]
TS=np.mean(TS,axis=0)
TS[np.isnan(Land)]=np.nan
TS=TS.data


T0=287.1968543493361

path2='/home/kgemmell/weathering/DeepMIP/3xCO2/3xCO2PrecDeepMIP.nc'
d2=Dataset(path2,'r')
Pr=d2.variables['pr'][:,:,:]
Prec=np.mean(Pr,axis=0)
Prec[np.isnan(Land)]=np.nan
Prec=Prec.data

#Define Constants
rho=1000 #kg/m3, density of freshwater
Rc=8.314 #J/mol/K
B=0.65 #From Colbourn
Lvap=2.26*10**6 #J/kg
C0=280
#Making Activation Energy Array, Values from Beusen et al. (2005)
Ea=55
#Making bsil array from Hartmann et al. (2014)
bsil=0.02
co2=840

#Equations
Prec=(Prec/rho)*60*60*24*365.25*1000
R=0.018*Prec**1.441

Pterm=(co2/C0)**0.5
Tterm=((np.exp((((-Ea*1000)/Rc)*(1/TS-1/T0)))))
Rterm=R*bsil
F=Pterm*Tterm*Rterm
F[F<0]=0
Fg=F*LandMask


F9=np.load('9xDeepMip.npy')
#------------------------------------------------------------------------------#
#Conversion to M (10^6) mol C
MMconv=12.011/44.011 #MM of C / MM of CO2
ttomol=(10**6)/12.01 #(g/t)/(g/mol)of C = mol C

F=F*MMconv*ttomol*10**-6*LandMask
F=np.roll(F,72,axis=1)
#------------------------------------------------------------------------------#
#Plotting
#Is used to set all non zero values nan which is used to plot a 0 weathering flux field
F_copy=np.copy(F)
F_copy[F!=0]=np.nan

#plotting the 0 value weathering fluxes
fig=plt.contourf(lon,lat,F_copy,colors='silver')
#Plotting the main data with a logarithmic colorbar
fig=plt.contourf(lon,lat,F,norm=colors.LogNorm(),levels=np.logspace(-4,np.log10(np.nanmax(F9)+1),100),cmap='viridis')
#plt.title('DeepMIP 3x $CO_2$ Weathering (Homogeneous Lithology)')
#cbar=plt.colorbar(fig)
#cbar.set_label('Weathering Flux (M mol km$^{-2}$ a$^{-1}$)',fontsize=20)
#cbar.set_ticks([1E-4,1E-3,5E-3,1E-2,5E-2,1E-1,5E-1,1,5,1E+1,50,1E+2])
plt.show()

Total=Fg*A
final=np.nansum(Total,axis=1)
final=np.nansum(final)/10**12

print(final/44.01)

np.save('3xDeepMip',F)

