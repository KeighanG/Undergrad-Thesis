from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#------------------------------------------------------------------------------#
#Read in Data
path='/home/kgemmell/weathering/CESM/CESMcontrol/TS.nc' #Read in temperature data
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
TS=d.variables['ts'][:,:,:]
lon=np.roll(lon,144) #To match longitude in Hartmann map
TS=np.mean(TS,axis=0)

path1='/home/kgemmell/weathering/CESM/CESMcontrol/Runoff.nc' #Read in runoff data
d1=Dataset(path1,'r')
Runoff=d1.variables['mrro'][:,:,:]
Run=np.mean(Runoff,axis=0)
Run[np.isnan(Run)]=0


#Define Constants
rho=1000 #kg/m3, density of freshwater
Rc=8.314 #J/mol/K
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
bsil= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)

#Removing Oceans from Hartmann data
Land[Land<0]=np.nan
Land[Land>0]=1

#Making Activation Energy Array using homogeneous value
Ea[np.isnan(Land)]=np.nan
Ea=55

#Making bsil array using homogeneous value
bsil[np.isnan(Land)]=np.nan
bsil=0.02

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

regridR=griddata((X.flatten(),Y.flatten()), Run.flatten(), (XI,YI), method='cubic')
regridR[np.isnan(Land)]=np.nan

regridT=griddata((X.flatten(),Y.flatten()), TS.flatten(), (XI,YI), method='cubic')
regridT[np.isnan(Land)]=np.nan

T0=287.1968543493361

#Define variables
R=(regridR/rho)*60*60*24*365.25*1000 #mm/a
C0=280
T=regridT
co2=280


#Equations
Pterm=(co2/C0)**0.5
Tterm=((np.exp((((-Ea*1000)/Rc)*(1/regridT-1/T0)))))
Rterm=R*bsil
F=Pterm*Tterm*Rterm
F[F<=0]=0
Fg=F

#------------------------------------------------------------------------------#
#Conversion to mol C
MMconv=12.011/44.011 #MM of C / MM of CO2
ttomol=(10**6)/12.01 #(g/t)/(g/mol)of C = mol C

F=F*MMconv*ttomol*10**-6 #m^2/km^2, final units =mol km^{-2} yr^{-1}
Fc=np.load('4xLithW.npy') #For setting colourbar
#------------------------------------------------------------------------------#
#Plotting
#Is used to set all non zero values nan which is used to plot a 0 weathering flux field
F_copy=np.copy(F)
F_copy[F!=0]=np.nan

#plotting the 0 value weathering fluxes
fig=plt.contourf(lon1,lat1,F_copy,colors='silver')
#Plotting the main data with a logarithmic colorbar
fig=plt.contourf(lon1,lat1,F,norm=colors.LogNorm(),levels=np.logspace(-5,np.log10(np.nanmax(Fc)+1),100),cmap='viridis')
plt.title('CMIP-6 Control Homogeneous Lithology Weathering')
cbar=plt.colorbar(fig)
cbar.set_label('Weathering Flux (mol km$^{-2}$ yr$^{-1}$)',fontsize=20)
cbar.set_ticks([1E-5,1E-4,1E-3,5E-3,1E-2,5E-2,1E-1,5E-1,1,5,1E+1,50,1E+2])

plt.show()
plt.savefig('CMIPControlHomogeneousMap.png')

Total=Fg*A
final=np.nansum(Total,axis=1)
final=np.nansum(final)/10**12
print(final/44.011) #Prints total global weathering flux in mol/yr


np.save('ControlHomogeneousW',F)

