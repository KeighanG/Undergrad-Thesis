from scipy.interpolate import griddata
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from fastkde import fastKDE
import pylab as plt
from scipy.optimize import curve_fit

#Some threshold values
p_thresh=0.01 #probability density threshold for contour lines: 2 percent of difference between min and max
p_lev=8 #number of probability levels for contour plots
rho=1000 #kg/m3, density of freshwater

#------------------------------------------------------------------------------#
#Read CMIP-6 Control data
path1='/home/kgemmell/weathering/CESM/CESMcontrol/Runoff.nc'
d1=Dataset(path1,'r')
lat=d1.variables['lat'][:]
lon=d1.variables['lon'][:]
lon=np.roll(lon,144)
Runoff=d1.variables['mrro'][:,:,:]

path2='/home/kgemmell/weathering/CESM/CESMcontrol/Precip.nc'
d2=Dataset(path2,'r')
Precip=d2.variables['pr'][:,:,:]

R=np.mean(Runoff,axis=0)
P=np.mean(Precip,axis=0)


#------------------------------------------------------------------------------#
#Load in Hartmann model data
ascii_grid = np.loadtxt("/home/kgemmell/weathering/glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land = np.loadtxt("/home/kgemmell/weathering/glim_wgs84_0point5deg.txt.asc", skiprows=6)
#Removing Oceans from Hartmann data
Land[Land<0]=np.nan
Land[Land>0]=1

lat1=np.linspace(90,-90,360)
lon1=np.linspace(0,360,720)

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(ascii_grid.shape)
for i in range(0,len(lat1)-1):
  for j in range(0,len(lon1)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat1[i]*np.pi/180.)-np.sin(lat1[i+1]*np.pi/180.))*np.abs(lon1[j]-lon1[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

A_weighted_global=A/np.sum(A)
AL=A*Land #Gives area of Land in each grid box
A2=np.nansum(AL) #Gives total area of land
A_weighted=A/A2 #Gives weighting of land in each grid box

#------------------------------------------------------------------------------#
#Regrid the CMIP-6 data
#Define grid spacial resolutions, x1 indicates Hartmann grid
X, Y =np.meshgrid(lon,lat)
XI,YI=np.meshgrid(lon1,lat1)

regridP=griddata((X.flatten(),Y.flatten()), P.flatten(), (XI,YI), method='cubic')
regridR=griddata((X.flatten(),Y.flatten()), R.flatten(), (XI,YI), method='cubic')
regridR[np.isnan(Land)]=np.nan
regridP[np.isnan(Land)]=np.nan

regridP[regridP<0]=0
regridR[regridR<0]=0


Rflat=regridR.flatten()
Pflat=regridP.flatten()

Runfinal=(Rflat/rho)*60*60*24*365.25*1000 #mm/a
Precfinal=(Pflat/rho)*60*60*24*365.25*1000 #mm/a

Runfinala=Runfinal
Precfinala=Precfinal

Runfinal = Runfinal[np.logical_or(~np.isnan(Rflat),~np.isnan(Pflat))]
Precfinal = Precfinal[np.logical_or(~np.isnan(Rflat),~np.isnan(Pflat))]

PDFrun_=np.log10(Runfinal)
PDFprec_=np.log10(Precfinal)


PDFrun = PDFrun_[np.logical_and(np.isfinite(PDFrun_),np.isfinite(PDFprec_))]
PDFprec = PDFprec_[np.logical_and(np.isfinite(PDFrun_),np.isfinite(PDFprec_))]

#Do the self-consistent density estimate
myPDF,axes = fastKDE.pdf(PDFprec,PDFrun)
#Extract the axes from the axis list
v1,v2 = axes

#------------------------------------------------------------------------------#
#Plotting the PDF
plt.scatter(Precfinal,Runfinal,s=0.1,color='grey')
plt.xlabel('Precipitation (mm/a)',fontsize=14)
plt.ylabel('Runoff (mm/a)',fontsize=14)
plt.yscale('log')
plt.xscale('log')
plt.title('Runoff-Precip Probability Density Function from CMIP-6 piControl Model Run',fontsize=20)

#Plot contours of the PDF
plt.contour(10**v1,10**v2,myPDF,levels=np.linspace((1-p_thresh)*np.min(myPDF)+p_thresh*np.max(myPDF),(1-p_thresh)*np.max(myPDF)+p_thresh*np.min(myPDF),p_lev))

#-----------------------------Regression Model---------------------------------#
newX = np.logspace(0, 4, base=10) #Array of precipitation values to plot regression on

#Exponential function, line in log-log space
def myExpFunc(Precfinal, a, b):
    return a * np.power(Precfinal, b)
popt, pcov = curve_fit(myExpFunc, Precfinal, Runfinal)
plt.plot(newX, myExpFunc(newX, *popt), 'r-',
         label='0.0174x$^{1.442}$')
print ('Exponential Fit: y = (a*(x**b))')
print ('\ta = popt[0] = {0}\n\tb = popt[1] = {1}'.format(*popt))

plt.legend(loc='lower right')
plt.show()

