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
#Read in Data
path='/home/kgemmell/weathering/CESM/2xCO2/TS.nc'
d=Dataset(path,'r')
lat=d.variables['lat'][:]
lon=d.variables['lon'][:]
TS=d.variables['ts'][:,:,:]
lon[lon>180]=lon[lon>180]-360 #To match longitude in Hartmann map

TS=np.mean(TS,axis=0)

path1='/home/kgemmell/weathering/CESM/2xCO2/Runoff.nc'
d1=Dataset(path1,'r')
Runoff=d1.variables['mrro'][:,:,:]
Run=np.mean(Runoff,axis=0)


#Define Constants
rho=1000 #kg/m3, density of freshwater
Rc=8.314 #J/mol/K
Lvap=2.26*10**6 #J/kg

path='/home/kgemmell/weathering/Archive/Lithology/Soil/MU_GLOBAL.nc4'
d=Dataset(path,'r')
lat1=d.variables['lat'][:]
lon1=d.variables['lon'][:]
lon1=np.roll(lon1,144)
data1=d.variables['MU_GLOBAL'][:,:]
data=ndarray.tolist(data1)
Land=np.array(data)

Land[Land<0]=np.nan
Land[Land>0]=1

#-------------------------Size of Each Gridbox---------------------------------#
r=6371000
A=np.zeros(Land.shape)
for i in range(0,len(lat1)-1):
  for j in range(0,len(lon1)-1):
    A[i,j]=np.pi/180.*r**2*np.abs(np.sin(lat1[i]*np.pi/180.)-np.sin(lat1[i+1]*np.pi/180.))*np.abs(lon1[j]-lon1[j+1])
A[-1,:]=A[0,:]
A[:,-1]=A[:,0]

A_weighted_global=A/np.sum(A)

A=A*Land #Gives area of Land in each grid box
A2=np.nansum(A) #Gives total area of land
A_weighted=A/A2 #Gives weighting of land in each grid box

#Define grid spacial resolutions, x1 indicates soil grid, x2 indicates Hartmann grid
X, Y =np.meshgrid(lon,lat)
XI,YI=np.meshgrid(lon1,lat1)

lat2=np.linspace(90,-90,360)
lon2=np.linspace(-180,180,720)
Xv, Yv=np.meshgrid(lon2,lat2)

#Load in Hartmann model data
ascii_grid = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Ea= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bsil= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bcarb= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
Land1 = np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)


Eas= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bsils= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
bcarbs= np.loadtxt("glim_wgs84_0point5deg.txt.asc", skiprows=6)
#Removing Oceans from Hartmann data
Land1[Land1<0]=np.nan
Land1[Land1>0]=1

#Making Activation Energy Array, Values from Beusen et al. (2005)
Ea[np.isnan(Land1)]=np.nan
Ea=griddata((Xv.flatten(),Yv.flatten()), Ea.flatten(), (XI,YI), method='nearest')
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
bsil[np.isnan(Land1)]=np.nan
bsil=griddata((Xv.flatten(),Yv.flatten()), bsil.flatten(), (XI,YI), method='nearest')
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


#Making Activation Energy Array, Values from Beusen et al. (2005)
Eas[np.isnan(Land1)]=np.nan
Eas[Eas==1]=60
Eas[Eas==2]=50
Eas[Eas==3]=60
Eas[Eas==4]=50
Eas[Eas==5]=60
Eas[Eas==6]=0
Eas[Eas==7]=60
Eas[Eas==8]=60
Eas[Eas==9]=60
Eas[Eas==10]=55 #Volcanic intermediate not listed, put in between acid and basic
Eas[Eas==11]=0 #Water bodies, no wEasthering
Eas[Eas==12]=46 #From Hartmann et al. (2014)
Eas[Eas==13]=55
Eas[Eas==14]=0
Eas[Eas==15]=0 #No data
Eas[Eas==16]=0 #Ice and Glaciers, no wEasthering


#Making bsils array from Hartmann et al. (2014)
bsils[np.isnan(Land1)]=np.nan
bsils[bsils==1]=0.021333
bsils[bsils==2]=0.04054
bsils[bsils==3]=0.019699
bsils[bsils==4]=0.04054
bsils[bsils==5]=0.021806
bsils[bsils==6]=0
bsils[bsils==7]=0.020762
bsils[bsils==8]=0.011918
bsils[bsils==9]=0.019307
bsils[bsils==10]=0.04054
bsils[bsils==11]=0
bsils[bsils==12]=0.076876
bsils[bsils==13]=0.019307
bsils[bsils==14]=0
bsils[bsils==15]=0
bsils[bsils==16]=0

#--------------------Load in list data for Soil--------------------------------#
MUlist=[]
FAO90=[]
shares=[]
Sym90=[]
FAO74=[]
Sym74=[]
with open("soiltable.csv", "r") as csv_file:
    csv_reader = csv.DictReader(csv_file, delimiter=',')
    for lines in csv_reader:
        MUlist.append(int(lines['MU_GLOBAL']))
        shares.append(float(lines['SHARE']))
        if lines['SU_CODE90']== '':
            FAO90.append(np.nan)
        else:
            FAO90.append(int(lines['SU_CODE90']))
        if lines['SU_SYM90']=='':
            Sym90.append(np.nan)
        else:
            Sym90.append(lines['SU_SYM90'])
        if lines['SU_CODE74']== '':
            FAO74.append(np.nan)
        else:
            FAO74.append(int(lines['SU_CODE74']))
        if lines['SU_SYM74']=='':
            Sym74.append(np.nan)
        else:
            Sym74.append(lines['SU_SYM74'])

MUarr=np.array(MUlist)
sharesarr=np.array(shares)
FAO90arr=np.array(FAO90)
FAO74arr=np.array(FAO74)
#-------------------------Making Unique Array----------------------------------#
MUunique= []
for i in MUlist:
    if i not in MUunique:
        MUunique.append(i)
MUuniquearr=np.array(MUunique)

#-------------------Sublists of Shares and Soil Classes------------------------#
Allshares=[]
for i in range(len(MUuniquearr)):
    Allshares.append(sharesarr[MUarr==MUuniquearr[i]])

allFAO90=[]
for i in range(len(MUuniquearr)):
    allFAO90.append(FAO90arr[MUarr==MUuniquearr[i]])

allFAO74=[]
for i in range(len(MUuniquearr)):
    allFAO74.append(FAO74arr[MUarr==MUuniquearr[i]])

#-------------------------Making Sheild Dict-----------------------------------#
sheild_dict90={}
for i in range(len(FAO90arr)):
    sheild_dict90[FAO90arr[i]] = FAO90arr[i]

for i in sheild_dict90:
        sheild_dict90[i]=1

#FAO90 sheildings
#GLEYSOLS
sheild_dict90[9.0]=0.1
sheild_dict90[10.0]=0.1
sheild_dict90[11.0]=0.1
sheild_dict90[12.0]=0.1
sheild_dict90[14.0]=0.1
sheild_dict90[15.0]=0.1
sheild_dict90[14.0]=0.1
sheild_dict90[16.0]=0.1
sheild_dict90[17.0]=0.1

#ACRISOLS
sheild_dict90[19.0]=0.1
sheild_dict90[20.0]=0.1
sheild_dict90[21.0]=0.1
sheild_dict90[22.0]=0.1
sheild_dict90[23.0]=0.1

#FERRALSOLS
sheild_dict90[71.0]=0.1
sheild_dict90[72.0]=0.1
sheild_dict90[73.0]=0.1
sheild_dict90[74.0]=0.1
sheild_dict90[75.0]=0.1
sheild_dict90[76.0]=0.1
sheild_dict90[77.0]=0.1

#HISTOSOLS
sheild_dict90[86.0]=0.1
sheild_dict90[87.0]=0.1
sheild_dict90[88.0]=0.1
sheild_dict90[89.0]=0.1
sheild_dict90[90.0]=0.1
sheild_dict90[91.0]=0.1

#LIXISOLS
sheild_dict90[115.0]=0.1
sheild_dict90[116.0]=0.1
sheild_dict90[117.0]=0.1
sheild_dict90[118.0]=0.1
sheild_dict90[119.0]=0.1
sheild_dict90[120.0]=0.1

#NITISOLS
sheild_dict90[121.0]=0.1
sheild_dict90[122.0]=0.1
sheild_dict90[123.0]=0.1
sheild_dict90[124.0]=0.1

sheild_dict74={}
for i in range(len(FAO74arr)):
    sheild_dict74[FAO74arr[i]] = FAO74arr[i]

for i in sheild_dict74:
        sheild_dict74[i]=1
#FAO74 Sheildings
#GLEYSOLS
sheild_dict74[1.0]=0.1
sheild_dict74[2.0]=0.1
sheild_dict74[3.0]=0.1
sheild_dict74[4.0]=0.1
sheild_dict74[5.0]=0.1
sheild_dict74[6.0]=0.1
sheild_dict74[7.0]=0.1
sheild_dict74[8.0]=0.1

#ACRISOLS
sheild_dict74[104.0]=0.1
sheild_dict74[105.0]=0.1
sheild_dict74[106.0]=0.1
sheild_dict74[107.0]=0.1
sheild_dict74[108.0]=0.1
sheild_dict74[109.0]=0.1

#NITISOLS
sheild_dict74[110.0]=0.1
sheild_dict74[111.0]=0.1
sheild_dict74[112.0]=0.1
sheild_dict74[113.0]=0.1

#FERRALSOLS
sheild_dict74[114.0]=0.1
sheild_dict74[115.0]=0.1
sheild_dict74[116.0]=0.1
sheild_dict74[117.0]=0.1
sheild_dict74[118.0]=0.1
sheild_dict74[119.0]=0.1
sheild_dict74[120.0]=0.1

#HISTOSOLS
sheild_dict74[121.0]=0.1
sheild_dict74[122.0]=0.1
sheild_dict74[123.0]=0.1
sheild_dict74[124.0]=0.1


#-------------------Apply sheilding dict to FAO90 Sublist----------------------#
allsheildings=[]
for i in range(0,len(allFAO90)):
    temparray=np.zeros(allFAO90[i].shape,type(allFAO90[i]))
    for j in range(0,len(allFAO90[i])):
        if not np.isnan(allFAO90[i][j]):
            temparray[j]=sheild_dict90[allFAO90[i][j]]
        else:
            temparray[j]=sheild_dict74[allFAO74[i][j]]
    allsheildings.append(temparray)

#---------------------Sheilding Equation to apply------------------------------#
#Apply to the sublists to get a single number to represent the effective
#sheilding in each grid box.

sheildingeff=[]
for i in range(0,len(allsheildings)):
    sheildingeff.append(np.sum((Allshares[i]/100)*allsheildings[i]))

#-------------------Making Effective Sheilding Dict----------------------------#
sheildeff_dict={}
for i in range(len(MUunique)):
    sheildeff_dict[MUunique[i]]=sheildingeff[i]

#-------------------------Load in actual data----------------------------------#
rows=len(lat1)
cols=len(lon1)

arr=[[-1 for i in range(cols)] for j in range(rows)]
for i in range(0,len(data)):
    for j in range(0,len(data[i])):
        if data[i][j] > 0:
            arr[i][j] = sheildeff_dict[data[i][j]]

soilarr=np.array(arr)
soilarr[soilarr<0]=np.nan

regridR=griddata((X.flatten(),Y.flatten()), Run.flatten(), (XI,YI), method='cubic')
regridR[np.isnan(Land)]=np.nan

regridT=griddata((X.flatten(),Y.flatten()), TS.flatten(), (XI,YI), method='cubic')
regridT0=griddata((X.flatten(),Y.flatten()), TS.flatten(), (XI,YI), method='cubic')
regridT[np.isnan(Land)]=np.nan

T0=287.1968543493361


#Define variables
R=(regridR/rho)*60*60*24*365.25*1000 #mm/a
C0=280
T=regridT
co2=560

#Equations
Pterm=(co2/C0)**0.5
Tterm=((np.exp((((-Ea*1000)/Rc)*(1/regridT-1/T0)))))
Tterm[Tterm<0.000001]=0 #To remove ice and glaciers
Rterm=R*bsil
F=Pterm*Tterm*Rterm*soilarr
F[F<0]=0
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
fig=plt.contourf(lon1,lat1,F,norm=colors.LogNorm(),levels=np.logspace(-3.5,np.log10(np.nanmax(Fc)+1),200),cmap='viridis')
plt.title('CESM Abrubt 2x $CO_2$ Weathering with Soil Shielding')
#cbar=plt.colorbar(fig)
#cbar.set_label('Weathering Flux (M mol km$^{-2}$ a$^{-1}$)',fontsize=20)
#cbar.set_ticks([1E-4,1E-3,5E-3,1E-2,5E-2,1E-1,5E-1,1,5,1E+1,50,1E+2])

plt.show()
#plt.savefig('2xSoilMap.png')

Flux=np.nanmean(Fg,axis=1)
final=np.nansum(Fg*A_weighted)
GlobalW=final*A2/10**15 #conversion from m^2 to km^2, and to Gt
#Conversion from Gt CO2 to mol C
MMconv=12.011/44.011
Gttomol=(10**15)/12.011 #(g/Gt)/(g/mol)=mol
GlobalW=GlobalW*MMconv*Gttomol/10**12

print(GlobalW)

np.save('2xSoilW',F)

