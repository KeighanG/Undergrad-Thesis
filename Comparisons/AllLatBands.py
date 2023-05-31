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
#Load in Present Day Latitudes
lat1=np.linspace(90,-90,360)

#Load in Eocene Latitudes
#Read in Data
path2='/home/kgemmell/weathering/DeepMIP/landmask.nc'
d2=Dataset(path2,'r')
lat2=d2.variables['lat'][:]
lat2=lat2.data
lat2=np.flip(lat2)

#------------------------------------------------------------------------------#
#CMIP-6 Lat Bands
CMIP115=67.88259164672945
CMIP130=20.96119379408244
CMIP190=11.156214559188104

CMIP215=67.46090334495976
CMIP230=20.995638602121776
CMIP290=11.543458052918472

CMIP415=60.75871412175332
CMIP430=23.351418568019216
CMIP490=15.889867310227466

#------------------------------------------------------------------------------#
#S/S0 Lat Bands
SS0815=64.92950652718855
SS0830=24.957245436638278
SS0890=10.113248036173164

SS0915=64.83134884206343
SS0930=25.943159089513607
SS0990=9.22549206842297

SS0115=64.69511151181196
SS0130=25.857003468946647
SS0190=9.447885019241403

SS01115=67.35162798208505
SS01130=25.45773264487814
SS01190=7.190639373036819

#------------------------------------------------------------------------------#
#Eocene Lat Bands
E115=67.76513536078286
E130=15.706565441007267
E190=16.52829919820988

E315=65.92052279879135
E330=13.597148186042793
E390=20.48232901516585

E615=66.91434989633655
E630=11.72690771363261
E690=17.59134766139259

E915=73.79703939696405
E930=9.156441733899765
E990=17.046518869136175

#------------------------------------------------------------------------------#
#Temperature Gradients
CMIPT1=61.3
CMIPT2=59.05
CMIPT4=53.2

SS0T8=50.95
SS0T9=55.25
SS0T1=59
SS0T11=63.45

ET1=42.5
ET3=34.625
ET6=32.375
ET9=28.5

#------------------------------------------------------------------------------#
#Plotting the data

CMIPTemps=np.array([CMIPT1,CMIPT2,CMIPT4])
SS0Temps=np.array([SS0T8,SS0T9,SS0T1,SS0T11])
ETemps=np.array([ET1,ET3,ET6,ET9])

CMIPfifteen=np.array([CMIP115,CMIP215,CMIP415])
SS0fifteen=np.array([SS0815,SS0915,SS0115,SS01115])
Efifteen=np.array([E115,E315,E615,E915])

CMIPthirty=np.array([CMIP130,CMIP230,CMIP430])
SS0thirty=np.array([SS0830,SS0930,SS0130,SS01130])
Ethirty=np.array([E130,E330,E630,E930])

CMIPninty=np.array([CMIP190,CMIP290,CMIP490])
SS0ninty=np.array([SS0890,SS0990,SS0190,SS01190])
Eninty=np.array([E190,E390,E690,E990])

plt.scatter(CMIPTemps,CMIPfifteen,color='red')
plt.scatter(SS0Temps,SS0fifteen,color='green')
plt.scatter(ETemps,Efifteen,color='blue')

line1,=plt.plot(CMIPTemps,CMIPfifteen,'-',color='black',label='Equator-15')
line2,=plt.plot(CMIPTemps,CMIPthirty,'--',color='black',label='15-30')
line3,=plt.plot(CMIPTemps,CMIPninty,':',color='black',label='30-Poles')
first_legend=plt.legend(handles=[line1,line2,line3],loc='lower center')

plt.plot(CMIPTemps,CMIPfifteen,'-',color='red',label='CMIP-6')
plt.plot(SS0Temps,SS0fifteen,'-',color='green',label='$S/S_0$')
plt.plot(ETemps,Efifteen,'-',color='blue',label='Eocene')

plt.scatter(CMIPTemps,CMIPthirty,color='red')
plt.scatter(SS0Temps,SS0thirty,color='green')
plt.scatter(ETemps,Ethirty,color='blue')
plt.plot(CMIPTemps,CMIPthirty,'--',color='red')
plt.plot(SS0Temps,SS0thirty,'--',color='green')
plt.plot(ETemps,Ethirty,'--',color='blue')

plt.scatter(CMIPTemps,CMIPninty,color='red')
plt.scatter(SS0Temps,SS0ninty,color='green')
plt.scatter(ETemps,Eninty,color='blue')
plt.plot(CMIPTemps,CMIPninty,':',color='red')
plt.plot(SS0Temps,SS0ninty,':',color='green')
plt.plot(ETemps,Eninty,':',color='blue')

leg2=plt.legend(fontsize=18,loc='center')

plt.xlabel('Equator-Pole Temp. Difference (K)',fontsize=18)
plt.ylabel('Cumulative Weathering Flux Contribution (%)',fontsize=18)

plt.show()

