import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
#------------------------------------------------------------------------------#
#Weathering Rates
CESM1=39.5221934361353
CESM2=67.59046715605632
CESM4=209.85702197258075

S08=430.1627701316608
S09=199.29218259945233
S01=57.30944275863661
S011=8.487303642047904

Eocene1=45.68907859599276
Eocene3=144.9780144946899
Eocene6=352.4012370573544
Eocene9=839.3008102350059

#Atmospheric CO2 values
#CESM
CO2CESM1=280
CO2CESM2=560
CO2CESM4=1120

#S/S0 Experimenets
S0CO28=31688
S0CO29=5254
S0CO21=368.9
S0CO211=4.5

#Eocene
CO2E1=280
CO2E3=840
CO2E6=1680
CO2E9=2520


CESMW=np.array([CESM1,CESM2,CESM4])
SS0W=np.array([S08,S09,S01,S011])
EoceneW=np.array([Eocene1,Eocene3,Eocene6,Eocene9])

CESMCO2=np.array([CO2CESM1,CO2CESM2,CO2CESM4])
SS0CO2=np.array([S0CO28,S0CO29,S0CO21,S0CO211])
EoceneCO2=np.array([CO2E1,CO2E3,CO2E6,CO2E9])

plt.scatter(CESMCO2,CESMW,color='red',label='CMIP-6')
plt.scatter(SS0CO2,SS0W,color='green',label='$S/S_0$')
plt.scatter(EoceneCO2,EoceneW,color='blue',label='Eocene')

plt.plot(CESMCO2,CESMW,color='red')
plt.plot(SS0CO2,SS0W,color='green')
plt.plot(EoceneCO2,EoceneW,color='blue')

plt.xscale('log')
plt.yscale('log')
plt.ylabel('Weathering Rate ($10^{12}$ mol $yr^-1$)',fontsize=14)
plt.xlabel('Atmospheric $CO_2$ Concentration (ppm)',fontsize=14)
plt.title('Weathering Rate for all Model Runs',fontsize=20)
plt.legend()
plt.show()

