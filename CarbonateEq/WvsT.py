import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
#------------------------------------------------------------------------------#
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

#Globally Averaged Temperatures
#CESM
TCESM1=287.1968773085694
TCESM2=289.57790268053
TCESM4=298.2636710030565

#S/S0 Experimenets
TS08=287.5366857056836
TS09=287.4008038806125
TS01=287.41047659254707
TS011=287.40237836824383

#Eocene
TE1=292.50103944519566
TE3=298.9591694295133
TE6=303.6137493797363
TE9=309.2195212402544

Weatherings=np.array([CESM1,CESM2,CESM4,S08,S09,S01,S011,Eocene1,Eocene3,Eocene6,Eocene9])
Temps=np.array([TCESM1,TCESM2,TCESM4,TS08,TS09,TS01,TS011,TE1,TE3,TE6,TE9])

CESMW=np.array([CESM1,CESM2,CESM4])
SS0W=np.array([S08,S09,S01,S011])
EoceneW=np.array([Eocene1,Eocene3,Eocene6,Eocene9])

CESMtemp=np.array([TCESM1,TCESM2,TCESM4])
SS0temp=np.array([TS08,TS09,TS01,TS011])
Eocenetemp=np.array([TE1,TE3,TE6,TE9])

plt.scatter(CESMtemp,CESMW,color='red',label='CMIP-6')
plt.scatter(SS0temp,SS0W,color='green',label='$S/S_0$')
plt.scatter(Eocenetemp,EoceneW,color='blue',label='Eocene')

plt.plot(CESMtemp,CESMW,color='red')
plt.plot(SS0temp,SS0W,color='green')
plt.plot(Eocenetemp,EoceneW,color='blue')

plt.ylabel('Weathering Rate ($10^{12}$ mol $yr^-1$)',fontsize=14)
plt.xlabel('Temperature (K)',fontsize=14)
#plt.title('Weathering Rate for all Model Runs',fontsize=20)
plt.legend()
plt.show()

