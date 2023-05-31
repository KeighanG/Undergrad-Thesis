import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
#------------------------------------------------------------------------------#
#Globally Averaged Temperatures
#CESM
TCESM1=287.1968773085694
TCESM2=289.57790268053
TCESM4=298.2636710030565

alkCESM1=2123.4533116873254
alkCESM2=2685.123510027255
alkCESM4=3469.7801167944053

#S/S0 Experimenets
TS08=287.5366857056836
TS09=287.4008038806125
TS01=287.41047659254707
TS011=287.40237836824383

alkS08=13190.063895634932
alkS09=6403.080071202527
alkS01=2324.0004969578463
alkS011=1076.4188107805016

#Eocene
TE1=292.50103944519566
TE3=298.9591694295133
TE6=303.6137493797363
TE9=309.2195212402544

alkE1=2123.4533116873254
alkE3=3112.9432146607555
alkE6=4059.2786262261425
alkE9=4765.149388512441


CESMtemp=np.array([TCESM1,TCESM2,TCESM4])
SS0temp=np.array([TS08,TS09,TS01,TS011])
Eocenetemp=np.array([TE1,TE3,TE6,TE9])

CESMalk=np.array([alkCESM1,alkCESM2,alkCESM4])
SS0alk=np.array([alkS08,alkS09,alkS01,alkS011])
Eocenealk=np.array([alkE1,alkE3,alkE6,alkE9])


plt.scatter(CESMtemp,CESMalk,color='red',label='CMIP-6')
plt.scatter(SS0temp,SS0alk,color='green',label='S/S0')
plt.scatter(Eocenetemp,Eocenealk,color='blue',label='Eocene')

plt.plot(CESMtemp,CESMalk,color='red')
plt.plot(SS0temp,SS0alk,color='green')
plt.plot(Eocenetemp,Eocenealk,color='blue')

plt.ylabel('Alkalinity (umol $kg^{-1}$)',fontsize=14)
plt.xlabel('Temperature (K)',fontsize=14)
plt.title('Alkalinity for all Model Runs',fontsize=20)
plt.legend()
plt.show()

