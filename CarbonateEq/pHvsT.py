import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
#------------------------------------------------------------------------------#
#Globally Averaged Temperatures and pH
#CESM
TCESM1=287.1968773085694
TCESM2=289.57790268053
TCESM4=298.2636710030565

pHCESM1=8.140464570949558
pHCESM2=7.9828026570636546
pHCESM4=7.822883503618185

#S/S0 Experimenets
TS08=287.5366857056836
TS09=287.4008038806125
TS01=287.41047659254707
TS011=287.40237836824383

pHS08=7.007475696381504
pHS09=7.45588727065919
pHS01=8.077987180655779
pHS011=9.0540886079346

#Eocene
TE1=292.50103944519566
TE3=298.9591694295133
TE6=303.6137493797363
TE9=309.2195212402544

pHE1=8.140464570949558
pHE3=7.889559669851588
pHE6=7.72810613500863
pHE9=7.632311156215022


CESMtemp=np.array([TCESM1,TCESM2,TCESM4])
SS0temp=np.array([TS08,TS09,TS01,TS011])
Eocenetemp=np.array([TE1,TE3,TE6,TE9])

CESMpH=np.array([pHCESM1,pHCESM2,pHCESM4])
SS0pH=np.array([pHS08,pHS09,pHS01,pHS011])
EocenepH=np.array([pHE1,pHE3,pHE6,pHE9])


plt.scatter(CESMtemp,CESMpH,color='red',label='CMIP-6')
plt.scatter(SS0temp,SS0pH,color='green',label='S/S0')
plt.scatter(Eocenetemp,EocenepH,color='blue',label='Eocene')

plt.plot(CESMtemp,CESMpH,color='red')
plt.plot(SS0temp,SS0pH,color='green')
plt.plot(Eocenetemp,EocenepH,color='blue')

plt.ylabel('pH',fontsize=14)
plt.xlabel('Temperature (K)',fontsize=14)
plt.title('Ocean pH for all Model Runs',fontsize=20)
plt.legend()
plt.show()

