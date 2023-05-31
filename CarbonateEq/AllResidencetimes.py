import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
#------------------------------------------------------------------------------#
CESM1=0.6325327438124191
CESM2=0.5002042944049496
CESM4=0.22030308529382642

S08=0.5748054454451945
S09=0.48301796670380176
S01=0.4911511305258789
S011=0.7320254203628456

Eocene1=0.5471566121238528
Eocene3=0.2797069107108757
Eocene6=0.1582603629993568
Eocene9=0.08042237215346709

#Globally Averaged Temperatures
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


Times=np.array([CESM1,CESM2,CESM4,S08,S09,S01,S011,Eocene1,Eocene3,Eocene6,Eocene9])
Temps=np.array([TCESM1,TCESM2,TCESM4,TS08,TS09,TS01,TS011,TE1,TE3,TE6,TE9])

CESMtime=np.array([CESM1,CESM2,CESM4])
SS0time=np.array([S08,S09,S01,S011])
Eocenetime=np.array([Eocene1,Eocene3,Eocene6,Eocene9])

CESMtemp=np.array([TCESM1,TCESM2,TCESM4])
SS0temp=np.array([TS08,TS09,TS01,TS011])
Eocenetemp=np.array([TE1,TE3,TE6,TE9])

CESMpH=np.array([pHCESM1,pHCESM2,pHCESM4])
SS0pH=np.array([pHS08,pHS09,pHS01,pHS011])
EocenepH=np.array([pHE1,pHE3,pHE6,pHE9])

plt.scatter(CESMtemp,CESMtime,color='red',label='CMIP-6')
plt.scatter(SS0temp,SS0time,color='green',label='$S/S_0$')
plt.scatter(Eocenetemp,Eocenetime,color='blue',label='Eocene')

plt.plot(CESMtemp,CESMtime,color='red')
plt.plot(SS0temp,SS0time,color='green')
plt.plot(Eocenetemp,Eocenetime,color='blue')

plt.ylabel('Residence Time ($10^5$ yr)',fontsize=14)
plt.xlabel('Temperature (K)',fontsize=14)
#plt.title('Residence Times for all Model Runs',fontsize=20)
plt.legend()
plt.show()

