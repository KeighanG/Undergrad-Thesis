import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
#------------------------------------------------------------------------------#
#Residence Times at fixed pH
pH91=8.080863375582929*10**5
pH92=9.448884437982963*10**5
pH94=6.086381193178489*10**5


pH81=0.4458577538672561*10**5
pH82=0.5213376586467101*10**5
pH84=0.335813157808118*10**5


pH71=0.05457326359206179*10**5
pH72=0.06381205041073683*10**5
pH74=0.04110373727895586*10**5


#Residence time using adaptive pH
Halevy1=0.6325327438124191*10**5
Halevy2=0.5002042944049496*10**5
Halevy4=0.22030308529382642*10**5


#Globally Averaged Temperatures
#CESM
TCMIP1=287.1968773085694
TCMIP2=289.57790268053
TCMIP4=298.2636710030565

pH9=np.array([pH91,pH92,pH94])
pH8=np.array([pH81,pH82,pH84])
pH7=np.array([pH71,pH72,pH74])
HalevypH=np.array([Halevy1,Halevy2,Halevy4])

CMIPtemp=np.array([TCMIP1,TCMIP2,TCMIP4])

plt.scatter(CMIPtemp,pH9,color='red',label='pH=9.0')
plt.scatter(CMIPtemp,pH8,color='green',label='pH=8.0')
plt.scatter(CMIPtemp,pH7,color='blue',label='pH=7.0')
plt.scatter(CMIPtemp,HalevypH,color='black',label='Halevy pH')

plt.plot(CMIPtemp,pH9,color='red')
plt.plot(CMIPtemp,pH8,color='green')
plt.plot(CMIPtemp,pH7,color='blue')
plt.plot(CMIPtemp,HalevypH,color='black')

plt.ylabel('Residence Time',fontsize=14)
plt.yscale('log')
plt.xlabel('Globally Averaged Temperature (K)',fontsize=14)
#plt.title('Residence Times for all Model Runs',fontsize=20)
plt.legend(loc='center right',fontsize=14)
plt.show()

