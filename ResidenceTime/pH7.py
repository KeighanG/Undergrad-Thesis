import PyCO2SYS as pyco2
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

pH=7
#------------------------------------------------------------------------------#
CO21=280

results1 = pyco2.sys(par1=pH, par2=CO21, par1_type=3, par2_type=4)
DIC1=results1['dic']
Alk1=results1['alkalinity']

#------------------------------------------------------------------------------#
CO22=560

results2 = pyco2.sys(par1=pH, par2=CO22, par1_type=3, par2_type=4)
DIC2=results2['dic']
Alk2=results2['alkalinity']
#------------------------------------------------------------------------------#
CO24=1120

results4 = pyco2.sys(par1=pH, par2=CO24, par1_type=3, par2_type=4)
DIC4=results4['dic']
Alk4=results4['alkalinity']
print(DIC4)

#------------------------------------------------------------------------------#
#Calculate reservoirs:
#Atmosphere
ConvCO2=2.124 #Gt of C/ppm from Friedlingstein (2019)
AtmC1=CO21*ConvCO2
AtmC2=CO22*ConvCO2
AtmC4=CO24*ConvCO2

#Ocean
rhosw=1030 #kg/m^3
conv1=10**6 #mol/umol
DICmm=12.01
conv2=10**9 #m^3/km^3
conv3=10**15 #Gt/g
Vocean=1.3324*10**9 #km^3 from Charlette and Smith 2010

OceanC1=((DIC1*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC2=((DIC2*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC4=((DIC4*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean

#Totals
TotalC1=OceanC1+AtmC1
TotalC2=OceanC2+AtmC2
TotalC4=OceanC4+AtmC4
print(TotalC1)
print(TotalC2)
print(TotalC4)

#Fluxes from weathering
Flux1=39.535191022627714*10**12*12.011/10**15 #GtC/yr
Flux2=67.62247528336282*10**12*12.011/10**15 #GtC/yr
Flux4=209.96284461413182*10**12*12.011/10**15 #GtC/yr

#Residence times
Residence1=(TotalC1/Flux1)*10**-5
Residence2=(TotalC2/Flux2)*10**-5
Residence4=(TotalC4/Flux4)*10**-5

print(Residence1)
print(Residence2)
print(Residence4)

#Plotting the Data
Conc=('280','560','1120')
Results=[Residence1,Residence2,Residence4]

plt.plot(Conc,Results,'--ko')
plt.ylabel('Residence Time (10^5 yr)',fontsize=14)
plt.xlabel('Atmospheric $CO_2$ Concentration (ppm)',fontsize=14)
plt.title('CO2 Ocean-Atmosphere Residence Times for Different Atmospheric $CO_2$ Concentrations',fontsize=18)
plt.ylim(bottom=0)
plt.show()

