import PyCO2SYS as pyco2
import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

#Defining constants needed in calculations
#Equation: 0=((2*nsat*Ksp*(H**4)))/(pCO2*KH*KA1*KA2)+(H**3)+(Kchnonca*(H**2))-((pCO2*KH*KA1+Kw)*H)-(2*pCO2*KH*KA1*KA2)
nsat=3
Ksp=10**-6.37
KH=10**-1.28
KA1=10**-6.06
KA2=10**-9.29
Kchnonca=-18.4*10**-3
Kw=10**-14.08
KHCO2=0.0324 #mol/kg/atm
P=1 #atm, assuming surface pressure
#------------------------------------------------------------------------------#
CO21=280
pCO21=CO21*10**-6*P #Assuming surface pressure

#Coefficients
Ca1=(2*nsat*Ksp)/(pCO21*KH*KA1*KA2)
Cb1=1
Cc1=Kchnonca
Cd1=-(pCO21*KH*KA1+Kw)
Ce1=-(2*pCO21*KH*KA1*KA2)

Coeffs1=[Ca1, Cb1, Cc1, Cd1, Ce1]
H1=np.roots(Coeffs1)
#Only first root is positive and physically meaningful for a H+ concentration
pH1=-np.log10(H1[0])

results1 = pyco2.sys(par1=pH1, par2=CO21, par1_type=3, par2_type=4)
DIC1=results1['dic']
Alk1=results1['alkalinity']

#------------------------------------------------------------------------------#
CO22=560
pCO22=CO22*10**-6*P
#Coefficients
Ca2=(2*nsat*Ksp)/(pCO22*KH*KA1*KA2)
Cb2=1
Cc2=Kchnonca
Cd2=-(pCO22*KH*KA1+Kw)
Ce2=-(2*pCO22*KH*KA1*KA2)

Coeffs2=[Ca2, Cb2, Cc2, Cd2, Ce2]
H2=np.roots(Coeffs2)
#Only first root is positive and physically meaningful for a H+ concentration
pH2=-np.log10(H2[0])

results2 = pyco2.sys(par1=pH2, par2=CO22, par1_type=3, par2_type=4)
DIC2=results2['dic']
Alk2=results2['alkalinity']

#------------------------------------------------------------------------------#
CO24=1120
pCO24=CO24*10**-6*P
#Coefficients
Ca4=(2*nsat*Ksp)/(pCO24*KH*KA1*KA2)
Cb4=1
Cc4=Kchnonca
Cd4=-(pCO24*KH*KA1+Kw)
Ce4=-(2*pCO24*KH*KA1*KA2)

Coeffs4=[Ca4, Cb4, Cc4, Cd4, Ce4]
H4=np.roots(Coeffs4)
#Only first root is positive and physically meaningful for a H+ concentration
pH4=-np.log10(H4[0])

results4 = pyco2.sys(par1=pH4, par2=CO24, par1_type=3, par2_type=4)
DIC4=results4['dic']
Alk4=results4['alkalinity']

print(Alk1)
print(Alk2)
print(Alk4)
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
Flux1=39.5221934361353*10**12*12.011/10**15 #GtC/yr
Flux2=67.59046715605632*10**12*12.011/10**15 #GtC/yr
Flux4=209.85702197258075*10**12*12.011/10**15 #GtC/yr

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

