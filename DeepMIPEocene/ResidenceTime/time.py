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
CO23=840
pCO23=CO23*10**-6*P
#Coefficients
Ca3=(2*nsat*Ksp)/(pCO23*KH*KA1*KA2)
Cb3=1
Cc3=Kchnonca
Cd3=-(pCO23*KH*KA1+Kw)
Ce3=-(2*pCO23*KH*KA1*KA2)

Coeffs3=[Ca3, Cb3, Cc3, Cd3, Ce3]
H3=np.roots(Coeffs3)
#Only first root is positive and physically meaningful for a H+ concentration
pH3=-np.log10(H3[0])

results3 = pyco2.sys(par1=pH3, par2=CO23, par1_type=3, par2_type=4)
DIC3=results3['dic']
Alk3=results3['alkalinity']

#------------------------------------------------------------------------------#
CO26=1680
pCO26=CO26*10**-6*P
#Coefficients
Ca6=(2*nsat*Ksp)/(pCO26*KH*KA1*KA2)
Cb6=1
Cc6=Kchnonca
Cd6=-(pCO26*KH*KA1+Kw)
Ce6=-(2*pCO26*KH*KA1*KA2)

Coeffs6=[Ca6, Cb6, Cc6, Cd6, Ce6]
H6=np.roots(Coeffs6)
#Only first root is positive and physically meaningful for a H+ concentration
pH6=-np.log10(H6[0])

results6 = pyco2.sys(par1=pH6, par2=CO26, par1_type=3, par2_type=4)
DIC6=results6['dic']
Alk6=results6['alkalinity']

#------------------------------------------------------------------------------#
CO29=2520
pCO29=CO29*10**-6*P
#Coefficients
Ca9=(2*nsat*Ksp)/(pCO29*KH*KA1*KA2)
Cb9=1
Cc9=Kchnonca
Cd9=-(pCO29*KH*KA1+Kw)
Ce9=-(2*pCO29*KH*KA1*KA2)

Coeffs9=[Ca9, Cb9, Cc9, Cd9, Ce9]
H9=np.roots(Coeffs9)
#Only first root is positive and physically meaningful for a H+ concentration
pH9=-np.log10(H9[0])

results9 = pyco2.sys(par1=pH9, par2=CO29, par1_type=3, par2_type=4)
DIC9=results9['dic']
Alk9=results9['alkalinity']

print(pH1)
print(pH3)
print(pH6)
print(pH9)
#------------------------------------------------------------------------------#
#Calculate reservoirs:
#Atmosphere
ConvCO2=2.124 #Gt of C/ppm from Friedlingstein (2019)
AtmC1=CO21*ConvCO2
AtmC3=CO23*ConvCO2
AtmC6=CO26*ConvCO2
AtmC9=CO29*ConvCO2

#Ocean
rhosw=1030 #kg/m^3
conv1=10**6 #mol/umol
DICmm=12.01
conv2=10**9 #m^3/km^3
conv3=10**15 #Gt/g
Vocean=1.3324*10**9 #km^3 from Charlette and Smith 2010

OceanC1=((DIC1*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC3=((DIC3*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC6=((DIC6*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC9=((DIC9*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean

#Totals
TotalC1=OceanC1+AtmC1
TotalC3=OceanC3+AtmC3
TotalC6=OceanC6+AtmC6
TotalC9=OceanC9+AtmC9

print(TotalC1)
print(TotalC3)
print(TotalC6)
print(TotalC9)

#Fluxes from weathering
Flux1=45.68907859599276*10**12*12.011/10**15 #GtC/yr
Flux3=144.9780144946899*10**12*12.011/10**15 #GtC/yr
Flux6=352.4012370573544*10**12*12.011/10**15 #GtC/yr
Flux9=839.3008102350059*10**12*12.011/10**15 #GtC/yr

#Residence times
Residence1=(TotalC1/Flux1)*10**-5
Residence3=(TotalC3/Flux3)*10**-5
Residence6=(TotalC6/Flux6)*10**-5
Residence9=(TotalC9/Flux9)*10**-5

print(Residence1)
print(Residence3)
print(Residence6)
print(Residence9)

#Plotting the Data
Conc=('280','840','1680','2520')
Results=[Residence1,Residence3,Residence6,Residence9]

plt.plot(Conc,Results,'--ko')
plt.ylabel('Residence Time (10^5 yr)',fontsize=14)
plt.xlabel('Atmospheric $CO_2$ Concentration (ppm)',fontsize=14)
plt.title('CO2 Ocean-Atmosphere Residence Times for Different Atmospheric $CO_2$ Concentrations',fontsize=18)
plt.ylim(bottom=0)
plt.show()

