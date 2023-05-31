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
#S/So=1.0
#Read in Data
path1='/home/czg/FYSP_clouds_archive/CAM4/cam4_10.nc'
d1=Dataset(path1,'r')
co2vmr1=d1.variables['co2vmr'][:]
co2vmr1=co2vmr1[0]
CO21=co2vmr1*10**6

#Convert ppm to pCO2
pCO21=co2vmr1*P

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

fCO21=pCO21
#Apply Henry's law to obtain [CO2(aq)]
CO2aq1=KHCO2*fCO21 #mol/kg
CO2aq1=CO2aq1*10**6 #umol/kg

results1 = pyco2.sys(par1=pH1, par2=CO21, par1_type=3, par2_type=4)
DIC1=results1['dic']
Alk1=results1['alkalinity']

#------------------------------------------------------------------------------#
#S/So=1.1
#Read in Data
path11='/home/czg/FYSP_clouds_archive/CAM4/cam4_11.nc'
d11=Dataset(path11,'r')
co2vmr11=d11.variables['co2vmr'][:]
co2vmr11=co2vmr11[0]
CO211=co2vmr11*10**6

#Convert ppm to pCO2
pCO211=co2vmr11*P

#Coefficients
Ca11=(2*nsat*Ksp)/(pCO211*KH*KA1*KA2)
Cb11=1
Cc11=Kchnonca
Cd11=-(pCO211*KH*KA1+Kw)
Ce11=-(2*pCO211*KH*KA1*KA2)

Coeffs11=[Ca11, Cb11, Cc11, Cd11, Ce11]
H11=np.roots(Coeffs11)
#Only first root is positive and physically meaningful for a H+ concentration
pH11=-np.log10(H11[0])

fCO211=pCO211
#Apply Henry's law to obtain [CO2(aq)]
CO2aq11=KHCO2*fCO211 #mol/kg
CO2aq11=CO2aq11*10**6 #umol/kg

results11 = pyco2.sys(par1=pH11, par2=CO211, par1_type=3, par2_type=4)
DIC11=results11['dic']
Alk11=results11['alkalinity']

#------------------------------------------------------------------------------#
#S/So=0.9
#Read in Data
path9='/home/czg/FYSP_clouds_archive/CAM4/cam4_09.nc'
d9=Dataset(path9,'r')
co2vmr9=d9.variables['co2vmr'][:]
co2vmr9=co2vmr9[0]
CO29=co2vmr9*10**6

#Convert ppm to pCO2
pCO29=co2vmr9*P

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

fCO29=pCO29
#Apply Henry's law to obtain [CO2(aq)]
CO2aq9=KHCO2*fCO29 #mol/kg
CO2aq9=CO2aq9*10**6 #umol/kg

results9 = pyco2.sys(par1=pH9, par2=CO29, par1_type=3, par2_type=4)
DIC9=results9['dic']
Alk9=results9['alkalinity']

#------------------------------------------------------------------------------#
#S/So=0.8
#Read in Data
path8='/home/czg/FYSP_clouds_archive/CAM4/cam4_08.nc'
d8=Dataset(path8,'r')
co2vmr8=d8.variables['co2vmr'][:]
co2vmr8=co2vmr8[0]
CO28=co2vmr8*10**6

#Convert ppm to pCO2
pCO28=co2vmr8*P

#Coefficients
Ca8=(2*nsat*Ksp)/(pCO28*KH*KA1*KA2)
Cb8=1
Cc8=Kchnonca
Cd8=-(pCO28*KH*KA1+Kw)
Ce8=-(2*pCO28*KH*KA1*KA2)

Coeffs8=[Ca8, Cb8, Cc8, Cd8, Ce8]
H8=np.roots(Coeffs8)
#Only first root is positive and physically meaningful for a H+ concentration
pH8=-np.log10(H8[0])

fCO28=pCO28
#Apply Henry's law to obtain [CO2(aq)]
CO2aq8=KHCO2*fCO28 #mol/kg
CO2aq8=CO2aq8*10**6 #umol/kg

results8 = pyco2.sys(par1=pH8, par2=CO28, par1_type=3, par2_type=4)
DIC8=results8['dic']
Alk8=results8['alkalinity']

print(Alk8)
print(Alk9)
print(Alk1)
print(Alk11)
#------------------------------------------------------------------------------#
#Calculate reservoirs:
#Atmosphere
ConvCO2=2.124 #Gt of C/ppm from Friedlingstein (2019)
AtmC1=CO21*ConvCO2
AtmC11=CO211*ConvCO2
AtmC9=CO29*ConvCO2
AtmC8=CO28*ConvCO2

#Ocean
rhosw=1030 #kg/m^3
conv1=10**6 #mol/umol
DICmm=12.01
conv2=10**9 #m^3/km^3
conv3=10**15 #Gt/g
Vocean=1.3324*10**9 #km^3 from Charlette and Smith 2010

OceanC1=((DIC1*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC11=((DIC11*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC9=((DIC9*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean
OceanC8=((DIC8*rhosw)/conv1)*((DICmm*conv2)/conv3)*Vocean

#Totals
TotalC1=OceanC1+AtmC1
TotalC11=OceanC11+AtmC11
TotalC9=OceanC9+AtmC9
TotalC8=OceanC8+AtmC8

#Fluxes from weathering
Flux8=430.1627701316608*10**12*12.011/10**15 #GtC/yr
Flux9=199.29218259945233*10**12*12.011/10**15 #GtC/yr
Flux1=57.30944275863661*10**12*12.011/10**15  #GtC/yr
Flux11=8.487303642047904*10**12*12.011/10**15  #GtC/yr

#Residence times
Residence1=(TotalC1/Flux1)*10**-5
Residence11=(TotalC11/Flux11)*10**-5
Residence9=(TotalC9/Flux9)*10**-5
Residence8=(TotalC8/Flux8)*10**-5

print(Residence8)
print(Residence9)
print(Residence1)
print(Residence11)

#Plotting the Data
Solars=('0.8','0.9','1.0','1.1')
Results=[Residence8,Residence9,Residence1,Residence11]

plt.plot(Solars,Results,'--ko')
plt.ylabel('Residence Time (10^5 yr)',fontsize=14)
plt.xlabel('Solar Constant',fontsize=14)
plt.title('CO2 Ocean-Atmosphere Residence Times for Different Solar Constants',fontsize=18)

plt.show()

