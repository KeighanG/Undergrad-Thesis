import PyCO2SYS as pyco2
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

pyco2_kws = {} #create empty dictionary
# Define the seawater conditions and add them to the dictionary
pyco2_kws["salinity"] = 0  # practical salinity
pyco2_kws["temperature"] = 25  # lab temperature
pyco2_kws["temperature_out"] = 25  # in-situ temperature
pyco2_kws["pressure"] = 0  # lab pressure
pyco2_kws["pressure_out"] = 0  # in-situ pressure
pyco2_kws["total_silicate"]  = 0  # total silicate in μmol/kg-sw
pyco2_kws["total_phosphate"] = 0  # total phosphate in μmol/kg-sw
pyco2_kws["total_ammonia"] = 0  # total ammonia in μmol/kg-sw
pyco2_kws["total_sulfide"] = 0  # total sulfide in μmol/kg-sw

# Define PyCO2SYS settings and add them to the dictionary
pyco2_kws["opt_pH_scale"] = 1
pyco2_kws["opt_k_carbonic"] = 15
pyco2_kws["opt_k_bisulfate"] = 1
pyco2_kws["opt_total_borate"] = 1

pH_start = 0.5
pH_end = 14  # maximum pH value
n_steps = 100
pH_increasing = np.linspace(pH_start, pH_end, n_steps) # generates an array of 100 pH values from 0 to 14


pyco2_kws["par1"] = 2100
pyco2_kws["par2"] = pH_increasing  # array of increasing pH
pyco2_kws["par1_type"] = 2
pyco2_kws["par2_type"] = 3

CO2dict_response = pyco2.sys(**pyco2_kws, k_water=1e-14, k_carbonic_1=10**(-6.3), k_carbonic_2=10**(-10.3))


# offline calculations made using the example numbers in thesis
HCO3 = -10.63+pH_increasing
CO3 = -19.76+2*pH_increasing
H2CO3 = np.full((100),-4.69)
proton = np.log10(CO2dict_response["hydrogen_free"]/1e6)
OH = np.log10(CO2dict_response["OH"]/1e6)
H=CO2dict_response["hydrogen_free"]/1e6

H2=H**2

DIC=np.log10(np.full((100),10**-4.69)+((10**-10.63)/H)+((10**-19.76)/H2))


# Plotting the data, using seaborn style
sns.set_style("whitegrid")
sns.set_style("ticks", {"xtick.bottom": True, "xtick.top": False})

fig, ax = plt.subplots()
ax.plot(pH_increasing, HCO3,label='HCO$_3^-$')
ax.plot(pH_increasing, CO3,label='CO$_3^{2-}$')
ax.plot(pH_increasing, H2CO3,label='H$_2$CO$_3^*$')
ax.plot(pH_increasing, proton,label='H$^+$')
ax.plot(pH_increasing, OH,label='OH$^-$')
ax.plot(pH_increasing, DIC,'k',label='DIC')

ax.tick_params(axis='both', which='major', labelsize=10, width=1, length=7)

plt.legend(title='Species', loc='lower center',fontsize=14)
#plt.title('Log Concentration vs pH Diagram Typically Observed in the Ocean',fontsize=22)

plt.xlabel('pH', fontsize=18)
plt.ylabel('Log Concentration',fontsize=18)

plt.show()

