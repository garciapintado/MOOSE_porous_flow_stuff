# generation of water thermodynamic tables via CoolProp [IAPWS 1995 formulation]
# script sent by Chris Green, CSIRO at 19-05-2022
# python3 water95.py

import CoolProp.CoolProp as cp
import pandas as pd
import numpy as np

# Range of pressure and temperature
#pressure = np.arange(1e7, 1.8e8+1, 5e6)      by Chris Green
#temperature = np.arange(273, 1300+1, 50)      "   "
atm2pascal = 101325     # [Pa/atm]
mH2O2pascal = 9806.65   # [Pa/mH2O]
Pmax = 45000*mH2O2pascal  # don't expect to ever go beyond 45 km mH2O depth 

#pressure = np.r_[9e5, np.arange(1e6, Pmax, 5e6)]  # Pa
pressure = np.arange(1e6, Pmax, 5e6)  # Pa
temperature = np.r_[273.17,np.arange(300, 1300+1, 50)]  # K

data = np.zeros((pressure.size * temperature.size, 6))

for i in np.arange(pressure.size):
    for j in np.arange(temperature.size):
        data[i * temperature.size + j, 0] = pressure[i]
        data[i * temperature.size + j, 1] = temperature[j]
        data[i * temperature.size + j, 2] = cp.PropsSI('D', 'P', pressure[i], 'T', temperature[j], 'Water') # [kg.m-3] mass density
        data[i * temperature.size + j, 3] = cp.PropsSI('V', 'P', pressure[i], 'T', temperature[j], 'Water') # [Pa.s] dynamic viscosity
        data[i * temperature.size + j, 4] = cp.PropsSI('H', 'P', pressure[i], 'T', temperature[j], 'Water') # [J.kg-1] mass specific enthalpy
        data[i * temperature.size + j, 5] = cp.PropsSI('U', 'P', pressure[i], 'T', temperature[j], 'Water') # [J.kg-1] mass specific internal energy
        # 'C' -> 'cp' mass specific constant pressure specific heat [J.kg-1.K-1]
        # 'Cvmass' -> 'cv' mass specific constant volume specific heat [J.kg-1.K-1]
        # 'L' -> thermal conductivity [W.m-1.K-1]
        # 'S' -> mass specific entropy [J.kg-1]

df = pd.DataFrame(data, columns = ['pressure', 'temperature', 'density', 'viscosity', 'enthalpy', 'internal_energy'])
df.to_csv('water_IAPWS95.csv', index=False)
