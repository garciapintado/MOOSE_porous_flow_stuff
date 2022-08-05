import matplotlib.pyplot as plt
#import pandas as pd
import numpy as np
from numpy import genfromtxt

#data = pd.read_table("water_IAPWS95_extrap.csv", delimiter=",")
data = genfromtxt('water_IAPWS95_extrap.csv', delimiter=',', names=True)
pressure = np.unique(data['pressure'])
temperature = np.unique(data['temperature'])
nP = len(pressure)
nT = len(temperature)

rhofit = data['density'].reshape(nP,nT)
plt.contourf(temperature, pressure * 1e-06, rhofit)
plt.colorbar()
plt.xlabel('Temperature(K)')
plt.ylabel('Pressure (MPa)')
plt.savefig('water_IAPWS95_extrap.png', dpi=600)
plt.close()

