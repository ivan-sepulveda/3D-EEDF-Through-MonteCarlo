import numpy as np
import scipy
from scipy.interpolate import interp1d #to interpolate data
from scipy import integrate  #to integrate the data
from math import cos, exp, pi
import pandas as pd
import glob
import os


#read EEDF function
data = pd.read_excel('EnergyDistribution.xlsx',header=None)
energy = data[0]   #energies
eedf = data[1]		#eedfs
print(data)
#energy[len(data)] = 1000.0
#eedf[len(data)] = 0.0
print(sum(eedf)*0.5, len(data))

y = interp1d(energy, energy)#*eedf/(sum(eedf)*0.5), kind='cubic')

#function we want to integrate
def f(x):
	#return np.cos(-2 * x * pi) + 3.2
	return y(x)

#integrate.quad(f, 0, np.inf)
integral,err = integrate.quad(f, energy[0], energy[len(data)-1])
print(integral*2/3, err)
simpson = integrate.simps(2*energy*eedf/(3*sum(eedf)*0.5), energy)
print(simpson)
romberg = integrate.romberg(f, energy[0], energy[len(data)-1])
print(romberg*2/3)
