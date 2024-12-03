import pandas as pd
import os
from matplotlib import pyplot as plt

path = 'C:\\Users\\petee\\Documents\\GitHub\\GaAs_VCSEL_Sim\\VCSEL Data\\3QW GaAs 3V sweep(high reflection).dat'
data = pd.read_csv(path, skiprows=1, header=None)

voltage = data[0]
current = data[4]
Leftpower = data[1]

plt.plot(current,Leftpower)
plt.xlabel('Current (A/cm2)')
plt.ylabel('output power(mW)')
plt.yscale('log')
plt.xscale('log')
plt.xlim(10,10000)
plt.ylim(0.0000001,10)

#%%
import math
#values for 0V
Ec = -.1
Efn = 0
Ev = -.02
Efp = 0


fc = 1/(1+math.exp((Ec-Efn)/(300*8.617E-5)))
fv = 1/(1+math.exp((Ev-Efp)/(300*8.617E-5)))

alphr = -math.log(0.995*0.95)/(2*0.241E-4)

gain = 1200*3*(fc-fv)

print(f'fc:{fc}')
print(f'fv:{fv}')
print(f'gain:{gain}')
print(f'mirror loss:{alphr}')