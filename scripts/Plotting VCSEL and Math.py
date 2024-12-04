import pandas as pd
import os
from matplotlib import pyplot as plt

user = 'Srulo'
if user == 'Pete':
    path = 'C:\\Users\\petee\\Documents\\GitHub\\GaAs_VCSEL_Sim\\VCSEL Data\\'
else:
    path = 'I:\\Documentos\GaAs_VCSEL_Sim\\'
    
data = pd.read_csv(path + '3QW GaAs 3V sweep(high reflection).dat', skiprows=1, header=None)

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
Ec = -5.06
Efn = -5.0
Ev = -7.03
Efp = -7.0


fc = 1/(1+math.exp((Ec-Efn)/(300*8.617E-5)))
fv = 1/(1+math.exp((Ev-Efp)/(300*8.617E-5)))


alphr = -math.log(0.995*0.99)/(2*1.04E-4)

gain = 5760*(fc-fv)

print(fc-fv)
# print(f'fc:{fc}')
# print(f'fv:{fv}')
# print(f'gain:{gain}')
# print(f'mirror loss:{alphr}')

#%%
#no bias band diagram
import pandas as pd
import os
from matplotlib import pyplot as plt
data = pd.read_csv(path + 'Band Diagram 0 Bias.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 6, 3
Position = data[0]
Ec = data[1]
Efn = data[2]
Ev = data[3]
Efp = data[4]
plt.plot(Position,Ec,color='red')
plt.plot(Position,Efn,color='blue')
plt.plot(Position,Ev,color='maroon')
plt.plot(Position,Efp,color='navy')
plt.xlabel('Position')
plt.ylabel('Energy (eV)')
plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
plt.xlim(1.7,2.4)
# plt.ylim(0.0000001,10)


#%%
#biased band diagram
import pandas as pd
import os
from matplotlib import pyplot as plt
data = pd.read_csv(path + 'Band Diagram T Bias.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 6, 3
Position = data[0]
Ec = data[1]
Efn = data[2]
Ev = data[3]
Efp = data[4]
plt.plot(Position,Ec,color='red')
plt.plot(Position,Efn,color='blue')
plt.plot(Position,Ev,color='maroon')
plt.plot(Position,Efp,color='navy')
plt.xlabel('Position')
plt.ylabel('Energy (eV)')
plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
plt.xlim(1.7,2.4)
# plt.ylim(0.0000001,10)

print(Ec)

#%%
#Output power from stimulated recombination old 1 QW
Ustim = 2.1e24 #1/(cm3 s)
A = 38.484e-8
T = 8e-7
V = A*T #in cm3
E850 = 3e8*6.626e-34/(850e-9) #in J
power = E850*V*Ustim*1000 #in mW
print(power)

#%%

freq850nm = 3.53e14
h = 1.055e-34
elec = 1.602*10e-19 
A = 38.484e-8

power=freq850nm*h/elec*0.52*1000*A
print(power)

#%%
#Output power from stimulated recombination 5 QW
Ustim = 2.4e25 #1/(cm3 s)
A = 38.484e-8
T = 8e-7
V = A*T #in cm3
E850 = 3e8*6.626e-34/(850e-9) #in J
power = E850*V*Ustim*1000 #in mW
print(power)

#%%
#biased band diagram
import pandas as pd
import os
from matplotlib import pyplot as plt
data = pd.read_csv(path + 'New Device Band Diagram T Bias.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 6, 3
Position = data[0]
Ec = data[1]
Efn = data[2]
Ev = data[3]
Efp = data[4]
plt.plot(Position,Ec,color='red')
plt.plot(Position,Efn,color='blue')
plt.plot(Position,Ev,color='maroon')
plt.plot(Position,Efp,color='navy')
plt.xlabel('Position')
plt.ylabel('Energy (eV)')
plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
plt.xlim(1.7,2.4)
# plt.ylim(0.0000001,10)

print(Ec)

#%%
data = pd.read_csv(path + 'VCSEL_Mod1_Recombination.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 6, 3
Position = data[0]
total = data[1]
stim = data[2]
plt.plot(Position,total,color='red')
plt.plot(Position,stim,color='blue')
plt.xlabel('Position', fontsize=14)
plt.ylabel('Recombination (1/cm^3s)', fontsize=14)
plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
plt.xlim(3.2,3.5)
#%%
data = pd.read_csv(path + 'VCSEL_Opt_Recombination.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 6, 3
Position = data[0]
total = data[1]
stim = data[2]
plt.plot(Position,total,color='red')
plt.plot(Position,stim,color='blue')
plt.xlabel('Position', fontsize=14)
plt.ylabel('Recombination (1/cm^3s)', fontsize=14)
#plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
plt.xlim(3.1,3.3)

#%%
data = pd.read_csv(path + 'VCSEL_Mod1_MACRO_Output.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 4, 4
Voltage = data[0]
Power = data[1]
Current = data[3]
plt.loglog(Current,Power,color='red')
plt.xlabel('Current density [A/cm^2]', fontsize=14)
plt.ylabel('Optical power (mW)', fontsize=14)
plt.grid()
#%%
data = pd.read_csv(path + 'VCSEL_Opt_MACRO_Output.dat', skiprows=1, header=None)

plt.rcParams['figure.figsize'] = 4, 4
Voltage = data[0]
Power = data[1]
Current = data[3]
plt.loglog(Current,Power,color='blue')
plt.xlabel('Current density [A/cm^2]', fontsize=14)
plt.ylabel('Optical power (mW)', fontsize=14)
plt.grid()
