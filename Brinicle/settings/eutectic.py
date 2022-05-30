#Graph of fusion point at given temperature and salinity
import numpy as np
import itertools
import matplotlib.pyplot as plt

with open('settings/parameters.txt','r') as file:
    for l in file.readlines():
        if 'Initial_temperature' in l:
            Initial_temperature = float(l.split(' ',1)[0])
        if 'Inflow_temperature' in l:
            Inflow_temperature = float(l.split(' ',1)[0])
        if 'Initial_salinity' in l:
            Initial_salinity = float(l.split(' ',1)[0])
        if 'Inflow_salinity' in l:
            Inflow_salinity = float(l.split(' ',1)[0])

with open('code/header.h', 'r') as header:
    for l in header.readlines():
        if 'FusionPoint_a' in l:
            FusionPoint_a = float((l.split(' ')[-1]).split(';')[0])
        if 'FusionPoint_b' in l:
            FusionPoint_b = float((l.split(' ')[-1]).split(';')[0])
            break

T_ref = Inflow_temperature - Initial_temperature
T0_ref = Initial_temperature
S_ref = Inflow_salinity - Initial_salinity
S0_ref = Initial_salinity

ZeroTemperature = T0_ref/T_ref
ZeroSalinity = S0_ref/S_ref

FusionPoint_a *= S_ref/T_ref
FusionPoint_b *= np.power(S_ref,3)/T_ref

T, S = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))

T_ = T + ZeroTemperature
S_ = S + ZeroSalinity
fp_function = (1-2*np.signbit(T_ref))*(T_ - (FusionPoint_a*S_ + FusionPoint_b*np.power(S_, 3)))

FusionPoint= np.where(fp_function<0, 0, fp_function)

max = np.amax(FusionPoint)
fig, ax = plt.subplots(layout='constrained')
pc = ax.pcolormesh(T, S, FusionPoint, vmin=-max, vmax=max, cmap='RdBu_r')
fig.colorbar(pc, ax=ax)
ax.set_title('Eutectic Curve')
ax.set_xlabel('Temperature')
ax.set_ylabel('Salinity')
plt.savefig('results/eutectic_curve.png')
