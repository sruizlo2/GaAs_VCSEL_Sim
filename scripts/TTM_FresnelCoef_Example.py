import numpy as np
from matplotlib import pyplot as plt
from TTM import Transfer_matrix_system as TM
from TTM import Trans_and_reflec as TR
from TTM import Transmitted_Angle as TA

# Example usage
n_air = 1;
n_glass = 1.5
layers = [
    (n_air, 0),         # Layer 0: n=1.0, air
    (n_glass, 0)
]

angles_inc = np.linspace(0, np.pi / 2, 100)
system_TE_trans = np.empty_like(angles_inc)
system_TE_reflect = np.empty_like(angles_inc)
system_TM_trans = np.empty_like(angles_inc)
system_TM_reflect =np.empty_like(angles_inc)

wavelength = 850e-9  # 500 nm
for i in range(len(angles_inc)):
    transfer_matrix_TE_system, transfer_matrix_TM_system = TM(layers, wavelength, angles_inc[i])
    system_TE_trans_, system_TE_reflect_ = TR(transfer_matrix_TE_system)
    system_TM_trans_, system_TM_reflect_ = TR(transfer_matrix_TM_system)
    system_TE_trans[i] = system_TE_trans_
    system_TE_reflect[i] = system_TE_reflect_
    system_TM_trans[i] = system_TM_trans_
    system_TM_reflect[i] = system_TM_reflect_

# Angle of transmitted light, using Snell's law
angles_trans = TA(n_air, n_glass, angles_inc)

# Factor to go from |t|^2 to |T|^2
t_factor = n_glass / n_air * np.cos(angles_trans) / np.cos(angles_inc)

# Create the plot
plt.plot(np.rad2deg(angles_inc), t_factor * system_TE_trans, label='T_TE')
plt.plot(np.rad2deg(angles_inc), t_factor * system_TM_trans, label='T_TM')
plt.plot(np.rad2deg(angles_inc), system_TE_reflect, label='R_TE')
plt.plot(np.rad2deg(angles_inc), system_TM_reflect, label='R_TM')

# Add labels and title
plt.title('Fresnel Coefficients')
# Display the legend
plt.legend()
# Set the x and y axis limits
plt.xlim(0, 90)
plt.ylim(0, 1)

# Show the plot
plt.show()
# Compare to https://en.wikipedia.org/wiki/Fresnel_equations#/media/File:Fresnel_power_air-to-glass.svg

# GaAs, index of refraction @ 850 nm = 3.569
# AlGaAs x = 45% index of refraction @ 850 nm = 3.299
