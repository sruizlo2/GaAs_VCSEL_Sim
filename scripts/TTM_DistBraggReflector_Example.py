# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 21:05:11 2024

@author: user
"""

import numpy as np
from matplotlib import pyplot as plt
from TTM import Transfer_matrix_system as TM
from TTM import Trans_and_reflec as TR
from TTM import Transmitted_Angle as TA

# Example usage
n_pairs = 15; # Number of pairs of GaAs|AlGaAs layers
#index_1 = 3.569 # For GaAs
#index_2 = 3.299 # For AlGaAs
delta_index = 0.55 #index_1 - index_2 # Delta n
average_index = 3.3 #(index_2 + index_1) / 2 # Average n
index_1 = 3.3 - 0.55/2 # For GaAs
index_2 = 3.3 + 0.55/2 # For AlGaAs
wavelength = 875e-9 # Nanometers
thickness = wavelength / 2 / average_index # lambda / 4 condition (per layer)

# Enter from GaAs
layers = [(index_2, wavelength)]
#layers = []

for i in range(n_pairs):
    layer1 = (index_1, thickness / 2)
    layer2 = (index_2, thickness / 2)
    layers.append(layer1)
    layers.append(layer2)

# Exit to air
layers.append((1, 0))

angle_inc = 0
transfer_matrix_TE_system, transfer_matrix_TM_system = TM(layers, wavelength, angle_inc)
system_TE_trans, system_TE_reflect = TR(transfer_matrix_TE_system)
system_TM_trans, system_TM_reflect = TR(transfer_matrix_TM_system)

print("Total Transfer Matrix of the system:")
print(transfer_matrix_TE_system)
print(f"Transmittance TE: {system_TE_trans}, Reflectance TE:  {system_TE_reflect}.")
print(f"Transmittance TM: {system_TM_trans}, Reflectance TM:  {system_TM_reflect}.")

DBG_reflectance = np.tanh(abs(2j * delta_index / wavelength * n_pairs * thickness)) ** 2
print(DBG_reflectance)
# GaAs, index of refraction @ 850 nm = 3.569
# AlGaAs x = 45% index of refraction @ 850 nm = 3.299
