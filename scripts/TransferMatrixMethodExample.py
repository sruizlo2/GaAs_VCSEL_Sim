from TTM import transfer_matrix_system as TM
from TTM import trans_and_reflec as TR

# Example usage
layers = [
    (1, 0),         # Layer 0: n=1.0, air
    (1.5, 200e-9),  # Layer 1: n=1.5, d=200nm
    (2.0, 100e-9),  # Layer 2: n=2.0, d=100nm
    (1.4, 150e-9)   # Layer 3: n=1.4, d=150nm
]

wavelength = 850e-9  # 500 nm
angle_inc = 0  # Normal incidence

transfer_matrix_TE_system, transfer_matrix_TM_system = TM(layers, wavelength, angle_inc)
system_TE_trans, system_TE_reflect = TR(transfer_matrix_TE_system)
system_TM_trans, system_TM_reflect = TR(transfer_matrix_TM_system)

print("Total Transfer Matrix of the system:")
print(transfer_matrix_TE_system)
print(f"Transmittance TE: {system_TE_trans}, Reflectance TE:  {system_TE_reflect}.")
print(f"Transmittance TM: {system_TM_trans}, Reflectance TM:  {system_TM_reflect}.")


# GaAs, index of refraction @ 850 nm = 3.569
# AlGaAs x = 45% index of refraction @ 850 nm = 3.299
