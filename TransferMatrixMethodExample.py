import numpy as np
from scipy import constants

# Constants
speed_of_light = constants.c  # speed of light in vacuum (m/s)

def layer_matrix(index_of_refraction_in, index_of_refraction_out, thickness, wavelength, angle_in=0):
    """
    Compute the transfer matrix for a single layer.

    Parameters:
    - index_of_refraction_in: Refractive index of the current layer
    - index_of_refraction_out: Refractive index of the next layer
    - thickness: Thickness of the layer (in meters)
    - wavelength: Wavelength of the incident light (in meters)
    - theta_inc: Angle of incidence (in radians, default is normal incidence)

    Returns:
    - Transfer matrix of the layer (2x2 matrix)
    """
    # Calculate the wave vector in the medium (propagation constant)
    wavenumber = 2 * np.pi / wavelength  # wave number in vacuum
    wavenumber_x = wavenumber * np.sin(angle_in)   # transverse wave vector component
    wavenumber_z = np.sqrt((wavenumber * index_of_refraction_out) ** 2 - wavenumber_x ** 2)  # longitudinal wave vector component in the medium
    
    # Phase shift in the layer due to its thickness
    path_length = wavenumber_z * thickness
    # Propagation matrix
    prop_matrix = np.array([[np.exp(-1j * path_length), 0],
                            [0, np.exp(1j * path_length)]], dtype=complex)

    # Ration of indices of refraction
    indices_ratio = index_of_refraction_out / index_of_refraction_in;
    # Angle of transmitted light, using Snell's law
    angle_out = np.arcsin(np.sin(angle_in) / indices_ratio)

    # Cosines ratio
    cosines_ratio = np.cos(angle_out) / np.cos(angle_in)

    # Interface matrix (TE)
    interface_matrix_TE = 0.5 * np.array([[1.0 + (cosines_ratio * indices_ratio),
                                         1.0 - (cosines_ratio * indices_ratio)],
                                         [1.0 - (cosines_ratio * indices_ratio),
                                         1.0 + (cosines_ratio * indices_ratio)]], dtype=complex)
    # Interface matrix (TM)
    interface_matrix_TM = 0.5 * np.array([[indices_ratio + cosines_ratio,
                                         indices_ratio - cosines_ratio],
                                         [indices_ratio - cosines_ratio,
                                         indices_ratio + cosines_ratio]], dtype=complex)
    
    # Transfer_matrix (TE)
    transfer_matrix_TE = np.matmul(interface_matrix_TE, prop_matrix)
    # Transfer_matrix (TM)
    transfer_matrix_TM = np.matmul(interface_matrix_TM, prop_matrix)
    print(index_of_refraction_in)
    print(index_of_refraction_out)

    # Transfer matrix for the layer
    #transfer_matrix_TE = np.array([[np.cos(path_length), 1j / (index_of_refraction_out * np.cos(angle_in)) * np.sin(path_length)],
    #                                [1j * index_of_refraction_out * np.cos(angle_in) * np.sin(path_length), np.cos(path_length)]], dtype=complex)
    
    return transfer_matrix_TE, transfer_matrix_TM, angle_out

def transfer_matrix_system(layers, wavelength, angle_in=0):
    """
    Compute the total transfer matrix for a multilayer system.

    Parameters:
    - layers: List of tuples where each tuple contains (n, d) for each layer.
              n is the refractive index and d is the thickness of the layer.
              Example: [(n1, d1), (n2, d2), ...]
    - wavelength: Wavelength of the incident light (in meters)
    - angle_in: Angle of incidence (in radians, default is normal incidence)

    Returns:
    - Total transfer matrix of the system (2x2 matrix)
    """
    # Start with the identity matrix
    transfer_matrix_TE_system = np.eye(2, dtype=complex)
    transfer_matrix_TM_system = np.eye(2, dtype=complex)

    # Previous index of refraction
    index_of_refraction_in = 1 # Start at air

    # Multiply the transfer matrices of each layer
    for index_of_refraction_out, thickness in layers:
        # Layer
        transfer_matrix_TE_layer, transfer_matrix_TM_layer, angle_in = layer_matrix(index_of_refraction_in, index_of_refraction_out, thickness, wavelength, angle_in)
        # System
        transfer_matrix_TE_system = np.matmul(transfer_matrix_TE_system, transfer_matrix_TE_layer)
        transfer_matrix_TM_system = np.matmul(transfer_matrix_TM_system, transfer_matrix_TM_layer)
        # Index of refraction
        index_of_refraction_in = index_of_refraction_out

    return transfer_matrix_TE_system, transfer_matrix_TM_system

def trans_and_reflec(transfer_matrix):
    """
    Computes the transmittance and reflectance from a transfer matrix

    Parameters:
    - transfer_matrix: 2x2 transfer matrix

    Returns:
    - transmittance: System's transmittance 1 / abs(M_11) ^ 2
    - reflectance: System's reflectance abs(M_11 / M_21) ^ 2
    """
    transmittance = 1 / abs(transfer_matrix[0,0]) ** 2
    reflectance = abs(transfer_matrix[1,0] / transfer_matrix[0,0]) ** 2
    return transmittance, reflectance

# Example usage
layers = [
    (1, 0),         # Layer 0: n=1.0, air
    (1.5, 200e-9),  # Layer 1: n=1.5, d=200nm
    (2.0, 100e-9),  # Layer 2: n=2.0, d=100nm
    (1.4, 150e-9)   # Layer 3: n=1.4, d=150nm
]

wavelength = 850e-9  # 500 nm
angle_inc = 0  # Normal incidence

transfer_matrix_TE_system, transfer_matrix_TM_system = transfer_matrix_system(layers, wavelength, angle_inc)
system_TE_trans, system_TE_reflect = trans_and_reflec(transfer_matrix_TE_system)
system_TM_trans, system_TM_reflect = trans_and_reflec(transfer_matrix_TM_system)

print("Total Transfer Matrix of the system:")
print(transfer_matrix_TE_system)
print(f"Transmittance TE: {system_TE_trans}, Reflectance TE:  {system_TE_reflect}.")
print(f"Transmittance TM: {system_TM_trans}, Reflectance TM:  {system_TM_reflect}.")
print(system_TE_trans * n_glass + system_TE_reflect)

# GaAs, index of refraction @ 850 nm = 3.569
# AlGaAs x = 45% index of refraction @ 850 nm = 3.299
