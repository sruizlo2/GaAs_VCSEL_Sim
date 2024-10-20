import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np

"""
TMM Functions
"""

def Transmitted_Angle(index_of_refraction_in, index_of_refraction_out, angle_in=0):
    # Ration of indices of refraction
    indices_ratio = index_of_refraction_in / index_of_refraction_out;
    # Angle of transmitted light, using Snell's law
    angle_out = np.arcsin(np.sin(angle_in) * indices_ratio)
    
    return angle_out

def Layer_matrix(index_of_refraction_in, index_of_refraction_out, thickness, wavelength, angle_in=0):
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
    angle_out = Transmitted_Angle(index_of_refraction_in, index_of_refraction_out, angle_in)
    
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

    # Transfer matrix for the layer
    #transfer_matrix_TE = np.array([[np.cos(path_length), 1j / (index_of_refraction_out * np.cos(angle_in)) * np.sin(path_length)],
    #                                [1j * index_of_refraction_out * np.cos(angle_in) * np.sin(path_length), np.cos(path_length)]], dtype=complex)
    
    return transfer_matrix_TE, transfer_matrix_TM, angle_out

def Transfer_matrix_system(layers, wavelength, angle_in=0):
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
        transfer_matrix_TE_layer, transfer_matrix_TM_layer, angle_in = Layer_matrix(index_of_refraction_in, index_of_refraction_out, thickness, wavelength, angle_in)
        # System
        transfer_matrix_TE_system = np.matmul(transfer_matrix_TE_system, transfer_matrix_TE_layer)
        transfer_matrix_TM_system = np.matmul(transfer_matrix_TM_system, transfer_matrix_TM_layer)
        # Index of refraction
        index_of_refraction_in = index_of_refraction_out

    return transfer_matrix_TE_system, transfer_matrix_TM_system

def Trans_and_reflec(transfer_matrix):
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

TM = Transfer_matrix_system
TR = Trans_and_reflec
TA = Transmitted_Angle

#loading files for pete
#path = 'C:/Users\petee\Dropbox (MIT)\classes\\3.46'
#os.chdir(path)

#Loading files for sebastian
path = '..\materials'
os.chdir(path)

#Bragg Mirror Imputs
n_pairs = 25; # Number of pairs of GaAs|AlGaAs layers
mat2 = "AlGaAs-X=0.452"
mat1 = "GaAs"

#Range of Interest
wavelength_min = 350 #in nm
wavelength_max = 1350 #in nm
wavelength_range = wavelength_max - wavelength_min
wavelength_spectra = np.linspace(wavelength_min,wavelength_max,(wavelength_range+1))

#refractive index and extinction coeff of each material
n_mat1 = []
k_mat1 = []
n_mat2 = []
k_mat2 = []

tally = 0 #i is not working in the loop for some reason so need this here to count the loops
for i in wavelength_spectra:
     Target_Wavelength = wavelength_spectra[tally]
     tally = tally + 1
     
     for filename in [
                      mat1,
                      mat2,
                      ]:
         materialinfo = pd.read_csv(f'{filename}.csv',sep=',',skiprows=0)
         currentmat = filename
        
         #divide wavelength, refractiveindex and extinctioncoeff into three matricies
         materialinfo = materialinfo.rename(columns={materialinfo.columns[0]:'wavelength',
                                                       materialinfo.columns[1]:'refractiveindex',
                                                       materialinfo.columns[2]:'extinctioncoeff',})
         wavelength = materialinfo['wavelength']*1000 #wavelength in nm
         refractiveindex = materialinfo['refractiveindex']
         extinctioncoeff = materialinfo['extinctioncoeff']

         #data interpolatinon
         interpolation_granularity = 100
         #wavelength vs refractive index
         wavelength_int = [] #interpolated wavelength for first material
         refractiveindex_int = [] #interpolated refractive index for first material
         for i in range(len(wavelength)-1):
             slope = (refractiveindex[i+1] - refractiveindex[i]) / (wavelength[i+1] - wavelength[i])
             intercept = refractiveindex[i] - slope * wavelength[i]
             wavelength_int_temp = np.linspace(wavelength[i], wavelength[i+1], interpolation_granularity)
             refractiveindex_int_temp = slope * wavelength_int_temp + intercept

             wavelength_int=np.append(wavelength_int,wavelength_int_temp)
             refractiveindex_int=np.append(refractiveindex_int,refractiveindex_int_temp)

         #wavelength vs extinction coeff
         wavelength_int = [] #interpolated wavelength for first material
         extinctioncoeff_int = [] #interpolated refractive index for first material
         for i in range(len(wavelength)-1):
             slope = (extinctioncoeff[i+1] - extinctioncoeff[i]) / (wavelength[i+1] - wavelength[i])
             intercept = extinctioncoeff[i] - slope * wavelength[i]
             wavelength_int_temp = np.linspace(wavelength[i], wavelength[i+1], interpolation_granularity)
             extinctioncoeff_int_temp = slope * wavelength_int_temp + intercept
             
             wavelength_int=np.append(wavelength_int,wavelength_int_temp)
             extinctioncoeff_int=np.append(extinctioncoeff_int,extinctioncoeff_int_temp)
             
         #searching function
         #value searching function for material 1
         def find_nearest_value(refractiveindex, extinctioncoeff, datawavelengths, Target_Wavelength):
             #calculated difference between target wavelength and all wavelengths
             abs_diff = np.abs(datawavelengths - Target_Wavelength)
             #find the value with minimum difference
             min_index = np.argmin(abs_diff)
             #return values
             return refractiveindex[min_index], extinctioncoeff[min_index], datawavelengths[min_index]
         
         refractiveindex,extinctioncoeff,wavelength = find_nearest_value(refractiveindex_int, extinctioncoeff_int, 
                                                        wavelength_int, Target_Wavelength)
         
         #Permanent assignment of variables
         filename = filename.replace(".csv", " ", 1)
                
         if currentmat == mat1:
             n_mat1 = np.append(n_mat1,refractiveindex)
             k_mat1 = np.append(k_mat1,extinctioncoeff)
             
         if currentmat == mat2:
             n_mat2 = np.append(n_mat2,refractiveindex)
             k_mat2 = np.append(k_mat2,extinctioncoeff)


# thickness = 850e-9 / 2 / 3.3 #layer thickness calculated at 850 nm
# tally1 = 0 #i is not working in the loop for some reason so need this here to count the loops
# for i in wavelength_spectra:
#      Target_Wavelength = wavelength_spectra[tally1]
#      delta_index = n_mat1[tally1] - n_mat2[tally1]
#      index_1 = n_mat1[tally1]
#      index_2 = n_mat2[tally1]
     
#      layers = []
#      for i in range(n_pairs):
#          layer1 = (index_1, thickness / 2)
#          layer2 = (index_2, thickness / 2)
#          layers.append(layer1)
#          layers.append(layer2)
     
#      tally1 = tally1 + 1
     
#      layers.append((1, 0))
#      angle_inc = 0
#      transfer_matrix_TE_system, transfer_matrix_TM_system = TM(layers, Target_Wavelength, angle_inc)
#      system_TE_trans, system_TE_reflect = TR(transfer_matrix_TE_system)
#      system_TM_trans, system_TM_reflect = TR(transfer_matrix_TM_system)
reflec_spectra = []
thickness = 850e-9 / 2 / 3.4  # lambda / 4 condition (per layer)

tallytest = 0
for i in wavelength_spectra:
     Target_Wavelength = wavelength_spectra[tallytest]
     
     delta_index = n_mat1[tallytest] - n_mat2[tallytest]
     print(delta_index)
     index_1 = n_mat2[tallytest]
     index_2 = n_mat1[tallytest]
     average_index = (index_1 + index_2 )/2

     # delta_index = 0.55
     # index_1 = 3.3 - 0.55/2
     # index_2 = 3.3 + 0.55/2
     # average_index = 3.3

     #layers = [(index_2, wavelength / index_2)]
     layers = []

     for i in range(n_pairs):
         layer1 = (index_1, thickness / 2)
         layer2 = (index_2, thickness / 2)
         layers.append(layer1)
         layers.append(layer2)

     # Exit to air
     layers.append((1, 0))

     angle_inc = 0
     transfer_matrix_TE_system, transfer_matrix_TM_system = TM(layers, Target_Wavelength*10**(-9), angle_inc)
     system_TE_trans, system_TE_reflect = TR(transfer_matrix_TE_system)
     
     # print("Total Transfer Matrix of the system:")
     # print(transfer_matrix_TE_system)
     # print(f"Transmittance TE: {system_TE_trans}, Reflectance TE:  {system_TE_reflect}.")
     # print(f"Transmittance TM: {system_TM_trans}, Reflectance TM:  {system_TM_reflect}.")
     
     reflec_spectra = np.append(reflec_spectra,system_TE_reflect)
     
     tallytest = tallytest + 1

plt.plot(wavelength_spectra,reflec_spectra)
plt.show
