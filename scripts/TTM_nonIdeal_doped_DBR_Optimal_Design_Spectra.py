# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:08:14 2024

@author: srulo
"""

import numpy as np
from matplotlib import pyplot as plt
from TTM import Transfer_matrix_system as TM
from TTM import Trans_and_reflec as TR
from MatFileLoader import Read_mat_file as ReadMatFile
from MatFileLoader import Relative_permitivity_free_charge_carrier as RelPermFCC
from MatFileLoader import Read_mat_file_Add_DSMOdel as ReadMatFileFCC
from CreateDBR import Create_DBR as DBR
import os

# Plotting options
plt.rcParams.update({'font.size': 18})

# Define spectrum
wavelength_design = 850e-9; # meters
wavelength_range = 40e-9; # meters
wavelegnth_step = 1e-9; # meters
# Vector with testing wavelengths
wavelengths = np.arange(wavelength_design - wavelength_range / 2,
                        wavelength_design + wavelength_range / 2, wavelegnth_step)
wavelength_design_index = np.argwhere(abs(wavelengths - wavelength_design) < 1e-16)[0][0]

# Path to materials files
path_to_mat_files = '../materials'
mat_file_GaAs = 'GaAs'
Al_concentration = 0.097 # 0.097
mat_file_AlGaAs = 'AlGaAs-X=' + str(Al_concentration)
mat_file_AlAs = 'AlAs'

# Design parameters
target_reflectance_top = 0.95
target_reflectance_bottom = 0.995
thickness_error = 0.05
sigmaFactor = 2

# Load materials parameters, n and K for all wavelength
index_Mat_GaAs_undoped, loss_Mat_GaAs_undoped = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_GaAs), False)
index_Mat_AlGaAs_undoped, loss_Mat_AlGaAs_undoped = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_AlGaAs))
index_Mat_AlAs_undoped, loss_Mat_AlAs_undoped = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_AlAs))

# Doping properties
charge_density_e = 2.5e18 * 1e6 # Electrons
charge_density_h = 5e18 * 1e6 # Holes
effective_mass_e_GaAs = 0.063 # Electrons
effective_mass_h_GaAs = 0.57 # Holes
effective_mass_e_AlAs = 0.15 # Electrons
effective_mass_h_AlAs = 0.76 # Holes
effective_mass_e_AlGaAs = 0.063 + (0.087 * Al_concentration) # Electrons
effective_mass_h_AlGaAs = 0.57 + (0.19 * Al_concentration) # Holes
damping_freq_e_GaAs = 1.03e13
damping_freq_h_GaAs = 1.71e13
damping_freq_e_AlGaAs = 9.24e12 # ?????
damping_freq_h_AlGaAs = 1.76e13
damping_freq_e_AlAs = 6.51e13
damping_freq_h_AlAs = 2.31e13

# Total optical properties
index_Mat_GaAs_p_doped, loss_Mat_GaAs_p_doped = ReadMatFileFCC(wavelengths, os.path.join(path_to_mat_files, mat_file_GaAs), 
                                               charge_density_h, effective_mass_h_GaAs, damping_freq_h_GaAs, False)
index_Mat_AlGaAs_p_doped, loss_Mat_AlGaAs_p_doped = ReadMatFileFCC(wavelengths, os.path.join(path_to_mat_files, mat_file_AlGaAs), 
                                               charge_density_h, effective_mass_h_AlGaAs, damping_freq_h_AlGaAs, False)
index_Mat_AlAs_p_doped, loss_Mat_AlAs_p_doped = ReadMatFileFCC(wavelengths, os.path.join(path_to_mat_files, mat_file_AlAs), 
                                               charge_density_h, effective_mass_h_AlAs, damping_freq_h_AlAs, False)
index_Mat_GaAs_n_doped, loss_Mat_GaAs_n_doped = ReadMatFileFCC(wavelengths, os.path.join(path_to_mat_files, mat_file_GaAs), 
                                               charge_density_e, effective_mass_e_GaAs, damping_freq_e_GaAs, False)
index_Mat_AlGaAs_n_doped, loss_Mat_AlGaAs_n_doped = ReadMatFileFCC(wavelengths, os.path.join(path_to_mat_files, mat_file_AlGaAs), 
                                               charge_density_e, effective_mass_e_AlGaAs, damping_freq_e_AlGaAs, False)
index_Mat_AlAs_n_doped, loss_Mat_AlAs_n_doped = ReadMatFileFCC(wavelengths, os.path.join(path_to_mat_files, mat_file_AlAs), 
                                               charge_density_e, effective_mass_e_AlAs, damping_freq_e_AlAs, False)

### TOP MIRROR
# Vector with testing number of layer pairs
n_pairs = np.array([15]) #np.arange(2, 41, 1) #
# Set indices of refraction for input layer
index_Mat_in = index_Mat_AlGaAs_undoped
loss_Mat_in = loss_Mat_AlGaAs_undoped
# Set indices of refraction for alternating layers
index_Mat_1 = index_Mat_AlAs_p_doped
loss_Mat_1 = loss_Mat_AlAs_p_doped
index_Mat_2 = index_Mat_AlGaAs_p_doped
loss_Mat_2 = loss_Mat_AlGaAs_p_doped
# Set indices of refraction for output layer
index_Mat_out = np.ones(wavelengths.shape[0]) # Air
loss_Mat_out = np.zeros(wavelengths.shape[0]) # Air

# Set parameters for design wavelength
index_Mat_1_design = index_Mat_1[wavelength_design_index]
loss_Mat_1_design = loss_Mat_1[wavelength_design_index]
index_Mat_2_design = index_Mat_2[wavelength_design_index]
loss_Mat_2_design = loss_Mat_2[wavelength_design_index]

# Average index of refraction
average_index = (index_Mat_1_design + index_Mat_2_design) / 2
thickness_ideal = wavelength_design / 4 / average_index # lambda / 4 condition (per layer)
# Normal incidence
angle_inc = 0
# Number of runs to emulate non-ideal thicknesses
num_runs = 10

# Iterate over design configurations
# Initialize variables
DBR_TE_trans = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TE_reflect = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TM_trans = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TM_reflect = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TE_mean_trans_design = []
DBR_TE_mean_reflect_design = []
DBR_TM_mean_trans_design = []
DBR_TM_mean_reflect_design = []
DBR_TE_std_trans_design = []
DBR_TE_std_reflect_design = []
DBR_TM_std_trans_design = []
DBR_TM_std_reflect_design = []

n_pair_index = 0;
for n_pair in n_pairs:
    wavelength_index = 0;
    for wavelength in wavelengths:
        # Assing index of refractions and losses for this wavelength
        index_in = index_Mat_in[wavelength_index] + 1j * loss_Mat_in[wavelength_index]
        index_1 = index_Mat_1[wavelength_index] + 1j * loss_Mat_1[wavelength_index]
        index_2 = index_Mat_2[wavelength_index] + 1j * loss_Mat_2[wavelength_index]
        index_out = index_Mat_out[wavelength_index] + 1j * loss_Mat_out[wavelength_index]
        # Create DBR
        DBR_TE_trans_run = []
        DBR_TE_reflect_run = []
        DBR_TM_trans_run = []
        DBR_TM_reflect_run = []
        for run in range(num_runs):
            DBR_ = DBR(n_pair, index_in, index_1, index_2, thickness_ideal, thickness_error, index_out)
            # Transfer matrix (wavelegnth in meters!)
            transfer_matrix_TE_DBR, transfer_matrix_TM_DBR = TM(DBR_, wavelength, angle_inc)
            # Transmittance and Reflectance
            DBR_TE_trans_, DBR_TE_reflect_ = TR(transfer_matrix_TE_DBR)
            DBR_TM_trans_, DBR_TM_reflect_ = TR(transfer_matrix_TM_DBR)
            DBR_TE_trans_run.append(DBR_TE_trans_)
            DBR_TE_reflect_run.append(DBR_TE_reflect_)
            DBR_TM_trans_run.append(DBR_TM_trans_)
            DBR_TM_reflect_run.append(DBR_TM_reflect_)
            
        # Mean values of runs
        DBR_TE_trans[n_pair_index, wavelength_index, 0] = np.mean(DBR_TE_trans_run)
        DBR_TE_reflect[n_pair_index, wavelength_index, 0] = np.mean(DBR_TE_reflect_run)
        DBR_TM_trans[n_pair_index, wavelength_index, 0] = np.mean(DBR_TM_trans_run)
        DBR_TM_reflect[n_pair_index, wavelength_index, 0] = np.mean(DBR_TM_reflect_run)
        # Standard deviation over runs
        DBR_TE_trans[n_pair_index, wavelength_index, 1] = np.std(DBR_TE_trans_run)
        DBR_TE_reflect[n_pair_index, wavelength_index, 1] = np.std(DBR_TE_reflect_run)
        DBR_TM_trans[n_pair_index, wavelength_index, 1] = np.std(DBR_TM_trans_run)
        DBR_TM_reflect[n_pair_index, wavelength_index, 1] = np.std(DBR_TM_reflect_run)
        # Next wavelength
        wavelength_index += 1
    
    # Assign reflectance for design wavelength
    DBR_TE_mean_trans_design.append(DBR_TE_trans[n_pair_index, wavelength_design_index, 0])
    DBR_TE_mean_reflect_design.append(DBR_TE_reflect[n_pair_index, wavelength_design_index, 0])
    DBR_TM_mean_trans_design.append(DBR_TM_trans[n_pair_index, wavelength_design_index, 0])
    DBR_TM_mean_reflect_design.append(DBR_TM_reflect[n_pair_index, wavelength_design_index, 0])
    DBR_TE_std_trans_design.append(DBR_TE_trans[n_pair_index, wavelength_design_index, 1])
    DBR_TE_std_reflect_design.append(DBR_TE_reflect[n_pair_index, wavelength_design_index, 1])
    DBR_TM_std_trans_design.append(DBR_TM_trans[n_pair_index, wavelength_design_index, 1])
    DBR_TM_std_reflect_design.append(DBR_TM_reflect[n_pair_index, wavelength_design_index, 1])
    
    # Plot spectral reflectance for this design
    show_spectral_reflectance = True
    if show_spectral_reflectance:
        plt.figure()
        plt.plot(wavelengths*1e9, DBR_TE_reflect[n_pair_index, :, 0],'k',label='R')
        plt.plot(wavelengths*1e9, DBR_TE_trans[n_pair_index, :, 0],'#87aadeff',label='T')
        absorpt = 1 - DBR_TE_reflect[n_pair_index, :, 0] - DBR_TE_trans[n_pair_index, :, 0]
        plt.plot(wavelengths*1e9, absorpt,'#ff5555ff',label='A')
        plt.axvline(x = wavelength_design*1e9, linestyle = '--', color = 'k')
        plt.axhline(y = target_reflectance_top, linestyle = '--', color = 'k')
        plt.xlim(wavelengths[0]*1e9, wavelengths[-1]*1e9)
        plt.ylim(0, 1.05)
        plt.grid()
        plt.xlabel('Wavelegnth (nm)')
        plt.ylabel('Optical parameters')
        plt.title('Tom DBR mirror')
        plt.legend(loc='upper left')
        plt.show()
    
    # Print this reflectance
    this_DBR_TE_reflect = DBR_TE_mean_reflect_design[n_pair_index]
    this_DBR_TE_std_reflect = DBR_TE_std_reflect_design[n_pair_index]
    this_DBR_TE_reflect_min = this_DBR_TE_reflect - sigmaFactor * this_DBR_TE_std_reflect
    this_DBR_TE_reflect_max = this_DBR_TE_reflect + sigmaFactor * this_DBR_TE_std_reflect
    print(f"N layers: {n_pair} | Reflectance: {this_DBR_TE_reflect} +- {this_DBR_TE_std_reflect} = [{this_DBR_TE_reflect_min} - {this_DBR_TE_reflect_max}]")
    # Next N
    n_pair_index += 1

# Reflectance versus number of layers
DBR_TE_error_reflect_design = [sigmaFactor * x for x in DBR_TE_std_reflect_design]
plt.figure()
plt.errorbar(n_pairs, DBR_TE_mean_reflect_design, yerr=DBR_TE_error_reflect_design, fmt='k.', label='R_TE')
#plt.plot(n_pairs, DBR_TE_mean_reflect_design,'k.',label='R_TE')
plt.axhline(y = target_reflectance_top, linestyle = '--', color = 'k', label = 'Top DBR target reflectance')
plt.title(' Top DBR mirror')
plt.xlabel('Number of alternating layers pairs')
plt.ylabel('Reflectance')
plt.ylim(0, 1.05)
plt.xlim(0, 42)
plt.grid()
plt.show()

# %%
### BOTTOM MIRROR
# Vector with testing number of layer pairs
n_pairs = np.array([30.5]) # np.arange(2, 41, 1) + 0.5 # 
# Set indices of refraction for input layer
index_Mat_in = index_Mat_AlGaAs_undoped
loss_Mat_in = loss_Mat_AlGaAs_undoped
# Set indices of refraction for alternating layers
index_Mat_1 = index_Mat_AlAs_n_doped
loss_Mat_1 = loss_Mat_AlAs_n_doped
index_Mat_2 = index_Mat_AlGaAs_n_doped
loss_Mat_2 = loss_Mat_AlGaAs_n_doped
# Set indices of refraction for output layer
index_Mat_out = index_Mat_GaAs_n_doped # Substrate
loss_Mat_out = index_Mat_GaAs_n_doped # Substrate

# Set parameters for design wavelength
index_Mat_1_design = index_Mat_1[wavelength_design_index]
loss_Mat_1_design = loss_Mat_1[wavelength_design_index]
index_Mat_2_design = index_Mat_2[wavelength_design_index]
loss_Mat_2_design = loss_Mat_2[wavelength_design_index]

# Average index of refraction
average_index = (index_Mat_1_design + index_Mat_2_design) / 2
thickness_ideal = wavelength_design / 4 / average_index # lambda / 4 condition (per layer)
# Normal incidence
angle_inc = 0

# Iterate over design configurations
# Initialize variables
DBR_TE_trans = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TE_reflect = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TM_trans = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TM_reflect = np.zeros((n_pairs.shape[0], wavelengths.shape[0], 2))
DBR_TE_mean_trans_design = []
DBR_TE_mean_reflect_design = []
DBR_TM_mean_trans_design = []
DBR_TM_mean_reflect_design = []
DBR_TE_std_trans_design = []
DBR_TE_std_reflect_design = []
DBR_TM_std_trans_design = []
DBR_TM_std_reflect_design = []

n_pair_index = 0;
for n_pair in n_pairs:
    wavelength_index = 0;
    for wavelength in wavelengths:
        # Assing index of refractions and losses for this wavelength
        index_in = index_Mat_in[wavelength_index] + 1j * loss_Mat_in[wavelength_index]
        index_1 = index_Mat_1[wavelength_index] + 1j * loss_Mat_1[wavelength_index]
        index_2 = index_Mat_2[wavelength_index] + 1j * loss_Mat_2[wavelength_index]
        index_out = index_Mat_out[wavelength_index] + 1j * loss_Mat_out[wavelength_index]
        # Create DBR
        DBR_TE_trans_run = []
        DBR_TE_reflect_run = []
        DBR_TM_trans_run = []
        DBR_TM_reflect_run = []
        
        for run in range(num_runs):
            DBR_ = DBR(n_pair, index_in, index_1, index_2, thickness_ideal, thickness_error, index_out)
            # Transfer matrix (wavelegnth in meters!)
            transfer_matrix_TE_DBR, transfer_matrix_TM_DBR = TM(DBR_, wavelength, angle_inc)
            # Transmittance and Reflectance
            DBR_TE_trans_, DBR_TE_reflect_ = TR(transfer_matrix_TE_DBR)
            DBR_TM_trans_, DBR_TM_reflect_ = TR(transfer_matrix_TM_DBR)
            DBR_TE_trans_run.append(DBR_TE_trans_)
            DBR_TE_reflect_run.append(DBR_TE_reflect_)
            DBR_TM_trans_run.append(DBR_TM_trans_)
            DBR_TM_reflect_run.append(DBR_TM_reflect_)
            
        # Mean values of runs
        DBR_TE_trans[n_pair_index, wavelength_index, 0] = np.mean(DBR_TE_trans_run)
        DBR_TE_reflect[n_pair_index, wavelength_index, 0] = np.mean(DBR_TE_reflect_run)
        DBR_TM_trans[n_pair_index, wavelength_index, 0] = np.mean(DBR_TM_trans_run)
        DBR_TM_reflect[n_pair_index, wavelength_index, 0] = np.mean(DBR_TM_reflect_run)
        # Standard deviation over runs
        DBR_TE_trans[n_pair_index, wavelength_index, 1] = np.std(DBR_TE_trans_run)
        DBR_TE_reflect[n_pair_index, wavelength_index, 1] = np.std(DBR_TE_reflect_run)
        DBR_TM_trans[n_pair_index, wavelength_index, 1] = np.std(DBR_TM_trans_run)
        DBR_TM_reflect[n_pair_index, wavelength_index, 1] = np.std(DBR_TM_reflect_run)
        # Next wavelength
        wavelength_index += 1
    
    # Assign reflectance for design wavelength
    DBR_TE_mean_trans_design.append(DBR_TE_trans[n_pair_index, wavelength_design_index, 0])
    DBR_TE_mean_reflect_design.append(DBR_TE_reflect[n_pair_index, wavelength_design_index, 0])
    DBR_TM_mean_trans_design.append(DBR_TM_trans[n_pair_index, wavelength_design_index, 0])
    DBR_TM_mean_reflect_design.append(DBR_TM_reflect[n_pair_index, wavelength_design_index, 0])
    DBR_TE_std_trans_design.append(DBR_TE_trans[n_pair_index, wavelength_design_index, 1])
    DBR_TE_std_reflect_design.append(DBR_TE_reflect[n_pair_index, wavelength_design_index, 1])
    DBR_TM_std_trans_design.append(DBR_TM_trans[n_pair_index, wavelength_design_index, 1])
    DBR_TM_std_reflect_design.append(DBR_TM_reflect[n_pair_index, wavelength_design_index, 1])
    
    # Plot spectral reflectance for this design
    show_spectral_reflectance = True
    print(show_spectral_reflectance)
    if show_spectral_reflectance:
        plt.figure()
        plt.plot(wavelengths*1e9, DBR_TE_reflect[n_pair_index, :, 0],'k',label='R')
        plt.plot(wavelengths*1e9, DBR_TE_trans[n_pair_index, :, 0],'#87aadeff',label='T')
        absorpt = 1 - DBR_TE_reflect[n_pair_index, :, 0] - DBR_TE_trans[n_pair_index, :, 0]
        plt.plot(wavelengths*1e9, absorpt,'#ff5555ff',label='A')
        plt.axvline(x = wavelength_design*1e9, linestyle = '--', color = 'k')
        plt.axhline(y = target_reflectance_bottom, linestyle = '--', color = 'k')
        plt.xlim(wavelengths[0]*1e9, wavelengths[-1]*1e9)
        plt.ylim(0, 1.05)
        plt.grid()
        plt.xlabel('Wavelegnth (nm)')
        plt.ylabel('Optical parameters')
        plt.title('Bottom DBR mirror')
        plt.legend(loc='upper left')
        plt.show()
    
    # Print this reflectance
    this_DBR_TE_reflect = DBR_TE_mean_reflect_design[n_pair_index]
    this_DBR_TE_std_reflect = DBR_TE_std_reflect_design[n_pair_index]
    this_DBR_TE_reflect_min = this_DBR_TE_reflect - sigmaFactor * this_DBR_TE_std_reflect
    this_DBR_TE_reflect_max = this_DBR_TE_reflect + sigmaFactor * this_DBR_TE_std_reflect
    print(f"N layers: {n_pair} | Reflectance: {this_DBR_TE_reflect} +- {this_DBR_TE_std_reflect} = [{this_DBR_TE_reflect_min} - {this_DBR_TE_reflect_max}]")
    # Next N
    n_pair_index += 1
    
# Total tickness
total_thickness = thickness_ideal * (15 + 30.5) * 2 + wavelength_design / index_Mat_2_design
print(f"Total thcikness: {total_thickness}")

# Reflectance versus number of layers
DBR_TE_error_reflect_design = [sigmaFactor * x for x in DBR_TE_std_reflect_design]
plt.figure()
plt.errorbar(n_pairs, DBR_TE_mean_reflect_design, yerr=DBR_TE_error_reflect_design, fmt='k.', label='R_TE')
#plt.plot(n_pairs, DBR_TE_mean_reflect_design,'k.',label='R_TE')
plt.axhline(y = target_reflectance_bottom, linestyle = '--', color = 'k', label = 'Bottom DBR target reflectance')
plt.title('Bottom DBR mirror')
plt.xlabel('Number of layers pairs')
plt.ylabel('Reflectance')
plt.ylim(0, 1.05)
plt.xlim(0, 42)
plt.grid()
plt.show()