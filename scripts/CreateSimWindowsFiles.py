# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:41:56 2024

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
import shutil

# Tunable variables
wavelength_design = 850e-9; # meters
n_pairs_top = 15;
n_pairs_bottom = 30.5;
charge_density_e = 2.5e18 * 1e6 # Electrons
charge_density_h = 5e18 * 1e6 # Holes
DeltaDoping = True
charge_density_delta_e = 2.5e19 * 1e6 # Electrons
charge_density_delta_h = 5e19 * 1e6 # Holes
delta_doping_thickness_factor = 1 / 20
quantum_well_thickness = 8e-3 # microns
n_quantum_wells = 1
substrate_thickness = 0.1; # microns
coarse_points_per_nm = 1
fine_points_per_nm = 10
target_reflectance_top = 0.95
target_reflectance_bottom = 0.995
cavity_radius = 3.5 # microns
# Plotting options
plt.rcParams.update({'font.size': 18})

# Define spectrum
wavelength_range = 4e-9; # meters
wavelegnth_step = 1e-9; # meters
# Vector with testing wavelengths
wavelengths = np.arange(wavelength_design - wavelength_range / 2,
                        wavelength_design + wavelength_range / 2, wavelegnth_step)
wavelength_design_index = np.argwhere(abs(wavelengths - wavelength_design) < 1e-16)[0][0]

# Path to materials files
path_to_mat_files = '../materials'
mat_file_GaAs = 'GaAs'
Al_concentration = 0.097 # 0.219
mat_file_AlGaAs = 'AlGaAs-X=' + str(Al_concentration)
mat_file_AlAs = 'AlAs'

# Design parameters
thickness_error = 0.00
sigmaFactor = 2

# Load materials parameters, n and K for all wavelength
index_Mat_GaAs_undoped, loss_Mat_GaAs_undoped = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_GaAs), False)
index_Mat_AlGaAs_undoped, loss_Mat_AlGaAs_undoped = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_AlGaAs))
index_Mat_AlAs_undoped, loss_Mat_AlAs_undoped = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_AlAs))

# Doping properties
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
n_pairs = np.array([n_pairs_top]) # np.arange(2, 41, 1) #np.array([15])
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
num_runs = 1

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
    show_spectral_reflectance = False
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
n_pairs = np.array([n_pairs_bottom]) # np.arange(2, 41, 1) + 0.5
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
index_Mat_3_design = index_Mat_out[wavelength_design_index]
loss_Mat_3_design = loss_Mat_out[wavelength_design_index]

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
    show_spectral_reflectance = False
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

# %% Compute parameters of interest
# Number of quantum barriers
n_quantum_barriers = n_quantum_wells - 1
n_quantum_layers = n_quantum_wells + n_quantum_barriers
# Average refractive indices
average_index_top = (index_Mat_AlAs_p_doped[wavelength_design_index] + index_Mat_AlGaAs_p_doped[wavelength_design_index]) / 2
average_index_bottom = (index_Mat_AlAs_n_doped[wavelength_design_index] + index_Mat_AlGaAs_n_doped[wavelength_design_index]) / 2
confinement_index =  index_Mat_AlGaAs_undoped[wavelength_design_index]
quantum_well_index =  index_Mat_GaAs_undoped[wavelength_design_index]
# Ideal thicknesses
digit_prec = 8
ideal_layer_thickness_top = round(wavelength_design / 4 / average_index_top * 1e6, digit_prec) # lambda / 4 condition (per layer)
ideal_thickness_top = round(ideal_layer_thickness_top * n_pairs_top * 2, digit_prec)
ideal_layer_thickness_bottom = round(wavelength_design / 4 / average_index_bottom * 1e6, digit_prec) # lambda / 4 condition (per layer)
ideal_thickness_bottom = round(ideal_layer_thickness_bottom * n_pairs_bottom * 2, digit_prec)
ideal_thickness_confinement = round((wavelength_design - (quantum_well_index * quantum_well_thickness * n_quantum_layers) * 1e-6) / 2 / confinement_index * 1e6, digit_prec)
ideal_thickness_bulk_top = round(ideal_thickness_top + ideal_thickness_confinement, digit_prec)
ideal_thickness_bulk_bottom = round(ideal_thickness_bottom + ideal_thickness_confinement + substrate_thickness, digit_prec)
ideal_thickness_cavity = round(2 * ideal_thickness_confinement + quantum_well_thickness * n_quantum_layers, digit_prec)
total_thickness = round(ideal_layer_thickness_top * 2 * n_pairs_top +
                  ideal_layer_thickness_bottom * 2 * (n_pairs_bottom) + 
                  ideal_thickness_cavity, digit_prec)
total_thickness = round(ideal_thickness_top + ideal_thickness_cavity + ideal_thickness_bottom + substrate_thickness, digit_prec)
# Cavity area surface
cavity_area = round(np.pi * cavity_radius ** 2, digit_prec)

# Refractive indices
digit_prec = 8
# GaAs
index_GaAs_undoped = round(index_Mat_GaAs_undoped[wavelength_design_index], digit_prec)
index_GaAs_n_doped = round(index_Mat_GaAs_n_doped[wavelength_design_index], digit_prec)
index_GaAs_p_doped = round(index_Mat_GaAs_p_doped[wavelength_design_index], digit_prec)
# AlGaAs
index_AlGaAs_undoped = round(index_Mat_AlGaAs_undoped[wavelength_design_index], digit_prec)
index_AlGaAs_n_doped = round(index_Mat_AlGaAs_n_doped[wavelength_design_index], digit_prec)
index_AlGaAs_p_doped = round(index_Mat_AlGaAs_p_doped[wavelength_design_index], digit_prec)
# AlAs
index_AlAs_undoped = round(index_Mat_AlAs_undoped[wavelength_design_index], digit_prec)
index_AlAs_n_doped = round(index_Mat_AlAs_n_doped[wavelength_design_index], digit_prec)
index_AlAs_p_doped = round(index_Mat_AlAs_p_doped[wavelength_design_index], digit_prec)

# Extinction coffecicients
# GaAs
loss_GaAs_undoped = loss_Mat_GaAs_undoped[wavelength_design_index]
loss_GaAs_n_doped = loss_Mat_GaAs_n_doped[wavelength_design_index]
loss_GaAs_p_doped = loss_Mat_GaAs_p_doped[wavelength_design_index]
# AlGaAs
loss_AlGaAs_undoped = loss_Mat_AlGaAs_undoped[wavelength_design_index]
loss_AlGaAs_n_doped = loss_Mat_AlGaAs_n_doped[wavelength_design_index]
loss_AlGaAs_p_doped = loss_Mat_AlGaAs_p_doped[wavelength_design_index]
# AlAs
loss_AlAs_undoped = loss_Mat_AlAs_undoped[wavelength_design_index]
loss_AlAs_n_doped = loss_Mat_AlAs_n_doped[wavelength_design_index]
loss_AlAs_p_doped = loss_Mat_AlAs_p_doped[wavelength_design_index]

# Absorption coefficients
def ExtinctioToAbsoprtion(extinction, wavelength):
    return 4 * np.pi * extinction / wavelength * 1e-2
# GaAs
absorption_GaAs_undoped = round(ExtinctioToAbsoprtion(loss_Mat_GaAs_undoped[wavelength_design_index], wavelength_design), digit_prec)
absorption_GaAs_n_doped = round(ExtinctioToAbsoprtion(loss_Mat_GaAs_n_doped[wavelength_design_index], wavelength_design), digit_prec)
absorption_GaAs_p_doped = round(ExtinctioToAbsoprtion(loss_Mat_GaAs_p_doped[wavelength_design_index], wavelength_design), digit_prec)
# AlGaAs
absorption_AlGaAs_undoped = round(ExtinctioToAbsoprtion(loss_Mat_AlGaAs_undoped[wavelength_design_index], wavelength_design), digit_prec)
absorption_AlGaAs_n_doped = round(ExtinctioToAbsoprtion(loss_Mat_AlGaAs_n_doped[wavelength_design_index], wavelength_design), digit_prec)
absorption_AlGaAs_p_doped = round(ExtinctioToAbsoprtion(loss_Mat_AlGaAs_p_doped[wavelength_design_index], wavelength_design), digit_prec)
# AlAs
absorption_AlAs_undoped = round(ExtinctioToAbsoprtion(loss_Mat_AlAs_undoped[wavelength_design_index], wavelength_design), digit_prec)
absorption_AlAs_n_doped = round(ExtinctioToAbsoprtion(loss_Mat_AlAs_n_doped[wavelength_design_index], wavelength_design), digit_prec)
absorption_AlAs_p_doped = round(ExtinctioToAbsoprtion(loss_Mat_AlAs_p_doped[wavelength_design_index], wavelength_design), digit_prec)

# Print info
print(f"Total thcikness: {total_thickness} um")
print(f"Top DBR layer thickness: {ideal_layer_thickness_top} um")
print(f"Top DBR thickness: {ideal_thickness_top} um")
print(f"Bottom DBR layer thickness: {ideal_layer_thickness_bottom} um")
print(f"Top DBR thickness: {ideal_thickness_bottom} um")
print(f"Confinement layer thickness: {ideal_thickness_confinement} um")
print(f"Top DBR + confinement layer thickness: {ideal_thickness_bulk_top} um")
print(f"Bottom DBR + confinement layer thickness + substrate: {ideal_thickness_bulk_bottom} um")
print(f"Total cavity thickness: {ideal_thickness_cavity} um")
print(f"Total VCSEL thickness: {total_thickness} um")

# %% Create device file
def WriteGrid(file, length, n_pints):
    file.write(f"grid length={length} points={n_pints:.0f}\n")
    return

def WriteRegion(file, length, region_type):
    file.write(f"region {region_type} length={length}\n")
    return

def WriteLayer(file, thickness, refractive_index, absorption, alloy, concentration):
    #if refractive_index > 0:
    #    file.write(f"Refractive_Index length={thickness} value={refractive_index}\n")
    #if absorption > 0:
    #    file.write(f"Absorption length={thickness} value={0}\n")
    file.write(f"structure material=gaas alloy={alloy} length={thickness} conc={concentration}\n")
    return

def WriteDoping(file, thickness, donors_conc, acceptors_conc):
    file.write(f"doping length={thickness}")
    if donors_conc > 0:
        file.write(f" Nd={donors_conc} Nd_deg=2 Nd_level=0.005")
    if acceptors_conc > 0:
        file.write(f" Na={acceptors_conc} Na_deg=4 Na_level=0.026")
    file.write("\n")
    return

def WriteCavity(file, radius, area, length, top_relfectance, bottom_reflectance, total_thickness):
    file.write(f"radius={radius}\n")
    file.write(f"cavity surface area={area} length={length}\n")
    file.write(f"mirror metal position=0 ref={top_relfectance}\n")
    file.write(f"mirror metal position={total_thickness} ref={bottom_reflectance}\n")
    return

device_file_name = 'VCSEL.DEV'
with open(f"../devices/{device_file_name}", "w") as file:
    # Write laser cavity parameters
    WriteCavity(file, cavity_radius, cavity_area, ideal_thickness_cavity, target_reflectance_top, target_reflectance_bottom, total_thickness)
    file.write("\n")
    # Start with bulk region
    WriteRegion(file, ideal_thickness_bulk_top, "bulk")
    # Create top DBR grid
    WriteGrid(file, ideal_thickness_top, ideal_thickness_top * 1e3 * coarse_points_per_nm)
    file.write("\n")
    # Top DBR
    for i in range(n_pairs_top):
        WriteLayer(file, ideal_layer_thickness_top, index_AlGaAs_p_doped, absorption_AlGaAs_p_doped, "p-AlGaAs", Al_concentration)
        WriteLayer(file, ideal_layer_thickness_top, index_AlAs_p_doped, absorption_AlGaAs_p_doped, "p-AlAs", 1)
        if DeltaDoping:
            WriteDoping(file, round(ideal_layer_thickness_top * delta_doping_thickness_factor, digit_prec), 0, charge_density_delta_h * 1e-6)
            WriteDoping(file, round(ideal_layer_thickness_top * (1 - 2 * delta_doping_thickness_factor), digit_prec), 0, charge_density_h * 1e-6)
            WriteDoping(file, round(2 * ideal_layer_thickness_top * delta_doping_thickness_factor, digit_prec), 0, charge_density_delta_h * 1e-6)
            WriteDoping(file, round(ideal_layer_thickness_top * (1 - 2 * delta_doping_thickness_factor), digit_prec), 0, charge_density_h * 1e-6)
            WriteDoping(file, round(ideal_layer_thickness_top * delta_doping_thickness_factor, digit_prec), 0, charge_density_delta_h * 1e-6)
    
    file.write("\n")
    # Top mirror is p-type
    if not DeltaDoping:
        WriteDoping(file, ideal_thickness_top, 0, charge_density_h * 1e-6)
    file.write("\n")
    # Create cavity grid
    WriteGrid(file, ideal_thickness_cavity, ideal_thickness_cavity * 1e3 * fine_points_per_nm)
    file.write("\n")
    # Top confinement layer
    WriteLayer(file, ideal_thickness_confinement, index_AlGaAs_undoped, absorption_AlGaAs_undoped, "AlGaAs", Al_concentration)
    file.write("\n")
    # Quantum well (change to qw region)
    for i in range(n_quantum_layers):
        if i % 2 == 0:
            WriteRegion(file, quantum_well_thickness, "qw")
            WriteLayer(file, quantum_well_thickness, index_GaAs_undoped, index_GaAs_undoped, "undoped", 0)
        else:
            WriteRegion(file, quantum_well_thickness, "bulk")
            WriteLayer(file, quantum_well_thickness, index_AlGaAs_undoped, index_AlGaAs_undoped, "AlGaAs", Al_concentration)
    file.write("\n")
    # Continue with bulk region
    WriteRegion(file, ideal_thickness_bulk_bottom, "bulk")
    # Bottom confinement layer
    WriteLayer(file, ideal_thickness_confinement, index_AlGaAs_undoped, absorption_AlGaAs_undoped, "AlGaAs", Al_concentration)
    file.write("\n")
    # Cavity is undoped
    #WriteDoping(file, ideal_thickness_confinement, 0, 1e17)
    WriteDoping(file, ideal_thickness_confinement, 0, 0)
    WriteDoping(file, quantum_well_thickness * n_quantum_layers, 0, 0)
    WriteDoping(file, ideal_thickness_confinement, 0, 0)
    #WriteDoping(file, ideal_thickness_cavity, 0, 0)
    # Create bottom DBR grid
    WriteGrid(file, ideal_thickness_bottom, ideal_thickness_bottom * 1e3 * coarse_points_per_nm)
    file.write("\n")
    # Bottom DBR
    for i in range(int(np.ceil(n_pairs_bottom))):
        WriteLayer(file, ideal_layer_thickness_bottom, index_AlAs_n_doped, absorption_AlAs_n_doped, "n-AlAs", 1)
        if DeltaDoping:
            WriteDoping(file, round(ideal_layer_thickness_bottom * delta_doping_thickness_factor, digit_prec), charge_density_delta_e * 1e-6, 0)
            WriteDoping(file, round(ideal_layer_thickness_bottom * (1 - 2 *delta_doping_thickness_factor), digit_prec), charge_density_e * 1e-6, 0)
            WriteDoping(file, round(ideal_layer_thickness_bottom * delta_doping_thickness_factor, digit_prec), charge_density_delta_e * 1e-6, 0)
        if i <= n_pairs - 1:
            WriteLayer(file, ideal_layer_thickness_bottom, index_AlGaAs_n_doped, absorption_AlGaAs_n_doped, "n-AlGaAs", Al_concentration)
            if DeltaDoping:
                WriteDoping(file, round(ideal_layer_thickness_bottom * delta_doping_thickness_factor, digit_prec), charge_density_delta_e * 1e-6, 0)
                WriteDoping(file, round(ideal_layer_thickness_bottom * (1 - 2 * delta_doping_thickness_factor), digit_prec), charge_density_e * 1e-6, 0)
                WriteDoping(file, round(ideal_layer_thickness_bottom * delta_doping_thickness_factor, digit_prec), charge_density_delta_e * 1e-6, 0)
    file.write("\n")
    # Create substrate grid
    WriteGrid(file, substrate_thickness, substrate_thickness * 1e3 * coarse_points_per_nm)
    file.write("\n")
    # Subtrate
    WriteLayer(file, substrate_thickness, index_GaAs_n_doped, absorption_GaAs_n_doped, "n-GaAs", 0)
    file.write("\n")
    # Bottom mirror is n-type
    if DeltaDoping:
        WriteDoping(file, round(substrate_thickness, digit_prec), charge_density_e * 1e-6, 0)
    else:
        WriteDoping(file, round(ideal_thickness_bottom + substrate_thickness, digit_prec), charge_density_e * 1e-6, 0)
        
# Create material file
# Define paths relative to the script's location
script_dir = os.path.dirname(__file__)  # Directory where the script is located
source_file = os.path.join(script_dir, '../materials/MATERIAL_raw.PRM')
destination_file = os.path.join(script_dir, '../materials/MATERIAL.PRM')
# Normalize paths for safety
source_file = os.path.abspath(source_file)
destination_file = os.path.abspath(destination_file)
# Copy raw file
shutil.copy(source_file, destination_file)

# Add text to specific lines
lines_to_edit = {187: f"\nREFRACTIVE_INDEX value={index_GaAs_undoped}\nABSORPTION value={0*absorption_GaAs_undoped}\n",
                 223: f"\nREFRACTIVE_INDEX value={index_AlAs_undoped}\nABSORPTION value={0*absorption_AlAs_undoped}\n",
                 277: f"\nREFRACTIVE_INDEX value={index_AlGaAs_undoped}\nABSORPTION value={0*absorption_AlGaAs_undoped}\n",
                 325: f"\nREFRACTIVE_INDEX value={index_GaAs_p_doped}\nABSORPTION value={0*absorption_GaAs_p_doped}\n",
                 361: f"\nREFRACTIVE_INDEX value={index_AlAs_p_doped}\nABSORPTION value={0*absorption_AlAs_p_doped}\n",
                 415: f"\nREFRACTIVE_INDEX value={index_AlGaAs_p_doped}\nABSORPTION value={0*absorption_AlGaAs_p_doped}\n",
                 463: f"\nREFRACTIVE_INDEX value={index_GaAs_n_doped}\nABSORPTION value={0*absorption_GaAs_n_doped}\n",
                 499: f"\nREFRACTIVE_INDEX value={index_AlAs_n_doped}\nABSORPTION value={0*absorption_AlAs_n_doped}\n",
                 553: f"\nREFRACTIVE_INDEX value={index_AlGaAs_n_doped}\nABSORPTION value={0*absorption_AlGaAs_n_doped}\n"}
# Read and modify lines
with open(destination_file, 'r', encoding='latin1') as file:
    lines = file.readlines()
    
for line_number, text_to_add in lines_to_edit.items():
    if 1 <= line_number <= len(lines):
        lines[line_number - 1] = lines[line_number - 1].rstrip() + text_to_add + '\n'
# Write the changes back to the file
with open(destination_file, 'w', encoding='latin1') as file:
    file.writelines(lines)
    
# Copy files to simWindows folder
simWindows_path = os.path.join('C:\SimWindows')
shutil.copy(destination_file, os.path.join(simWindows_path, 'MATERIAL.PRM'))
device_source_file = os.path.abspath(os.path.join(script_dir, f'../devices/{device_file_name}'))
shutil.copy(device_source_file, os.path.join(simWindows_path, 'VCSEL', device_file_name))
