# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 16:27:21 2024

@author: srulo
"""

import numpy as np
from matplotlib import pyplot as plt
from TTM import Transfer_matrix_system as TM
from TTM import Trans_and_reflec as TR
from MatFileLoader import Read_mat_file as ReadMatFile
from CreateDBR import Create_DBR as DBR
import os

# Plotting options
plt.rcParams.update({'font.size': 18})

# Define spectrum
wavelength_design = 850e-9; # meters
wavelength_range = 400e-9; # meters
wavelegnth_step = 1e-9; # meters
# Vector with testing wavelengths
wavelengths = np.arange(wavelength_design - wavelength_range / 2,
                        wavelength_design + wavelength_range / 2, wavelegnth_step)
wavelength_design_index = np.argwhere(abs(wavelengths - wavelength_design) < 1e-16)[0][0]

# Path to materials files
path_to_mat_files = '../materials'
mat_file_GaAs = 'GaAs'
mat_file_AlGaAs = 'AlGaAs-X=0.219'
mat_file_AlAs = 'AlAs'

# Design parameters
target_reflectance_top = 0.95;
target_reflectance_bottom = 0.995;
thickness_error = 0.05
# Materials parameters, n and K for design wavelegnth
#index_Mat_GaAs_design, loss_Mat_GaAs_design = ReadMatFile(wavelength_design, os.path.join(path_to_mat_files, mat_file_GaAs))
#index_Mat_AlGaAs_design, loss_Mat_AlGaAs_design = ReadMatFile(wavelength_design, os.path.join(path_to_mat_files, mat_file_AlGaAs))
#index_Mat_AlAs_design, loss_Mat_AlAs_design = ReadMatFile(wavelength_design, os.path.join(path_to_mat_files, mat_file_AlAs))

# Load materials parameters, n and K for all wavelength
index_Mat_GaAs, loss_Mat_GaAs = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_GaAs), False)
index_Mat_AlGaAs, loss_Mat_AlGaAs = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_AlGaAs))
index_Mat_AlAs, loss_Mat_AlAs = ReadMatFile(wavelengths, os.path.join(path_to_mat_files, mat_file_AlAs))

# Plot spectral reflectance for this design
#plt.figure()
#plt.plot(wavelengths,index_Mat_1,label='n')
#plt.plot(wavelengths,index_Mat_2,label='R_TM')
#plt.show

### TOP MIRROR
# Vector with testing number of layer pairs
n_pairs = np.arange(2, 41, 1) #np.array([15])
# Set indices of refraction for input layer
index_Mat_in = index_Mat_AlGaAs
loss_Mat_in = loss_Mat_AlGaAs
# Set indices of refraction for alternating layers
index_Mat_1 = index_Mat_AlAs
loss_Mat_1 = loss_Mat_AlAs
index_Mat_2 = index_Mat_AlGaAs
loss_Mat_2 = loss_Mat_AlGaAs
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
num_runs = 100

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
    this_DBR_TE_reflect_min = this_DBR_TE_reflect - 4 * this_DBR_TE_std_reflect
    this_DBR_TE_reflect_max = this_DBR_TE_reflect + 4 * this_DBR_TE_std_reflect
    print(f"N layers: {n_pair} | Reflectance: {this_DBR_TE_reflect} +- {this_DBR_TE_std_reflect} = [{this_DBR_TE_reflect_min} - {this_DBR_TE_reflect_max}]")
    # Next N
    n_pair_index += 1

# Reflectance versus number of layers
DBR_TE_error_reflect_design = [4 * x for x in DBR_TE_std_reflect_design]
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
n_pairs = np.arange(2, 41, 1) + 0.5 # np.array([26.5])
# Set indices of refraction for input layer
index_Mat_in = index_Mat_AlGaAs
loss_Mat_in = loss_Mat_AlGaAs
# Set indices of refraction for alternating layers
index_Mat_1 = index_Mat_AlAs
loss_Mat_1 = loss_Mat_AlAs
index_Mat_2 = index_Mat_AlGaAs
loss_Mat_2 = loss_Mat_AlGaAs
# Set indices of refraction for output layer
index_Mat_out = index_Mat_GaAs # Air
loss_Mat_out = index_Mat_GaAs # Air

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

# Total tickness
total_thickness = thickness_ideal * (15 + 26.5) * 2 + wavelength_design / index_Mat_2_design
print(f"Total thcikness: {total_thickness}")

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
    this_DBR_TE_reflect_min = this_DBR_TE_reflect - 4 * this_DBR_TE_std_reflect
    this_DBR_TE_reflect_max = this_DBR_TE_reflect + 4 * this_DBR_TE_std_reflect
    print(f"N layers: {n_pair} | Reflectance: {this_DBR_TE_reflect} +- {this_DBR_TE_std_reflect} = [{this_DBR_TE_reflect_min} - {this_DBR_TE_reflect_max}]")
    # Next N
    n_pair_index += 1
    
# Reflectance versus number of layers
DBR_TE_error_reflect_design = [4 * x for x in DBR_TE_std_reflect_design]
plt.figure()
plt.errorbar(n_pairs, DBR_TE_mean_reflect_design, yerr=DBR_TE_error_reflect_design, fmt='k.', label='R_TE')
#plt.plot(n_pairs, DBR_TE_mean_reflect_design,'k.',label='R_TE')
plt.axhline(y = target_reflectance_bottom, linestyle = '--', color = 'k', label = 'Bottom DBR target reflectance')
plt.title('Bottom DBR mirro')
plt.xlabel('Number of alternating layers pairs')
plt.ylabel('Reflectance')
plt.ylim(0, 1.05)
plt.xlim(0, 42)
plt.grid()
plt.show()