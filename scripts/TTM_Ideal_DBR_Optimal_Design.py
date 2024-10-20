# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 18:14:02 2024

@author: srulo
"""

import numpy as np
from matplotlib import pyplot as plt
from TTM import Transfer_matrix_system as TM
from TTM import Trans_and_reflec as TR
from TTM import Transmitted_Angle as TA
from MatFileLoader import Read_mat_file as ReadMatFile
from CreateDBR import Create_DBR as DBR
import os

# Define spectrum
wavelength_design = 850e-9; # meters
wavelength_range = 10e-9; # meters
wavelegnth_step = 1e-9; # meters
# Vector with testing wavelengths
wavelengths = np.arange(wavelength_design - wavelength_range / 2,
                        wavelength_design + wavelength_range / 2, wavelegnth_step)
wavelength_design_index = np.argwhere(abs(wavelengths - wavelength_design) < 1e-16)[0][0]

# Path to materials files
path_to_mat_files = '../materials'
mat_file_GaAs = 'GaAs'
mat_file_AlGaAs = 'AlGaAs-X=0.452'
mat_file_AlAs = 'AlAs'

# Design parameters
target_reflectance_top = 0.995;
target_reflectance_bottom = 0.95;
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

# Orthogonal incidence
angle_inc = 0

# Vector with testing number of layer pairs
n_pairs = np.arange(1, 50, 5)

# Initialize variables
n_pair_index = 0;
DBR_TE_trans = np.zeros((n_pairs.shape[0], wavelengths.shape[0]))
DBR_TE_reflect = np.zeros((n_pairs.shape[0], wavelengths.shape[0]))
DBR_TM_trans = np.zeros((n_pairs.shape[0], wavelengths.shape[0]))
DBR_TM_reflect = np.zeros((n_pairs.shape[0], wavelengths.shape[0]))
DBR_TE_trans_design = []
DBR_TE_reflect_design = []
DBR_TM_trans_design = []
DBR_TM_reflect_design = []

# Iterate over design configurations
for n_pair in n_pairs:
    wavelength_index = 0;
    for wavelength in wavelengths:
        # Assing index of refractions and losses for this wavelength
        index_in = index_Mat_in[wavelength_index] + 1j * loss_Mat_in[wavelength_index]
        index_1 = index_Mat_1[wavelength_index] + 1j * loss_Mat_1[wavelength_index]
        index_2 = index_Mat_2[wavelength_index] + 1j * loss_Mat_2[wavelength_index]
        index_out = index_Mat_out[wavelength_index] + 1j * loss_Mat_out[wavelength_index]
        # Average index of refraction
        #average_index = (index_Mat_1[wavelength_index] + index_Mat_2[wavelength_index]) / 2
        #thickness_ideal = wavelength_design / 4 / average_index # lambda / 4 condition (per layer)
        # Create DBR
        DBR_ = DBR(n_pair, index_in, index_1, index_2, thickness_ideal, index_out)
        # Transfer matrix (wavelegnth in meters!)
        transfer_matrix_TE_DBR, transfer_matrix_TM_DBR = TM(DBR_, wavelength, angle_inc)
        # Transmittance and Reflectance
        DBR_TE_trans_, DBR_TE_reflect_ = TR(transfer_matrix_TE_DBR)
        DBR_TM_trans_, DBR_TM_reflect_ = TR(transfer_matrix_TM_DBR)
        DBR_TE_trans[n_pair_index, wavelength_index] = DBR_TE_trans_
        DBR_TE_reflect[n_pair_index, wavelength_index] = DBR_TE_reflect_
        DBR_TM_trans[n_pair_index, wavelength_index] = DBR_TM_trans_
        DBR_TM_reflect[n_pair_index, wavelength_index] = DBR_TM_reflect_
        wavelength_index += 1
    
    # Assign reflectance for design wavelength
    DBR_TE_trans_design.append(DBR_TE_trans[n_pair_index, wavelength_design_index])
    DBR_TE_reflect_design.append(DBR_TE_reflect[n_pair_index, wavelength_design_index])
    DBR_TM_trans_design.append(DBR_TM_trans[n_pair_index, wavelength_design_index])
    DBR_TM_reflect_design.append(DBR_TM_reflect[n_pair_index, wavelength_design_index])    
    
    # Plot spectral reflectance for this design
    show_spectral_reflectance = False
    if show_spectral_reflectance:
        plt.figure()
        plt.plot(wavelengths*1e9, DBR_TE_reflect[n_pair_index, :],'k',label='R_TE')
        plt.plot(wavelengths*1e9, DBR_TM_reflect[n_pair_index, :],'b--',label='R_TM')
        plt.axvline(x = wavelength_design*1e9, linestyle = '--', color = 'k', label = 'Design wavelength')
        plt.axhline(y = target_reflectance_top, linestyle = '--', color = 'k', label = 'Top DBR target reflectance')
        plt.axhline(y = target_reflectance_bottom, linestyle = '--', color = 'k', label = 'Bottom DBR Target reflectance')
        plt.xlim(wavelengths[0]*1e9, wavelengths[-1]*1e9)
        plt.ylim(0, 1)
        plt.show()
    
    # Next N
    n_pair_index += 1
    
# Reflectance versus number of layers
plt.figure()
plt.plot(n_pairs, DBR_TE_reflect_design,'k.',label='R_TE')
plt.axhline(y = target_reflectance_top, linestyle = '--', color = 'k', label = 'Top DBR target reflectance')
plt.axhline(y = target_reflectance_bottom, linestyle = '--', color = 'k', label = 'Bottom DBR Target reflectance')
plt.show