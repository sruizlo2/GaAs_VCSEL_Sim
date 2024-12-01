# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:59:09 2024

@author: srulo
"""

import numpy as np

def WriteGrid(file, length, n_pints):
    file.write(f"grid length={length} points={n_pints:.0f}\n")
    return

def WriteRegion(file, length, region_type):
    file.write(f"region {region_type} length={length}\n")
    return

def WriteLayer(file, thickness, refractive_index, concentration):
    if concentration > 0:
        if refractive_index > 0:
            file.write(f"Refractive_Index length={thickness} value={refractive_index}\n")
        file.write(f"structure material=gaas alloy=al length={thickness} conc={concentration}\n")
    else:
        if refractive_index > 0:
            file.write(f"Refractive_Index length={thickness} value={refractive_index}\n")
        file.write(f"structure material=gaas length={thickness}\n")
    return

def WriteDoping(file, thickness, donors_conc, acceptors_conc):
    file.write(f"doping length={thickness}")
    if donors_conc > 0:
        file.write(f" Nd={donors_conc}")
    if acceptors_conc > 0:
        file.write(f" Na={acceptors_conc}")
    file.write("\n")
    return

class Layer:
  def __init__(self, thickness, concentrations, refractive_indices):
    self.thickness = thickness
    self.concentrations = concentrations
    self.refractive_indices = refractive_indices

def WreiteDeviceFile(filename,points_per_nm, thicknesses, concentrations, refractive_indices, n_pairs, charge_densities):
    total_thickness = thicknesses(0)
    ideal_layer_thickness_top = thicknesses(1)
    ideal_thickness_top = thicknesses(2)
    ideal_thickness_bulk_top = thicknesses(3)
    
    index_AlGaAs_p_doped = thicknesses(4)
    index_AlGaAs_p_doped = thicknesses(4)
    
    n_pairs_top = n_pairs(0)
    n_pairs_bottom = n_pairs(1)
    
    charge_density_e = charge_densities(0)
    charge_density_h = charge_densities(1)
    
    Al_concentration = concentrations(2)
    Al_concentration = concentrations(2)
    
    with open(f"{filename}.DEV", "w") as file:
        # Create device grid
        WriteGrid(file, total_thickness, total_thickness * 1e3 * points_per_nm)
        file.write("\n")
        # Start with bulk region
        WriteRegion(file, ideal_thickness_bulk_top, "bulk")
        # Top DBR
        for i in range(n_pairs_top):
            WriteLayer(file, ideal_layer_thickness_top, index_AlGaAs_p_doped, Al_concentration)
            WriteLayer(file, ideal_layer_thickness_top, index_AlAs_p_doped, 1)
        
        file.write("\n")
        # Top confinement layer
        WriteLayer(file, ideal_thickness_confinement, index_AlGaAs_undoped, Al_concentration)
        file.write("\n")
        # Quantum well (change to qw region)
        WriteRegion(file, quantum_well_thickness, "qw")
        WriteLayer(file, quantum_well_thickness, index_GaAs_undoped, 0)
        file.write("\n")
        # Continue with bulk region
        WriteRegion(file, ideal_thickness_bulk_bottom, "bulk")
        # Bottom confinement layer
        WriteLayer(file, ideal_thickness_confinement, index_AlGaAs_undoped, Al_concentration)
        file.write("\n")
        # Bottom DBR
        for i in range(int(np.ceil(n_pairs_bottom))):
            WriteLayer(file, ideal_layer_thickness_bottom, index_AlAs_n_doped, 1)
            if i <= n_pairs - 1:
                WriteLayer(file, ideal_layer_thickness_bottom, index_AlGaAs_n_doped, Al_concentration)
        file.write("\n")
        # Subtrate
        WriteLayer(file, substrate_thickness, index_GaAs_n_doped, 0)
        # Doping
        file.write("\n")
        # Top mirror is p-type
        WriteDoping(file, ideal_thickness_top, charge_density_h * 1e-6, 0)
        # Cavity is undoped
        WriteDoping(file, ideal_thcikness_cavity, 0, 0)
        # Bottom mirror is n-type
        WriteDoping(file, ideal_thickness_bottom + substrate_thickness, 0, charge_density_e * 1e-6)