# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 18:23:21 2024

@author: user
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy.constants as scc

def Read_mat_file(wavelengths_interp, filename, make_plot = 0):
    materialinfo = pd.read_csv(f'{filename}.csv',sep=',',skiprows=0)
    # Divide wavelength, refractiveindex and extinctioncoeff into three matricies
    materialinfo = materialinfo.rename(columns={materialinfo.columns[0]:'wavelength',
                                                materialinfo.columns[1]:'refractiveindex',
                                                materialinfo.columns[2]:'extinctioncoeff',})
    wavelengths = (materialinfo['wavelength']*1e-6).to_numpy() # wavelength in meter
    refractiveindex = materialinfo['refractiveindex'].to_numpy()
    extinctioncoeff = materialinfo['extinctioncoeff'].to_numpy()
    # Interpolate in desired spectral bandwidth
    refractiveindex_interp = np.interp(wavelengths_interp, wavelengths, refractiveindex)
    extinctioncoeff_interp = np.interp(wavelengths_interp, wavelengths, extinctioncoeff)
    
    # Plot spectral reflectance for this design
    if make_plot:
        plt.figure()
        plt.plot(wavelengths*1e9, refractiveindex, 'b-', label='n')
        plt.plot(wavelengths_interp*1e9, refractiveindex_interp, 'g--',label='n interp')
        plt.plot(wavelengths*1e9, extinctioncoeff, 'k-', label='K')
        plt.plot(wavelengths_interp*1e9, extinctioncoeff_interp, 'r--',label='K interp')
        plt.legend()
        plt.show()
    
    return refractiveindex_interp, extinctioncoeff_interp

def Read_mat_file_Add_DSMOdel(wavelengths_interp, filename,  charge_density, effective_mass, damping, make_plot = 0):
    # Relative permitivity of undoped material
    refractiveindex_base, extinctioncoeff_base = Read_mat_file(wavelengths_interp, filename, 0)
    relative_permitivity_base = (refractiveindex_base + 1j * extinctioncoeff_base) ** 2
    # Contribution of free charge carriers to the relative permitivity
    relative_permitivity_fcc = Relative_permitivity_free_charge_carrier(wavelengths_interp, charge_density, effective_mass, damping)
    # Total relative permitivity
    # Experimental data is eps_r = 1 + eps_r,b. So eps_r,b = eps_r - 1.
    # Then we want to compute eps_r = eps_r,e + eps_r,b:
    relative_permitivity_total = relative_permitivity_fcc + relative_permitivity_base - 1
    # Resulting complex-valued index of refraction
    refractiveindex = np.sqrt( (np.abs(relative_permitivity_total) + np.real(relative_permitivity_total) ) / 2 )
    extinctioncoeff = np.sqrt( (np.abs(relative_permitivity_total) - np.real(relative_permitivity_total) ) / 2 )
    
    if make_plot:
        plt.figure()
        plt.plot(wavelengths_interp*1e9, np.real(relative_permitivity_fcc), 'k-', label='FCC Real')
        plt.plot(wavelengths_interp*1e9, np.imag(relative_permitivity_fcc), 'k--', label='FCC Imag')
        plt.plot(wavelengths_interp*1e9, np.real(relative_permitivity_base), 'r-', label='Baseline Real')
        plt.plot(wavelengths_interp*1e9, np.imag(relative_permitivity_base), 'r--', label='Baseline Imag')
        plt.grid()
        plt.xlabel('Wavelegnth (nm)')
        plt.ylabel('Real and Imaginary parts')
        plt.title('Drude-Summerfeld relative permitivity')
        plt.legend(loc='center left')
        plt.show()
    
    return refractiveindex, extinctioncoeff

def Relative_permitivity_free_charge_carrier(wavelengths, charge_density, effective_mass, damping):
    # Plasma frequency
    plasma_freq = np.sqrt(charge_density * scc.elementary_charge ** 2 / 
                          (effective_mass * scc.electron_mass * scc.epsilon_0))
    freq = np.divide(2 * np.pi * scc.speed_of_light, wavelengths)
    relative_permitivity_DSModel = 1 - (np.divide(plasma_freq ** 2, (freq ** 2 + (1j * freq * damping) )))
    #relative_permitivity_DSModel = ( 1 - np.divide(plasma_freq ** 2, freq ** 2 + damping ** 2) ) + (1j * np.divide(damping * plasma_freq ** 2, freq * (freq ** 2 + damping ** 2) ))
    return relative_permitivity_DSModel