# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 18:23:21 2024

@author: user
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

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