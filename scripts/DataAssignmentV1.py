import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np

#loading files for pete
path = '..\materials'
os.chdir(path)

#Loading files for sebastian
# path = 'materials'
# os.chdir(path)

for filename in [
                 "AlGaAs-X=0.452",
                 "AlGaAs-X=0.342",
                 "AlGaAs-X=0.097",
                 "GaAs"
                 ]:
    materialinfo = pd.read_csv(f'{filename}.csv',sep=',',skiprows=0)
    
   
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
    Target_Wavelength = 850

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
    
    filename = filename.replace(".csv", " ", 1)
    f'refractiveindex{filename}' == 12
    
    print(refractiveindex,extinctioncoeff,wavelength)
    print(f'refractiveindex{filename}')
    
