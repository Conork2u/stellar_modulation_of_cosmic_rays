# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:47:39 2023

@author: conor
"""

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Imports ::::

import matplotlib.pyplot as plt
import numpy as np
import importlib
from scipy.integrate import solve_ivp
from scipy.integrate import quad
from astropy import units as u
from astropy import constants as const
#::| Throughout this code whenever the .value astropy fucntion is applied is 
#::| due to the calculation being unable to accept quantities with different 
#::| or no dimensions
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: Independant values ::::
print('----------------------------------')
print('Start of Stellar Values Calculation')
#::  Constants are placed into cgs to simplify subequent calculations
parsec_in_cm = const.pc.to('cm') # parsec in cm
print('parsec = ', f'{parsec_in_cm:.2E}')
AU = const.au.to('cm')
#print('AU = ', f'{AU:.2E}')
year_in_s = 1 * u.year.to('s') * u.s # year in seconds
#print('Year(cgs) = ', year_in_s)
k_B = 1.380649e-16 * (((u.cm**2 * u.g) / u.s**2) / u.K)
#print('Boltzmann Constant =',k_B) # Boltzmann Constant
n = 1 * u.cm**-3 
#print('Number Density =', n)
T = 1e4 * u.K 
#print('Temperature of Region II =', T)
m_p = 1.67262192e-30*u.gram  # Proton Mass in grams
#print('Proton Mass =', m_p)
solar_mass = const.M_sun.to('g')  #solar mass in grams
#print(f'Solar Mass = {solar_mass:.2E}')
solar_radius_cm = const.R_sun.to('cm') # Solar Radius in cm
print(f'Solar radius = {solar_radius_cm:.2E}')
solar_luminosity = (1 * u.Lsun).value # one SolLum
print(f'Solar luminosity = {solar_luminosity:.2E}')
teff_sun = 5780 # the effective teperature of sol
#::::::::::::::::::::::::::::::::::::::::::::::::::::::: Value Definitions ::::

#::| Values Dependant of the Variables associated with the star are inputted
#::| and converted to cgs, will be updated to include any F, G , K and M star
#::| For New entries, teff and luminosity(total bolometric) need to be inputted
#::| to calculate the HZ Limits and the model if only valid for F, G, K, and M stars

Stellar_Values = 'Sol'   # eg 'Sol' or 'Proxima_centauri'

if   Stellar_Values == 'Sol':
     mass_loss_rate = 6.804432124072167e-14 * solar_mass / year_in_s
     mass_loss_rate = mass_loss_rate.to('g/s')
     vw = 400 * u.km / u.s     
     vw = vw.to('cm/s')
     lw = 0.5 * mass_loss_rate * vw**2 
     lw = lw.to('erg/s') # Mechanical Luminosity
     stellar_age = 4.6e9 * year_in_s # stellar age in seconds  
     R_star = 1.0 * solar_radius_cm
     B0 = 1.0 # Gauss
     luminosity = solar_luminosity # SolLum
     teff = 5780 # Kelvin
elif Stellar_Values == 'Weaver':
     mass_loss_rate = 1e-6 *solar_mass / year_in_s
     mass_loss_rate = mass_loss_rate.to('g/s')
     vw = 2000. * u.km / u.s 
     vw = vw.to('cm/s')
     lw = 0.5 * mass_loss_rate * vw**2
     lw = lw.to('erg/s')
     R_star = 1.39E+12 *u.cm
elif Stellar_Values == 'Proxima Centauri': # estimations
     mass_loss_rate = 2e-13 * solar_mass / year_in_s
     mass_loss_rate = mass_loss_rate.to('g/s')
     vw = 50 * u.km / u.s     
     vw = vw.to('cm/s')
     lw = 0.5 * mass_loss_rate * vw**2 
     lw = lw.to('erg/s') # Mechanical Luminosity
     stellar_age = 4.6e9 * year_in_s # stellar age in seconds  
     R_star = 0.14 * solar_radius_cm
     B0 = 2000 # Gauss
     luminosity = solar_luminosity * 0.001567 # SolLum
     teff = 2992 # Kelvin
elif Stellar_Values == 'TRAPPIST-1': # estimations
     mass_loss_rate = 2e-15 * solar_mass / year_in_s
     mass_loss_rate = mass_loss_rate.to('g/s')
     vw = 50 * u.km / u.s     
     vw = vw.to('cm/s')
     lw = 0.5 * mass_loss_rate * vw**2 
     lw = lw.to('erg/s') # Mechanical Luminosity
     stellar_age = 4.6e9 * year_in_s # stellar age in seconds  
     R_star = 0.1192 * solar_radius_cm
     B0 = 600 # Gauss
     luminosity = solar_luminosity * 0.000553 # SolLum
     teff = 2566 # Kelvin     
     
else:
     raise ValueError('Invalid stellar value')
t6 = (10**6)
print(Stellar_Values, 'Mass loss rate =', f'{mass_loss_rate:.2E}')
print(Stellar_Values, 'Terminal Velocity (Vw) =', f'{vw:.2E}')
print(Stellar_Values, 'Mechanical luminosity (Lw) =',f'{lw:.2E}')
teff_value = teff
luminosity_value = luminosity

print('End of Stellar Values Calculation')
print('----------------------------------')














































