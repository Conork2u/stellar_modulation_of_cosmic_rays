
"""
Created on Thu Feb 22 13:08:12 2024

@author: conor
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import importlib
import sys
import math
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.integrate import quad
from astropy import units as u
from astropy import constants as const
import CRspectra as cr
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: Independant values ::::

#::  Constants are placed into cgs to simplify subequent calculations
parsec_in_cm = const.pc.to('cm') # parsec in cm
AU = const.au.to('cm')
year_in_s = 1 * u.year.to('s') * u.s # year in seconds
k_B = 1.380649e-16 * (((u.cm**2 * u.g) / u.s**2) / u.K)
n = 1 * u.cm**-3 
T = 1e4 * u.K 
m_p = 1.67262192e-24*u.gram  # Proton Mass in grams
solar_mass = const.M_sun.to('g')  #solar mass in grams
solar_radius_cm = const.R_sun.to('cm') # Solar Radius in cm
solar_luminosity = (1 * u.Lsun).value # one SolLum
teff_sun = 5780 # the effective teperature of sol
lightspeed = const.c.to('cm/s')
rho_ISM = n * m_p

#::::::::::::::::::::::::::::::::::::::::::::::::::::::: Dependant Values :::::
#SOL Test
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# plotting = False
# Spectral_Class = 'Sol'
# teff =  5780
# mass_loss_rate = 6.804432124072167e-14* solar_mass / year_in_s
# mass_loss_rate = mass_loss_rate.to('g/s')
# R_star =  1 * solar_radius_cm
# luminosity =  10**0 * solar_luminosity
# vw =  400* u.km / u.s
# vw = vw.to('cm/s')
# B0 =  1
# v_star =  22 * u.km / u.s
# v_star = v_star.to('cm/s')

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Spectral_Class = str(sys.argv[1])
teff =  float(sys.argv[2])
mass_loss_rate = float(sys.argv[3])* solar_mass / year_in_s
R_star =  float(sys.argv[4]) * solar_radius_cm
luminosity =  10**float(sys.argv[5])* solar_luminosity
vw =  float(sys.argv[6])* u.km / u.s
B0 =  float(sys.argv[7])
v_star =  float(sys.argv[8])* u.km / u.s
v_star = v_star.to('cm/s')
vw = vw.to('cm/s')
mass_loss_rate = mass_loss_rate.to('g/s')
plotting_input = sys.argv[9]

if plotting_input.lower() == "true":
    plotting = True
elif plotting_input.lower() == "false":
    plotting = False
else:
    # Handle invalid input, for example, setting a default value
    print("Invalid input. Defaulting to False.")
    plotting = False


#  check if oit one of the three values R star L or T eff
#teff = teff_value
#luminosity = luminosity_value
Stellar_Values = Spectral_Class
lw = 0.5 * mass_loss_rate * vw**2 
lw = lw.to('erg/s') # Mechanical Luminosity
R_star_pc = R_star / parsec_in_cm
R_star_au = R_star / AU


#::::::::::::::::::::::::::::::::::::::::::::::::: Hab Zone Calcualtions :::::
"""
 "Habitable Zones Around Main-Sequence Stars: New Estimates" by Kopparapu et al.(2013), Astrophysical Journal, 765, 131    arXiv link

 "Habitable Zones Around Main-Sequence Stars: Dependence on Planetary Mass" by Kopparapu et al.(2014), Astrophysical Journal Letters, 787, L29    arXiv
  link
  
 "Translated to python code by John Armstrong (jcarmstrong@weber.edu) 04 June 2014"
"""
#************************************************************************************
# This code calculates habitable zone 'fluxes' using the expression given in the 
# Kopparapu et al.(2014) paper. The corresponding output file is 'HZs.dat'. 
# It also generates a file 'HZ_coefficients.dat' that gives the coefficients for 
# the analytical expression.
#
# Ravi kumar Kopparapu April 19 2014
#
# 
#************************************************************************************
#************************************************************************************
# Output files.

print('Start of HZ Calculations')


if not os.path.isfile('./HZ_coefficients.dat'): # and not os.path.isfile('./HZs.dat'):
    #hzdat = open('HZs.dat', 'w')
    hzcoeff = open('HZ_coefficients.dat', 'w')
    
    #************************************************************************************
    # Coeffcients to be used in the analytical expression to calculate habitable zone flux 
    # boundaries
    
    seff = [0,0,0,0,0,0]
    seffsun  = [1.776,1.107, 0.356, 0.320, 1.188, 0.99] 
    a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
    b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
    c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
    d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]    
    #************************************************************************************
    # Writing coefficients into 'HZ_coefficients.dat' file

    hzcoeff.write('# The coefficients are as follows. The columns, i, are arranged according to\n')
    hzcoeff.write('# the HZ limits given in the paper.\n')
    hzcoeff.write('#\n')
    hzcoeff.write('# i = 1 --> Recent Venus\n')
    hzcoeff.write('# i = 2 --> Runaway Greenhouse\n')
    hzcoeff.write('# i = 3 --> Maximum Greenhouse\n')
    hzcoeff.write('# i = 4 --> Early Mars\n')
    hzcoeff.write('# i = 5 --> Runaway Greenhouse for 5 ME\n')
    hzcoeff.write('# i = 6 --> Runaway Greenhouse for 0.1 ME\n')
    
    hzcoeff.write('# First row: S_effSun(i)\n')
    hzcoeff.write('# Second row: a(i)\n')
    hzcoeff.write('# Third row:  b(i)\n')
    hzcoeff.write('# Fourth row: c(i)\n')
    hzcoeff.write('# Fifth row:  d(i)\n')
    hzcoeff.write('   '+'{:6.6E}'.format(seffsun[0]) + '  ' +
                  '{:6.6E}'.format(seffsun[1]) + '  ' +
                  '{:6.6E}'.format(seffsun[2]) + '  ' +
                  '{:6.6E}'.format(seffsun[3]) + '  ' +
                  '{:6.6E}'.format(seffsun[4]) + '  ' +
                  '{:6.6E}'.format(seffsun[5]) + '  ' +
                  '\n')
    hzcoeff.write('   '+'{:6.6E}'.format(a[0]) + '  ' +
                  '{:6.6E}'.format(a[1]) + '  ' +
                  '{:6.6E}'.format(a[2]) + '  ' +
                  '{:6.6E}'.format(a[3]) + '  ' +
                  '{:6.6E}'.format(a[4]) + '  ' +
                  '{:6.6E}'.format(a[5]) + '  ' +
                  '\n')
    hzcoeff.write('   '+'{:6.6E}'.format(b[0]) + '  ' +
                  '{:6.6E}'.format(b[1]) + '  ' +
                  '{:6.6E}'.format(b[2]) + '  ' +
                  '{:6.6E}'.format(b[3]) + '  ' +
                  '{:6.6E}'.format(b[4]) + '  ' +
                  '{:6.6E}'.format(b[5]) + '  ' +
                  '\n')
    hzcoeff.write('   '+'{:6.5e}'.format(c[0]) + '  ' +
                  '{:6.5e}'.format(c[1]) + '  ' +
                  '{:6.5e}'.format(c[2]) + '  ' +
                  '{:6.5e}'.format(c[3]) + '  ' +
                  '{:6.5e}'.format(c[4]) + '  ' +
                  '{:6.5e}'.format(c[5]) + '  ' +
                  '\n')
    hzcoeff.write('   '+'{:6.5e}'.format(d[0]) + '  ' +
                  '{:6.5e}'.format(d[1]) + '  ' +
                  '{:6.5e}'.format(d[2]) + '  ' +
                  '{:6.5e}'.format(d[3]) + '  ' +
                  '{:6.5e}'.format(d[4]) + '  ' +
                  '{:6.5e}'.format(d[5]) + '  ' +
                  '\n')

    #************************************************************************************
    # Calculating HZ fluxes for stars with 2600 K < T_eff < 7200 K. The output file is
    # 'HZ_fluxes.dat'
   # sys.path.append('C:/Users/conor/Stellar Bubbles/Final Model')

    # L = luminosity_value 
    # teff = teff_value
    # print('Effective Temperature =', teff)
    # print('Stellar Luminosity =',  L)
    # hzdat.write('#  Teff(K)        Recent        Runaway        Maximum        Early        5ME Runaway   0.1ME Runaway\n')
    # hzdat.write('#                 Venus         Greenhouse     Greenhouse     Mars         Greenhouse    Greenhouse\n')
    
    # starTemp = []
    # recentVenus = []
    # runawayGreenhouse = []
    # maxGreenhouse = []
    # earlyMars = []
    # fivemeRunaway = []
    # tenthmeRunaway = []
    
    # while(teff <= 7201.0): 
    #   tstar = teff - 5780.0
    #   for i in range(len(a)):
    #      seff[i] = seffsun[i] + a[i]*tstar + b[i]*tstar**2 + c[i]*tstar**3 + d[i]*tstar**4
    
    #   starTemp.append(teff)
    #   recentVenus.append(seff[0])
    #   runawayGreenhouse.append(seff[1])
    #   maxGreenhouse.append(seff[2])
    #   earlyMars.append(seff[3])
    #   fivemeRunaway.append(seff[4])
    #   tenthmeRunaway.append(seff[5])
         
    #   hzdat.write('   '+'{:6.0f}'.format(teff) + '      ' +
    #               '{:6.6E}'.format(seff[0]) + '      ' +
    #               '{:6.6E}'.format(seff[1]) + '     ' +
    #               '{:6.6E}'.format(seff[2]) + '   ' +
    #               '{:6.6E}'.format(seff[3]) + ' ' +
    #               '{:6.6E}'.format(seff[4]) + ' ' +
    #               '{:6.6E}'.format(seff[5]) + '  ' +
    #               '\n')
    #   teff = teff + 2
    '''
    print ('************************************************************')
    print ('')
    print ('The HZ coefficients are printed in HZ_coefficients.dat file.')
    print ('HZs for stars with 2600 K <= Teff <=7200 K is in HZs.dat file.')
    print ('')
    print ('************************************************************')
    '''
    #hzdat.close()
    hzcoeff.close()

# sys.exit(0)
#teff = teff_value
def compute_hz_boundaries(teff, L):
    Tssun = 5780 # Sun's effective temperature in K
    tstar = teff - Tssun
    
    # Reading coefficients from .dat file
    coefficients = np.loadtxt("HZ_coefficients.dat")

    # Map coefficients to each habitable zone category
    categories = [
        'recent_venus',
        'runaway_greenhouse',
        'max_greenhouse',
        'early_mars',
        'runaway_greenhouse_5ME',
        'runaway_greenhouse_0.1ME'
    ]

    # Calculate Seff and distance for each HZ category
    hz_boundaries = {}
    for i, category in enumerate(categories):
        seff = sum(coefficients[j][i] * tstar**j for j in range(len(coefficients)))
        distance = math.sqrt(L / seff)
        hz_boundaries[category] = distance

    return hz_boundaries


# Call the function and print the boundaries
boundaries = compute_hz_boundaries(teff,luminosity)
# for category, distance in boundaries.items():
#     print(f"{category} = {distance} AU")
# Convert boundaries from AU to parsec
hz_boundaries_pc = {category: (distance * u.AU).to(u.pc).value for category, distance in boundaries.items()}
#convert boundaries from AU to cm
hz_boundaries_cm = {category: (distance) * AU for category, distance in boundaries.items()}

hz_boundaries_au = {category: (distance)  for category, distance in boundaries.items()}
# Print each category with its corresponding distance in parsecs


#for category, distance in hz_boundaries_pc.items():
    #print(f"{category} = {distance} parsecs" )

#for category, distance in hz_boundaries_cm.items():
    #print(f"{category} = {distance:.2E}")


# Retrieve 'recent_venus' and 'early_mars' values from hz_boundaries_cm
recent_venus_cm = hz_boundaries_cm['recent_venus']
early_mars_cm = hz_boundaries_cm['early_mars']
# Calculate the mean distance in cm

hz_mean_cm = (recent_venus_cm + early_mars_cm) / 2

Mid_HZ = (hz_boundaries_au['recent_venus'] + hz_boundaries_au['early_mars']) /2 
print('Recent Venus=',hz_boundaries_au['recent_venus'])
print('Mid_HZ=',Mid_HZ)
print('Early Mars=',hz_boundaries_au['early_mars'])

print('End of HZ Calculations')
print('----------------------------------') 
  
#::::::::::::::::::::::::::::::::::::::::::::::::: Stellar profiles ::::::::::

#::|Pressure Calcualtions

def pRAM(rho_ISM, v_star):
    return (rho_ISM * v_star**2)
def pthermal(n, k_B, T):
    return(T * k_B * n)
def pISM(pthermal, pRAM):
    return (pthermal(n, k_B, T) + pRAM(rho_ISM, v_star))   

#::| The Radius of the termination shock (R1) and the contact discontiniuty(RC)
#::| are calculated and converted in to values suitable for scaling the plot
#::| Taking into account that the star is moving, the equation for R1 Changes
def R1(lw, vw, pISM):
    return np.sqrt(lw / (vw * (2 * np.pi) * pISM(pthermal, pRAM))).value

R1_Value = R1(lw, vw, pISM) * u.cm
R1_pc_Value = R1_Value.to('pc')
R1_au_Value = R1_Value.to('AU')
R1_Valueplus = R1_Value + 1*u.cm
R1_Valueminus = R1_Value - 1*u.cm
print('R1 = ',R1_au_Value)
 
#::|RC Contact discontinuity calcualtions
RC = (1+0.38) * R1_Value
RC_Value = (1+0.38) * R1_Value 
RC_pc_Value = RC_Value.to('pc')
RC_au_Value = RC_Value.to('au')
RC_Valueplus =  RC_Value + 0.1*u.cm
RC_Valueminus = RC_Value - 10000*u.cm

#::| The velocity profile is defined, always taking cm as the input unit
def v(radius):
    Zone_A = np.where(radius < R1_Value, vw, 0.) 
    Zone_B = np.where((radius >= R1_Value)*(radius <= RC_Value), (vw/4) * (R1_Value**2/radius**2).value, Zone_A)
    Zone_C = np.where(radius > RC_Value, 0. ,Zone_B )
    return Zone_C

#::| The Magnetic Field profile is defined
def B(radius):
    # 
    Zone_A = np.where((0 < radius) * (radius < R_star) , 0.,  0.)
    Zone_B = np.where((R_star <= radius) * (radius <= R1_Value), B0 * R_star / radius, Zone_A)
    Zone_C = np.where((R1_Value < radius)* (radius <= RC_Value), B0*4 * R_star / R1_Value, Zone_B)
    Zone_D = np.where(radius > RC_Value, 0., Zone_C)
    return Zone_D


#::| The Cosmic Ray Profile is defined
#::|Variables for the comsic ray function ::
c, d = -0.5, 0.5
N_0 = 1 #   initial cosmic ray density 
B_ISM = 5**-6 #Gauss
D_0 = 4e23 # cm^2/s ## Ideally ~~ 10 ^ 26  #4.0
E_0 = 0.01602 #10 Gev = 0.01602 ergs in cgs

def N_B(radius,E): 
    beta = (v(R1_Valueplus).value * RC_Value.value) / (D_0 * (B(RC_Valueminus)/B_ISM)**c    * (E/E_0)**d).value
    a , b = -2 , 0
    C = 1
    delta = (a - (b * c) +1)
    return N_0        *        np.exp(beta  * (-1 * (1./(delta) * (RC_Value.value/radius.value) ** delta) + (1./(delta) * (C) ** delta)))

def N_A(radius,E):
    alpha =  (vw.value  * R1_Value.value )    /      (D_0 * (B(R1_Valueminus)/B_ISM)**c   * (E/E_0)**d).value
    a , b = 0 , -1
    C = 1
    delta = (a - (b * c) +1) 
    return N_B(R1_Valueplus,E) * np.exp(alpha * (1 *(  1./(delta) * (radius.value /R1_Value.value) ** delta) - (1./(delta) * (C) ** delta)))


def N(radius,E): 
    Zone_ISM  = np.where(radius > RC_Value, 1.0, 1.0)
    Zone_B    = np.where((R1_Value < radius) * (radius < RC_Value), N_B(radius,E), Zone_ISM)
    Zone_A    = np.where((R_star < radius) * (radius < R1_Value), N_A(radius,E), Zone_B)
    Zone_Star = np.where((0 < radius) * (radius <= R_star), 0.0, Zone_A)
    return Zone_Star

#::| The Diffusion Coefficent Profile is defined

def D_c_B(radius, E):
    a , b = -2 , 0
    return  D_0 * ((B(RC_Valueminus)/B_ISM) * (radius.value / RC_Valueminus.value)**b)**c * (E_0 / E).value**d

def D_c_A(radius, E):
    a , b = 0 , -1
    return  D_0 * ((B(R1_Valueminus)/B_ISM) * (radius.value / R1_Valueminus.value)**b)**c * (E_0 / E).value**d
def D_c(radius, E):
    Zone_ISM  = np.where(radius > RC_Value, 1.0, 1.0)
    Zone_B    = np.where((R1_Value < radius) & (radius <= RC_Value), D_c_B(radius, E), 0.0)
    Zone_A    = np.where((R_star < radius) & (radius <= R1_Value), D_c_A(radius, E), Zone_B)
    Zone_Star = np.where(radius <= R_star, 0.0, Zone_A)
    return Zone_Star



#::::::::::::::::::::::::::::::::::::::::::::::::: Plotting ::::::::::::::::::
    
    

#::| Plotting parameters
E_min = 0.001602 / 100 * 1.6
E_max = 0.00162  * 1000
#10 Gev = 0.01602 ergs in cgs

energy = np.geomspace(E_min, E_max, 100) *u.erg
radius = np.geomspace(R_star.value, 1.3*RC_Value.value, 100000) * u.cm 
# radius = np.array([rhz_in rhz med rhz max])
radenergy = np.meshgrid(radius, energy)
radius = radenergy[0]
energy = radenergy[1]






if plotting:
#:::::::::::::::::::::: Wind Velocity Profile Plotting ::::::::::::::::::::::::

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    ax.plot((radius/AU)[0,:], v(radius[0,:]).to('km/s'), label="Velocity Profile")
    ax.axvline(x=R1_au_Value.value, color='purple', linestyle=':', label='Wind Termination Shock') 
    ax.axvline(x=RC_au_Value.value, color='green' , linestyle=':', label='Contact Discontinuity')
    ax.axvline(x=R_star_au.value, color='Black' , linestyle='-', label='Stellar Radius')
    ax.axhline(y=22, color='red', linestyle='--', label='Stellar Velocity ')
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['recent_venus'], hz_boundaries_au['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['runaway_greenhouse'], hz_boundaries_au['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
    #::| The x axis is set to a logarithmic scale as the difference between R1 and
    #::| RC is more than a factor of one thousand
    #scaling
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    #labels
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('Velocity[km/s]')
    ax.set_title('Wind Velocity profile for '+ Stellar_Values) 
    ax.legend()
    ax.grid(False)
    plt.show()


#:::::::::::::::::::::: Magnetic Field Profile Plotting :::::::::::::::::::::::

    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
    ax.plot((radius/AU)[0,:], B(radius[0,:]), label=" Magnetic Field Profile")
    ax.axvline(x=R1_au_Value.value, color='purple', linestyle=':', label='Wind Termination Shock') 
    ax.axvline(x=RC_au_Value.value, color='green' , linestyle=':', label='Contact Discontinuity')
    ax.axvline(x=R_star_au.value, color='Black' , linestyle='-', label='Stellar Radius')
    
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['recent_venus'], hz_boundaries_au['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['runaway_greenhouse'], hz_boundaries_au['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
    #scaling
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    #labels
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('Magentic Field (Gauss)')
    ax.set_title(' Magnetic Field Profile for '+ Stellar_Values )
    ax.legend()
    ax.grid(False)                                                                                                                                                                                                                                                           
    plt.show()

#::::::::::::::::::::::::: Cosmic Ray Profile Plotting ::::::::::::::::::::::::

    e_index = 70
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
    NrE =N(radius,energy)
    
    ax.plot((radius/AU)[e_index,:],NrE[e_index,:], label=" Cosmic Ray Profile")
    ax.axvline(x=R1_au_Value.value, color='purple', linestyle=':', label='Wind Termination Shock') 
    ax.axvline(x=RC_au_Value.value, color='green' , linestyle=':', label='Contact Discontinuity')
    ax.axvline(x=R_star_au.value, color='Black' , linestyle='-', label='Stellar Radius')
    
    # Hab Zone
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['recent_venus'], hz_boundaries_au['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['runaway_greenhouse'], hz_boundaries_au['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
    #scaling
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    
    #labels
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('N(r,E)')
    ax.set_title(' Cosmic Ray Profile for '+ Stellar_Values )
    ax.legend()
    ax.grid(False)
    plt.show()

#::::::::::::::::::::::: Mean Free Path plotting ::::::::::::::::::::::::::::::

    MFP = ((3 * D_c(radius, energy)) / lightspeed.value)/AU
    
    def read_and_validate_csv(file_path):
         try:
             data = pd.read_csv(file_path)
             if 'radius' not in data.columns or 'mean_free_path' not in data.columns:
                 raise ValueError("CSV file must contain 'radius' and 'mean_free_path' columns")
             return data
         except Exception as e:
             print(f"Error reading {file_path}: {e}")
             sys.exit(1)
    
    def plot_data(ax, data, label):
         ax.plot(data['radius'], data['mean_free_path'], label=label)
    
     #Read and validate the CSV data
    data1 = read_and_validate_csv(r"C:\Users\conor\Stellar Bubbles\Parallel MFP.csv")
    data2 = read_and_validate_csv(r"C:\Users\conor\Stellar Bubbles\Radial MFP.csv")
    
    e_index = 0
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
    
     #Plot the mean free path as a function of radius (in AU)
    ax.plot((radius / AU)[e_index, :], MFP[e_index, :], label=r'Parallel Mean Free Path $\lambda(r,E)$')
    plot_data(ax, data1, 'Parallel MFP Pei et al')
    plot_data(ax, data2, 'Radial MFP Pei et al')
    
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('Mean Free Path [AU]')
    ax.set_title('Mean Free Path '+ '(Stationary)')
    
    ax.grid(False)
    ax.legend()
    plt.show()

#::::::::::::::::::::::::: Cosmic Ray Spectra  Plotting :::::::::::::::::::::::

    r_index = np.argmin(abs(radius[0,:].value - (hz_mean_cm.value)))
    
    
    # Calculate residuals
    residuals = NrE[:, r_index] - cr.CR_modulation(energy[:, r_index].to(u.eV).value)
    
    # Remove plotting fucntions 
    # write NrE to file for each value of R_hab(min med Max),   assign unique name to values in text based on 
    # SpV B_field valiue and V-star Value 
    
    # Create a plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), dpi=600,sharex = True , gridspec_kw={'height_ratios': [3, 1], 'hspace':0})
    
    # Plot the cosmic ray spectra 
    ax1.plot(energy[:, r_index].to(u.GeV), NrE[:, r_index], label="Cosmic Ray Spectra")
    ax1.plot(energy[:, r_index].to(u.GeV), cr.CR_modulation(energy[:, r_index].to(u.eV).value), label='Rogers-lee 2020')
    
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('GeV')
    ax1.set_ylabel('N')
    ax1.set_title(f'Cosmic Ray Spectra for {Stellar_Values} at $D_{{\mathrm{{ref}}}} = {D_0:.1e}$')
    ax1.legend()
    ax1.grid(True)
    
    # Plot the residuals on the second subplot
    ax2.plot(energy[:, r_index].to(u.GeV), residuals, label="Residuals", color='red')
    ax2.set_xscale('log')
    ax2.set_xlabel('GeV')
    ax2.set_ylabel('Residual')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    
    plt.show()
else:
    print('No Plots')