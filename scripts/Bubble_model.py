
"""
Created on Thu Feb 22 13:08:12 2024

@author: conor
"""
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
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
plt.close('all')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: Independant values ::::
print("START SCRIPT")
#::  Constants are placed into cgs to simplify subequent calculations
parsec_in_cm = const.pc.to('cm') # parsec in cm
AU = const.au.to('cm')
year_in_s = 1 * u.year.to('s') * u.s # year in seconds
k_B = 1.380649e-16 * (((u.cm**2 * u.g) / u.s**2) / u.K)


# ISM
n = 1 * u.cm**-3 
T = 1e4 * u.K
"""
# LISM
n = 0.2 * u.cm**-3 
T = 7e3 * u.K
"""


m_p = 1.67262192e-24*u.gram  # Proton Mass in grams
solar_mass = const.M_sun.to('g')  #solar mass in grams
solar_radius_cm = const.R_sun.to('cm') # Solar Radius in cm
solar_luminosity = (1 * u.Lsun).value # one SolLum
teff_sun = 5780 # the effective teperature of sol
lightspeed = const.c.to('cm/s')

rho_ISM = n * m_p
# Corrected Rho for helium
#mu = 1.3 # Helium micture
#rho_ISM = n * m_p * mu


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::: Dependant Values :::::
# SOL Test
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

plotting = True
Spectral_Class = 'Sol'
teff =  5780
mass_loss_rate = 6.8E-14* solar_mass / year_in_s
mass_loss_rate = mass_loss_rate.to('g/s')
print(f'{mass_loss_rate:.2E}')
R_star =  1 * solar_radius_cm
luminosity =  10**0* solar_luminosity
vw =  400* u.km / u.s
vw = vw.to('cm/s')
B0 =  1
v_star_km =  22 * u.km / u.s
v_star = v_star_km.to('cm/s')


# F5V Test
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''
plotting = True
Spectral_Class = 'F5V'
teff =  6550
mass_loss_rate = 7.94E-13* solar_mass / year_in_s
mass_loss_rate = mass_loss_rate.to('g/s')
R_star = 1.473 * solar_radius_cm
luminosity =  10**0.56 * solar_luminosity
vw =  800* u.km / u.s
vw = vw.to('cm/s')
B0 =  50
v_star_km =  2 * u.km / u.s
v_star = v_star_km.to('cm/s')
'''
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''
#Simcode Run
Spectral_Class = str(sys.argv[1])
teff =  float(sys.argv[2])
mass_loss_rate = float(sys.argv[3])* solar_mass / year_in_s
R_star =  float(sys.argv[4]) * solar_radius_cm
luminosity =  10**float(sys.argv[5])* solar_luminosity
vw =  float(sys.argv[6])* u.km / u.s
B0 =  float(sys.argv[7])
v_star_km =  float(sys.argv[8])* u.km / u.s
v_star = v_star_km.to('cm/s')
vw = vw.to('cm/s')
mass_loss_rate = mass_loss_rate.to('g/s')
plotting = False #sys.argv[9]
'''
if len(sys.argv) >= 9:
    Spectral_Class = str(sys.argv[1])
    teff = float(sys.argv[2])
    mass_loss_rate = float(sys.argv[3]) * solar_mass / year_in_s
    R_star = float(sys.argv[4]) * solar_radius_cm
    luminosity = 10**float(sys.argv[5]) * solar_luminosity
    vw = float(sys.argv[6]) * u.km / u.s
    B0 = float(sys.argv[7])
    v_star_km = float(sys.argv[8]) * u.km / u.s
    v_star = v_star_km.to('cm/s')
    vw = vw.to('cm/s')
    mass_loss_rate = mass_loss_rate.to('g/s')
    plotting = False

##::::::::::::::::::::::::::::::::Teegarden's Star
'''
plotting = True
Spectral_Class = 'TGN'  # Teegarden's Star '
teff =  2680
mass_loss_rate = 6.8E-15* solar_mass / year_in_s
mass_loss_rate = mass_loss_rate.to('g/s')
R_star =  0.12 * solar_radius_cm
luminosity =  10**-3.19* solar_luminosity
vw =  80* u.km / u.s
vw = vw.to('cm/s')
B0 =  1200
v_star_km =  92 * u.km / u.s
v_star = v_star_km.to('cm/s')
'''
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
print('RC = ',RC_au_Value)
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
D_0 = 3.8e23 # cm^2/s ## Ideally ~~ 10 ^ 26  #4.0
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


e_conversion = const.e.value * 10**16
#::| Plotting parameters
E_min = e_conversion / 100 * 1.59
E_max = e_conversion  * 1000
pei_value = (100 * u.MeV).to(u.erg)
 # If running for single anaylsis set to 100, if running for multiple set to 500
# index_to_insert = np.searchsorted(energy, pei_value)
# energy = np.insert(energy, index_to_insert, pei_value)
# 10 Gev = 0.01602 ergs in cgs
print ('emin',E_min)
#E_minmev = E_min.to('MeV')
#print ('emin MeV',E_minmev)
print(E_max)

if plotting:
    radius = np.geomspace(R_star.value, 1.3*RC_Value.value, 100000) * u.cm # for single star plotting
    energy = np.geomspace(E_min, E_max, 100) *u.erg 
else:
    radius = np.array([recent_venus_cm.value,  hz_mean_cm.value,  early_mars_cm.value]) * u.cm 
    energy = np.geomspace(E_min, E_max, 500) *u.erg   # for multi star simulation    
print('Max E',np.argmax(energy))
print('e_conversion = ',e_conversion)
e_index = 17
print(energy[e_index].to(u.TeV))
print(energy[e_index].to(u.GeV))
print(energy[e_index].to(u.MeV))
print(pei_value)




radius_names = ['recent_venus', 'hz_mean', 'early_mars']
radius, energy = np.meshgrid(radius, energy)





print("BEFORE FIRST PLOT")
plt.rcParams.update({
    "axes.edgecolor": "black",
    "axes.linewidth": 2.0,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "lines.linewidth": 2.5,
    "legend.fontsize": 16
})
if plotting:
    NrE =N(radius,energy)
    MFP = ((3 * D_c(radius, energy)) / lightspeed.value)/AU
#:::::::::::::::::::::: Wind Velocity Profile Plotting ::::::::::::::::::::::::

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    ax.plot((radius/AU)[0,:], v(radius[0,:]).to('km/s'), label="r$V_w$ (Velocity Profile)")
    ax.axvline(x=R1_au_Value.value, color='purple', linestyle=':', label='$R_1$ (Wind Termination Shock)') 
    ax.axvline(x=RC_au_Value.value, color='green' , linestyle=':', label='$R_C$ (Contact Discontinuity)')
    ax.axvline(x=R_star_au.value, color='Black' , linestyle='-', label='$R_star$ (Stellar Radius)')
    ax.axhline(y=22, color='red', linestyle='--', label='\(R_{\*}\) (Stellar Velocity)')
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['recent_venus'], hz_boundaries_au['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
    ax.fill_betweenx(ax.get_ylim(), hz_boundaries_au['runaway_greenhouse'], hz_boundaries_au['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
    #::| The x axis is set to a logarithmic scale as the difference between R1 and
    #::| RC is more than a factor of one thousand
    #scaling
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    # Set custom y-ticks and labels
    # yticks = [1, 10, 100, 200, 300, 400]  # Example y-ticks
    # yticklabels = ['1', '10', '100', '200', '300', '400']  # Example y-tick labels
    # ax.set_yticks(yticks)
    # ax.set_yticklabels(yticklabels)
    #labels
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('Velocity[km/s]')

    ax.set_title('Wind Velocity profile for '+ Stellar_Values) 
    ax.legend()
    ax.grid(True)
    #plt.show()

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
    ax.grid(True)                                                                                                                                                                                                                                                           
    plt.show()

#::::::::::::::::::::::::: Cosmic Ray Profile Plotting ::::::::::::::::::::::::
    '''
    e_index = 17
    print('Energy Used for CR profile',energy[0])
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
     
    
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
    plotenergy = np.median(energy[e_index].to('MeV'))
    print('plotenergy',plotenergy)
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('N(r,E)')
    ax.set_title(f"Cosmic ray profile for {Stellar_Values} at {plotenergy:.0f}")
    ax.grid(True)
    plt.show()
    '''
#::::::::::::::::::::::: Mean Free Path plotting ::::::::::::::::::::::::::::::
    
    
    
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
    data1 = read_and_validate_csv("C:/Users/samir/Desktop/Stellar Bubbles/Parallel MFP.csv")
    data2 = read_and_validate_csv("C:/Users/samir/Desktop/Stellar Bubbles/Radial MFP.csv")
    '''
    e_index = 17
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
    
    #Plot the mean free path as a function of radius (in AU)
    ax.plot((radius / AU)[e_index, :], MFP[e_index, :], label=r'Parallel Mean Free Path $\lambda(r,E)$')
    plot_data(ax, data1, 'Parallel MFP Pei et al')
    plot_data(ax, data2, 'Radial MFP Pei et al')
    
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plotenergy = np.median(energy[e_index].to('MeV'))
    ax.set_xlabel('Radius [AU]')
    ax.set_ylabel('Mean Free Path [AU]')
    ax.set_title(f"Mean Free Path in {Stellar_Values} at {plotenergy:.0f}")    
    ax.grid(True)
    ax.legend()
    plt.show()
    '''
#::::::::::::::::::::::::: Cosmic Ray Spectra  Plotting :::::::::::::::::::::::
    print("BEFORE SECOND PLOT")

    #this finds the index of the value in the radius array that is clossest to hz_mean_cm.value
    r_index = np.argmin(abs(radius[0,:].value - (AU.value)))
    radius_in_au = radius[0, r_index].to(u.AU).value
    # Calculate residuals
    residuals = NrE[:, r_index] - cr.CR_modulation(energy[:, r_index].to(u.eV).value)
    print('radius_in_au =',radius_in_au)
     
    # write NrE to file for each value of R_hab(min med Max),   assign unique name to values in text based on 
    # SpV B_field valiue and V-star Value 
    
    # Create a plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=600,sharex = True , gridspec_kw={'height_ratios': [2, 1], 'hspace':0})
    
    # Plot the cosmic ray spectra 
    ax1.plot(energy[:, r_index].to(u.GeV), NrE[:, r_index], label="Cosmic Ray Spectra")
    ax1.plot(energy[:, r_index].to(u.GeV), cr.CR_modulation(energy[:, r_index].to(u.eV).value), label='Rogers-lee 2020')
    
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('GeV')
    ax1.set_ylabel('N')
    ax1.set_title(f'Cosmic Ray Spectra for {Stellar_Values} at $D_{{\mathrm{{ref}}}} =  {D_0:.1e}\\,\\mathrm{{cm^2\\,s^{{-1}}}}$ and R = 1 AU ')
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
    
    #:::::::::::::::::::::::::::::::::::::Teegarden plotting

    
    r_min = np.argmin(abs(radius[0,:].value - (recent_venus_cm.value)))
    r_max = np.argmin(abs(radius[0,:].value - (early_mars_cm.value)))
    Teegarden_b = 0.0252*AU
    Teegarden_c = 0.0443*AU
    Teegarden_b_r = np.argmin(abs(radius[0,:].value - Teegarden_b.value))
    Teegarden_c_r = np.argmin(abs(radius[0,:].value - Teegarden_c.value))
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('GeV ')
    ax.set_ylabel('Modulation')
    ax.set_title(f'Cosmic Ray Spectra for {Stellar_Values} at $D_{{\mathrm{{ref}}}} = {D_0:.1e}$ Over HZ')
    ax.fill_between(energy[:, r_min].to(u.GeV).value, NrE[:, r_min], NrE[:, r_max], color='gray', alpha=0.3, label='Between HZ Min and Max')
    ax.plot(energy[:,Teegarden_b_r].to(u.GeV), NrE[:, Teegarden_b_r],color = 'green', label = 'Teegarden_b')
    ax.plot(energy[:,Teegarden_c_r].to(u.GeV), NrE[:, Teegarden_c_r],color = 'blue' ,label = 'Teegarden_c')
    ax.grid(True)
    ax.legend()
    plt.show()
    
    '''
    '''
    #::::::::::::::::::::::::::::::::::: Diagram Plotting:::::::::::::::::::::
    
    # Convert to AU for plotting if necessary
    radius_au = radius.to(u.au).value
    
    # Initialize the figure and axis
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
    ax.set_aspect('equal')
    
   # Redefine the plot_circle function to accept a label
    def plot_circle(r, color, alpha=1.0, linestyle='-', label=None):
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(theta, np.full_like(theta, r), color=color, alpha=alpha, linestyle=linestyle, label=label)

    
    # Plot 'Sol' at the center with a star symbol
    sol = patches.RegularPolygon((0, 0), numVertices=8, radius=R_star_au*50, transform=ax.transData._b, color='yellow')
    ax.add_patch(sol)
    
    # Function to plot a circle in polar coordinates and fill between two radii
    def fill_between_circles(r_inner, r_outer, color, alpha=1.0, label=None):
        theta = np.linspace(0, 2*np.pi, 100)
        r_inner_circle = np.full_like(theta, r_inner)
        r_outer_circle = np.full_like(theta, r_outer)
        ax.fill_between(theta, r_inner_circle, r_outer_circle, color=color, alpha=alpha, label=label)
    
    # Plot the Habitable Zones as filled areas
    fill_between_circles(hz_boundaries_au['recent_venus'], hz_boundaries_au['early_mars'], 'green', alpha=0.2, label='Optimistic HZ')
    fill_between_circles(hz_boundaries_au['runaway_greenhouse'], hz_boundaries_au['max_greenhouse'], 'green', alpha=0.3, label='Pessimistic HZ')
    plot_circle((1 * u.AU).value, 'blue', linestyle='-',label='Earth Orbit')
    # Plot the Wind Termination Shock and Contact Discontinuity
    plot_circle(R1_au_Value.value, 'purple', linestyle=':',label='Wind Termination Shock')
    plot_circle(RC_au_Value.value, 'blue', linestyle=':',label='Contact Discontuinity')
    
    # Use the radius array to set the radial limits
    ax.set_rscale('log')
    ax.set_rlim(bottom=np.min(radius_au), top=np.max(radius_au))
    
    # Hide the polar grid and labels
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.spines['polar'].set_visible(False)
    ax.legend(loc='upper right', fontsize='small', fancybox=True)
    # Add a legend manually
    # Your legend code here...
    ax.set_title(Spectral_Class +' System Log scale Diagram', va='top',pad=20)
    #ax.set_xlabel('Distance (AU)')
    
    
    
    
    max_radius = np.max(np.abs(radius_au))
    # Add an inset axis for the bottom X-axis with symmetrical log scale tick marks
    # [left, bottom, width, height] in figure coordinate
    inset_axis = fig.add_axes([0.1, 0.05, 0.8, 0.02]) 
    inset_axis.set_xlim(-max_radius, max_radius)
    # Set up the bottom axis with symmetrical log scale using symlog function
    inset_axis.set_xscale('symlog', linthresh=0.03, linscale=0.03)
    
    # Set limits for the x-axis based on the max radius in your data
    max_radius = np.max(radius_au[radius_au > 0])
    inset_axis.set_xlim(-max_radius, max_radius)
    tick_vals = [0.1, 0.5, 1, 5, 10, 50, 100]*u.AU
    tick_vals = tick_vals.value
    # Define the tick positions on the log scale, making sure they are symmetrical around 0
    #tick_vals = np.geomspace(start=0.2, stop=max_radius, num=4)
    
    ticks = np.concatenate((-tick_vals[::-1], [0], tick_vals))  # symmetrical ticks around 0
    
    # Custom formatter for the ticks on the symlog scale
    # def custom_formatter(value, tick_number):
    #     # Using "{:.2g}" for compact numerical representation
    #     return "{:.4f} AU".format(value) if value != 0 else "0"
    # Define a new formatter function to display AU labels
    def custom_formatter(value, tick_number):
        return "{:.1f} AU".format(value)    
    
    
    
    # Set the ticks and the custom formatter ##Works
    inset_axis.set_xticks(ticks)
    inset_axis.tick_params(axis='x', labelsize=8,labelrotation=45)
    inset_axis.xaxis.set_major_formatter(plt.FuncFormatter(custom_formatter))
    
    # Hide the y-axis and the spines for the top, left, and right
    inset_axis.yaxis.set_visible(False)
    inset_axis.spines['top'].set_visible(False)
    inset_axis.spines['left'].set_visible(False)
    inset_axis.spines['right'].set_visible(False)
    plt.subplots_adjust(bottom=0.1, top=0.9,right=0.9,left=0.1)
    # Optionally, add grid lines to the inset_axis for better readability
    inset_axis.grid(True, which='both', axis='x', linestyle='--', linewidth=0.7)
    
    # Display the plot
    plt.show()
    
else:
    print('No Plots')
    NrE =N(radius,energy)
    print(np.shape(NrE))
    print(np.shape(radius))
    print(np.shape(energy))
    filename = "Data/" + Stellar_Values+"_"+'B'+str(B0)+"_"+'V'+str(v_star_km.value)
    print("ABOUT TO SAVE TO:")
    print(filename)
    print(energy.dtype)
    #radius.flatten, energy.flatten, NrE.flatten
    print(radius, radius.flatten())
    print(energy, energy.flatten())
    print(os.path.abspath(filename))
    np.savetxt(filename, np.array([radius.flatten(), energy.flatten(), NrE.flatten()]).transpose(), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ')
    print(np.shape(radius.flatten()))
    
print("END SCRIPT")

# # out = out.reshape((int(np.sqrt(out.shape[0])),)*2 + (out.shape[-1],))
# # out = out.swapaxes(0,-1).swapaxes





# e_index = 70
# r_index = np.argmin(abs(radius[0,:].value - (hz_mean_cm.value))) # possible issue as R index should be at 1 AU? 
# radius_in_au = radius[0, r_index].to(u.AU).value

# R_values= [recent_venus_cm, hz_mean_cm, early_mars_cm]
# for R in R_values:
#     r_index = np.argmin(abs(radius[0,:].value - R.value))
    
#     filename = f"CR_Spectra_{Stellar_Values}_B{B0}_V{v_star_km.value}_R{R.value}.txt"
#     with open(filename, 'w') as file:
#         for e, n in zip(energy[:, r_index], NrE[:, r_index]):
#             file.write(f"{e.value} {n}\n")

#np.reshape 
# print(f"Energy at e_index ({e_index}): {energy[e_index,0].value} erg")
# print(f"Radius at r_index ({r_index}): {radius[0,r_index].value} cm")
# print(f"Radius at r_index ({r_index}): {radius_in_au} AU")
