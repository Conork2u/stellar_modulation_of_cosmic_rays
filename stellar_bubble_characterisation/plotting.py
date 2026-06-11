# -*- coding: utf-8 -*-



"""
Created on Thu May 18 16:23:35 2023

@author: conor
"""

import matplotlib.pyplot as plt
from function_definitions import *
from hz_calculations import hz_boundaries_pc
from hz_calculations import hz_mean_cm
from stellar_parameters import parsec_in_cm
import numpy as np
import sys
from astropy import units as u
from astropy import constants as const
import CRspectra as cr

E_min = 0.001602 / 100 * 1.6
E_max = 0.00162  * 1000

energy = np.geomspace(E_min, E_max, 100) *u.erg

radius = np.geomspace(R_star.value, 1.3*RC_Value.value, 100000) * u.cm


print(energy.value)
#print(energy)

radenergy = np.meshgrid(radius, energy)
#print(np.shape(radenergy))
radius = radenergy[0]
energy = radenergy[1]
#print(energy)
#print(radius)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::: Wind Velocity Profile Plotting ::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
print('----------------------------------')
print('Start of Wind Velocity Plotting')    

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
ax.plot((radius/parsec_in_cm)[0,:], v(radius[0,:]).to('km/s'), label="Velocity Profile")
ax.axvline(x=R1_pc_Value.value, color='purple', linestyle=':', label='Wind Termination Shock') 
ax.axvline(x=RC_pc_Value.value, color='green' , linestyle=':', label='Contact Discontinuity')
ax.axvline(x=R_star_pc.value, color='Black' , linestyle='-', label='Stellar Radius')

ax.fill_betweenx(ax.get_ylim(), hz_boundaries_pc['recent_venus'], hz_boundaries_pc['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
ax.fill_betweenx(ax.get_ylim(), hz_boundaries_pc['runaway_greenhouse'], hz_boundaries_pc['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
#::| The x axis is set to a logarithmic scale as the difference between R1 and
#::| RC is more than a factor of one thousand
#scaling
ax.set_xscale('log')
#ax.set_yscale('log')

#labels
ax.set_xlabel('Radius [pc]')
ax.set_ylabel('Velocity[km/s]')
ax.set_title('Wind Velocity profile  for '+ Stellar_Values )
ax.legend()

plt.show()

print('End of Wind Velocity Plotting')
print('----------------------------------') 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::: Magnetic Field Profile Plotting :::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

print('----------------------------------')
print('Start of Magnetic Field Plotting')    

fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
ax.plot((radius/parsec_in_cm)[0,:], B(radius[0,:]), label=" Magnetic Field Profile")
ax.axvline(x=R1_pc_Value.value, color='purple', linestyle=':', label='Wind Termination Shock') 
ax.axvline(x=RC_pc_Value.value, color='green' , linestyle=':', label='Contact Discontinuity')
ax.axvline(x=R_star_pc.value, color='Black' , linestyle='-', label='Stellar Radius')

ax.fill_betweenx(ax.get_ylim(), hz_boundaries_pc['recent_venus'], hz_boundaries_pc['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
ax.fill_betweenx(ax.get_ylim(), hz_boundaries_pc['runaway_greenhouse'], hz_boundaries_pc['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
#scaling
ax.set_xscale('log')
ax.set_yscale('log')

#labels
ax.set_xlabel('Radius [pc]')
ax.set_ylabel('Magentic Field (Gauss)')
ax.set_title(' Magnetic Field Profile for '+ Stellar_Values )
ax.legend()
                                                                                                                                                                                                                                                                                    
plt.show()

print('End of Magnetic Field Plotting')
print('----------------------------------') 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::: Cosmic Ray Profile Plotting :::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

print('----------------------------------')
print('Start of Cosmic Ray Plotting')    

e_index = 70
print(energy[e_index,0])   
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
NrE =N(radius,energy)
ax.plot((radius/parsec_in_cm)[e_index,:],NrE[e_index,:], label=" Cosmic Ray Profile")
ax.axvline(x=R1_pc_Value.value, color='purple', linestyle=':', label='Wind Termination Shock') 
ax.axvline(x=RC_pc_Value.value, color='green' , linestyle=':', label='Contact Discontinuity')
ax.axvline(x=R_star_pc.value, color='Black' , linestyle='-', label='Stellar Radius')


ax.fill_betweenx(ax.get_ylim(), hz_boundaries_pc['recent_venus'], hz_boundaries_pc['early_mars'], color='green', alpha=0.2,label = 'Optimistic Habitable Zone')
ax.fill_betweenx(ax.get_ylim(), hz_boundaries_pc['runaway_greenhouse'], hz_boundaries_pc['max_greenhouse'], color='green', alpha=0.5,label = 'Pesimistic Habitable Zone')
#scaling
ax.set_xscale('log')
ax.set_yscale('log')


#labels
ax.set_xlabel('Radius [pc]')
ax.set_ylabel('N(r,E)')
ax.set_title(' Cosmic Ray Profile for '+ Stellar_Values )
ax.legend()
plt.show()
print(hz_boundaries_pc['recent_venus'])

print('End of Cosmic Ray Plotting')
print('----------------------------------') 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::: Cosmic Ray Spectra  Plotting :::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''
print('----------------------------------')
print('Start of Cosmic Ray Spectra Plotting')    

 
print(f"Habitable Zone mean:{hz_mean_cm}")
r_index = np.argmin(abs(radius[0,:].value - (hz_mean_cm.value)))  
print(r_index, radius[0, r_index])  

fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
ax.plot(energy[:,r_index].to(u.GeV), NrE[:,r_index] , label=" Cosmic Ray Spectra")
ax.plot(energy[:,r_index].to(u.GeV), cr.CR_modulation(energy[:,r_index].to(u.eV).value), label='Rogers-lee 2020')

#scaling
ax.set_xscale('log')
ax.set_yscale('log')
#10 Gev = 0.01602 ergs in cgs
#1 Gev = 0.00162 ergs

#labels
ax.set_xlabel('GeV')
ax.set_ylabel('N')
ax.set_title(' Cosmic Ray Spectra for '+ Stellar_Values)
ax.legend()
plt.show()
print('End of Cosmic Ray Spectra Plotting')
print('----------------------------------') 
'''

def equation_checks():
    print('________________________________________')
    print('Equation result checks:::::::::::::::')
    print('----------------Zone B-----------------')
    print(f"N_B  at RC    : {N_B(RC_Value)}")
    print(f"N_B  at -RC   : {N_B(RC_Valueminus)}")
    print(f"N_B  at 1/2RC : {N_B(RC_Value/2)}")
    print(f"N_B  at +R1!  : {N_B(R1_Valueplus)}")
    print(f"N_B  at R1    : {N_B(R1_Value)}")
    print('----------------Zone A-----------------')
    print(f"N_A at  R1    : {N_A(R1_Value)}")
    print(f"N_A at  -R1   : {N_A(R1_Valueminus)}")
    print(f"N_A at  1/2R1 : {N_A(R1_Value/2)}")
    print(f"N_A at R_star : {N_A(R_star)}")

    
