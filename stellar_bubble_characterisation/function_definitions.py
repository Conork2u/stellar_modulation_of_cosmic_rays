# -*- coding: utf-8 -*-
"""
Created on Thu May 18 16:23:33 2023

@author: conor
"""
import numpy as np
from stellar_parameters import *
from astropy import units as u
from astropy import constants as const
import sys
import sympy as sp
#:::::::::::::::::::::::::::::::::::::::::::::::::::: Function Definitions ::::
print('----------------------------------')
print('Start of Function Definitions')

#::| The Total Mass lossed by the star over its lifetime is calculated
def constant_mass_loss_rate(t):
    mass_loss_rate_value = mass_loss_rate.value
    return mass_loss_rate_value
start_time = 0.0  
end_time = stellar_age.value
result, _ = quad(constant_mass_loss_rate, start_time, end_time)
MII = result * u.g
#print(f'Integrated mass loss = {MII:.2E}')
#::| The Pressure is calculated, also = to the pressure in region I and II
def pISM(n, k_B, T):
    return(T * k_B * n)
#print('pISM= ', pISM(n, k_B, T))

#::| The Radius of the termination shock (R1) and the contact discontiniuty(RC)
#::| are calculated and converted in to values suitable for scaling the plot
def R1(lw, vw, pISM):
    return np.sqrt(lw / (vw * (2 * np.pi) * pISM(n, k_B, T))).value
R1_Value = R1(lw, vw, pISM) * u.cm
#print(Stellar_Values + ' R1 =', f'{R1_Value:.2E}')
R1_pc_Value = R1_Value.to('pc')
R1_au_Value = R1_Value.to('au')
#print(Stellar_Values + ' R1 =', f'{R1_au_Value}')
R1_Valueplus = R1_Value + 1*u.cm
R1_Valueminus = R1_Value - 1*u.cm

def RC(t, vw, R1_Value, m_p):
    return(np.cbrt(3/4. * (R1_Value**2).value * vw.value * t.value  + (R1_Value**3).value))
RC_Value = RC(stellar_age, vw, R1_Value, m_p) *u.cm
RC_pc_Value = RC_Value.to('pc')
RC_au_Value = RC_Value.to('au')
RC_Valueplus =  RC_Value + 0.1*u.cm
RC_Valueminus = RC_Value - 10000*u.cm
#print(Stellar_Values + ' RC =', f'{RC_Value:.2E}')

def rho_R1(R1_Value):
    return (mass_loss_rate.value / (4 * np.pi * R1_Value.value **2 *vw.value))/ m_p.value
#print(rho_R1(R1_Value))
R_star_pc = R_star / parsec_in_cm

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::: Velocity Profile ::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
#::| The velocity profile is defined, always taking cm as the input unit
def v(radius):
    Zone_A = np.where(radius < R1_Value, vw, 0.) 
    Zone_B = np.where((radius >= R1_Value)*(radius <= RC_Value), (vw/4) * (R1_Value**2/radius**2).value, Zone_A)
    Zone_C = np.where(radius > RC_Value, 0. ,Zone_B )
    return Zone_C
#print('V(RC+) =', v(RC_Value + 1000*u.cm))
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::: Magnetic Field Profile ::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def B(radius):
    # 
    Zone_A = np.where((0 < radius) * (radius < R_star) , 0.,  0.)
    Zone_B = np.where((R_star <= radius) * (radius <= R1_Value), B0 * R_star / radius, Zone_A)
    Zone_C = np.where((R1_Value < radius)* (radius <= RC_Value), B0*4 * R_star / R1_Value, Zone_B)
    Zone_D = np.where(radius > RC_Value, 0., Zone_C)
    return Zone_D    


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::: Cosmic Ray Profile ::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::|Variables for the comsic ray fucntion   ::
c, d = -0.5, 0.5
N_0 = 1 #   initial cosmic ray density 
B_ISM = 5**-6 #Gauss
D_0 = 1.8e26 # cm^2/s ## Ideally ~~ 10 ^ 26
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





def argument_checks(): ### exp(argumnet) checks
    ###          Zone B Argument checks
    print('__________Zone B Argument Checks_____')
    a , b = -2 , 0
    C = 1
    delta = (a - (b * c) +1)
    print("delta = ", delta)
    RC_test_Value = RC_Value/2
    print(RC_test_Value)
    print('RC/RC =',RC_Value.value/RC_Value.value)
    print('D = Delta = ',delta)                
    print('1/D =',1./(delta) )                     
    print('D*R*D =', (-1 * (1./(delta) * (RC_Value.value/RC_Value.value) ** delta)))
    print('D*C*D =',(1./(delta) * (C) ** delta))                                  
    print('D*R*D + D*C*D =', (-1 * (1./(delta) * (RC_Value.value/RC_Value.value) ** delta) + (1./(delta) * (C) ** delta)))
    print('D*R/2*D + D*C* , , D =', (-1 * (1./(delta) * (RC_test_Value.value/RC_Value.value) ** delta) + (1./(delta) * (C) ** delta)))
    print('exp(arg)', np.exp((-1 * (1./(delta) * (RC_Value.value/RC_Value.value) ** delta) + (1./(delta) * (C) ** delta))))
    print(N_B(R1_Valueplus))
    ###             Zone A Argument checks
    print('__________Zone A Argument Checks_____')
    a , b = 0 , -1
    C = 1
    delta = (a - (b * c) +1)
    print("delta = ", delta)
    test_Value = R1_Value /2
    print('R1/R1 =',R1_Value.value/R1_Value.value)
    print('D = Delta = ',delta)                
    print('1/D =',1./(delta) )                     
    print('D*R1*D =', (-1 * (1./(delta) * (R1_Value.value/R1_Value.value) ** delta)))
    print('D*C*D =',(1./(delta) * (C) ** delta))                                  
    print('D*R1*D + D*C*D =', (-1 * (1./(delta) * (R1_Value.value/R1_Value.value) ** delta) + (1./(delta) * (C) ** delta)))
    print('D*R test*D + D*C*D =', (-1 * (1./(delta) * (test_Value.value/R1_Value.value) ** delta) + (1./(delta) * (C) ** delta)))
    print('exp(arg) at R1', np.exp(alpha * (-1 * (1./(delta) * (R1_Value.value/R1_Value.value) ** delta) + (1./(delta) * (C) ** delta))))
    print('exp(arg) at test R ',np.exp(alpha * (-1 * (1./(delta) * (test_Value.value/R1_Value.value) ** delta) + (1./(delta) * (C) ** delta))) )

#argument_checks()
   

def prefactor_checks():### Prefactor and equation results checks
    print(f"  :{v(R1_Value+1*u.cm)}")
    print('________________________________________')
    print('Prefactor checks:::::::::::::::::::::')
    print(f"v(R1 Plus)  :{v(R1_Value+1*u.cm)}")
    print(f"v(R1 Minus) :{v(R1_Value-1*u.cm)}")
    print(f"B(R1 Minus) :{B(R1_Valueminus)}")
    print(f"B(R1 Plus)  :{B(R1_Valueplus)}")
    print(f"B(RC Minus) :{B(RC_Valueminus)}")
    print(f"B(R Star)   :  {B(R_star)}")
    #print('Beta        :', beta)
    #print('Alpha       :', alpha)
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

        

#prefactor_checks() 





def N(radius,E): 
    Zone_ISM  = np.where(radius > RC_Value, 1.0, 1.0)
    Zone_B    = np.where((R1_Value < radius) * (radius < RC_Value), N_B(radius,E), Zone_ISM)
    Zone_A    = np.where((R_star < radius) * (radius < R1_Value), N_A(radius,E), Zone_B)
    Zone_Star = np.where((0 < radius) * (radius <= R_star), 0.0, Zone_A)
    return Zone_Star

print('End of Function Definitions')
print('----------------------------------')















