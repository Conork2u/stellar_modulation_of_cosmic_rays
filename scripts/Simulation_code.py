# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 17:59:31 2024

@author: conor
"""

import itertools
import subprocess
import sys
import os
import astropy
from astropy import units as u
from astropy import constants as const

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)


if not os.path.exists('Data'):
    os.makedirs('Data')
    print('Data file created')
else:
    print('path exists')


# from Bouvier et al
# Define your parameters
spectral_classes = [
    ('F5V', 6550, 5E-14, 1.473, 0.56, 500),
    ('G5V', 5660, 2.4E-14, 0.977, -0.05, 400),
    ('K5V', 4440, 2.99E-15, 0.701, -0.76, 300),
    ('M5V', 3060, 1.92E-16, 0.196, -2.52, 200)
]


magnetic_fields = {
    'F5V': (0.8, 2.2, 3.3),
    'G5V': (0.2, 1, 47.9),
    'K5V': (1.50, 36.5, 78),
    'M5V': (16, 51, 406)
}

stellar_velocities = (10, 20, 30, 60, 120)

# Function to run the script with a set of parameters
def run_script(parameters):
    command = [sys.executable, "Bubble_model.py", *map(str, parameters)]
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode == 0:
        return result.stdout
    else:
        print(f"Error: {result.stderr}")
        return None
 
specific_systems = [
    ('Trappist1', 2570, 6.9e-17, 0.114, -3.28, 200, 600, 80.81),
    ('Proxima',   2930, 3.5e-17, 0.156, -2.79, 400, 600, 25.0),
    ('Teegarden', 2680, 6.8e-16, 0.120, -3.19, 200, 1200, 92.0),
]

for params in specific_systems:
    print(params)
    output = run_script(params)
    if output is not None:
        filename = 'Data/' + '_'.join(map(str, params))
        with open(filename, 'w') as file:
            file.write(output)

# Sol specific parameters
sol_parameters = ('Sol', 5780, 6.8E-14, 1, 0, 400, 1, 22)
print(sol_parameters)
output = run_script(sol_parameters)
if output is not None:
    filename = 'Data/' + '_'.join(map(str, sol_parameters)) + ".0"
    with open(filename, 'w') as file:
        file.write(output)


'''
# Iterate over each spectral class
for sp_class in spectral_classes:
    # Iterate over each magnetic field strength for the current spectral class
    for b_field in magnetic_fields[sp_class[0]]:
        # Iterate over each stellar velocity
        for v_star in stellar_velocities:
            # Combine the current set of parameters
            current_parameters = sp_class + (b_field, v_star)
            print(current_parameters)
            
            # Run the script with the current set of parameters
            output = run_script(current_parameters)
            # if output is not None:
            #     filename = 'Data/' + '_'.join(map(str, current_parameters)) + ".0"
            #     with open(filename, 'w') as file:
            #         file.write(output)}  simread.py{# -*- coding: utf-8 -*-
'''
