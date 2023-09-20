# stellar_modulation_of_cosmic_rays
Description
This project calculates stellar parameters and habitable zones for different main-sequence stars. It uses various scientific constants and mathematical models, including Kopparapu et al.'s work on habitable zones.

Dependencies
Python 3.x
Matplotlib
NumPy
SciPy
Astropy
Installation
To install all required packages, you can use pip:

bash
Copy code
pip install matplotlib numpy scipy astropy
Usage

To run the program:
Set your desired star by editing the Stellar_Values variable in stellar_parameters.py.
current options srae sol , proxima_centauri and TRAPPIST-1.
New stellar values need to be pulled from the literature for the model to work. Eventually, we will have every stellar type from F0 to M8. 
Once Stellar parameters are set, run the plotting file to produce a graph showing the modulation of cosmic rays. 

Files
stellar_parameters.py: Script to calculate stellar parameters such as mass loss rate, terminal velocity, and mechanical luminosity for a given star.
habitable_zone_calculation.py: Script that calculates the habitable zone based on the Kopparapu et al. model.
functions_definitions.py: Script that calculates the values for the stellar wind, magnetic field and cosmic ray density per radius and energy

Credits
The habitable zone calculations are based on:
"Habitable Zones Around Main-Sequence Stars: New Estimates" by Kopparapu et al.(2013), Astrophysical Journal, 765, 131
"Habitable Zones Around Main-Sequence Stars: Dependence on Planetary Mass" by Kopparapu et al.(2014), Astrophysical Journal Letters, 787, L29
 Habitable zone calculator Translated to Python by John Armstrong.
