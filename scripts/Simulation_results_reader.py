# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 16:39:04 2024

@author: conor
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy import units as u




def normalize(value, min_value, max_value):
    """ Normalize a value to a range between 0 and 1. """
    return (value - min_value) / (max_value - min_value)
def get_middle_value(items):
    return sorted(items)[len(items) // 2]

def extract_parameters_from_filename(filename):
    # Extracts spectral class, magnetic field, and v_star from a filename.
    parts = filename.split('_')
    spectral_class = parts[0]
    B_field = float(parts[1][1:])
    V_star = float(parts[2][1:-2])
    return spectral_class, B_field, V_star

def get_unique_parameters(data_dir):
    spectral_classes = set()
    magnetic_fields = {}
    stellar_velocities = set()

    for file in os.listdir(data_dir):
        if file.endswith('.0'):
            sp_class, B_field, V_star = extract_parameters_from_filename(file)
            spectral_classes.add(sp_class)
            magnetic_fields.setdefault(sp_class, set()).add(B_field)
            stellar_velocities.add(V_star)

    return spectral_classes, magnetic_fields, stellar_velocities

def plot_spectral_class_data(data_dir, spectral_classes, magnetic_fields, stellar_velocities):
    color_map = {
        'F5V': 'lightblue',
        'G5V': 'gold',
        'K5V': 'orange',
        'M5V': 'red',
        'Sol': 'black'
    }

    # Function to plot data for a specific spectral class
    def plot_data_for_class(ax, sp_class, color):
        for file in os.listdir(data_dir):
            if file.startswith(sp_class) and file.endswith('.0'):
                r, E, data = np.loadtxt(os.path.join(data_dir, file), unpack=True)
                NrE = np.reshape(data, (500, 3))
                r = np.reshape(r, (500, 3))
                E = np.reshape(E, (500, 3))
                Eplot = E[:, 0]*u.erg.to('GeV')
                Nmin = NrE[:, 0]
                Nmid = NrE[:, 1]
                Nmax = NrE[:, 2]
                ax.fill_between(Eplot, Nmin, Nmax, color=color, alpha=0.4)

    for sp_class in spectral_classes:
        fig, ax = plt.subplots()
        
        # Plot data for the specific spectral class
        color = color_map.get(sp_class, 'Black')
        plot_data_for_class(ax, sp_class, color)

        # Plot data for 'Sol', in every plot
        sol_color = color_map['Sol']
        if sp_class != 'Sol':  # Avoid duplicating 'Sol' data
            plot_data_for_class(ax, 'Sol', sol_color)

          # Create legend handles
        class_handle = plt.Line2D([0], [0], color=color, lw=4, label=sp_class)
        handles = [class_handle]
        if sp_class != 'Sol':  
            sol_handle = plt.Line2D([0], [0], color=sol_color, lw=4, label='Sol')
            handles.append(sol_handle)
        ax.legend(loc='center right', fontsize='small',handles=handles)
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy (GeV)')
        ax.set_ylabel('NrE (Modulation)')
        ax.set_title(f'{sp_class}: Energy vs NrE')
        ax.grid(True)
        plt.show()
def plot_varying_parameters(data_dir, spectral_classes, magnetic_fields, stellar_velocities):
    color_map = {
        'F5V': 'lightblue',
        'G5V': 'green',
        'K5V': 'orange',
        'M5V': 'red',
        'Sol': 'black'
    }

    def plot_data(ax, files, color, label, is_v_star=True):
        sorted_files = sorted(files, key=lambda f: float(f.split('_')[2][1:-2]))  # Sort files by V_star or B_field
        for i, file in enumerate(sorted_files):
            _, B_field, V_star = extract_parameters_from_filename(file)
            alpha = 1 - (i / len(sorted_files))  # Alpha increases with higher V_star or B_field
            label = f"V={V_star:.1f}" if is_v_star else f"B={B_field:.1f}"
            
            r, E, data = np.loadtxt(os.path.join(data_dir, file), unpack=True)
            NrE = np.reshape(data, (500, 3))
            r = np.reshape(r, (500, 3))
            E = np.reshape(E, (500, 3))
            Eplot = E[:, 0]*u.erg.to('GeV')
            Nmin = NrE[:, 0]
            Nmid = NrE[:, 1]
            Nmax = NrE[:, 2]
            ax.fill_between(Eplot, Nmin, Nmax, color=color, alpha=alpha, label=label)

    for sp_class in spectral_classes:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        color = color_map.get(sp_class, 'black')

        # Plot 1: Varying V star
        middle_B_field = get_middle_value(magnetic_fields[sp_class])
        B_field_files = [file for file in os.listdir(data_dir) if file.startswith(f"{sp_class}_B{middle_B_field:.1f}") and file.endswith('.0')]
        plot_data(ax1, B_field_files, color, f'{sp_class}', is_v_star=True)

        # Plot 2: Varying B field
        middle_V_star = get_middle_value(stellar_velocities)
        V_star_files = [file for file in os.listdir(data_dir) if file.startswith(sp_class) and f"_V{middle_V_star:.1f}" in file]
        plot_data(ax2, V_star_files, color, f'{sp_class}', is_v_star=False)

        # Sol plotting
        sol_files = [file for file in os.listdir(data_dir) if file.startswith('Sol') and file.endswith('.0')]
        sol_color = color_map['Sol']
        plot_data(ax1, sol_files, sol_color, 'Sol', is_v_star=True)
        plot_data(ax2, sol_files, sol_color, 'Sol', is_v_star=False)

        # Axes properties
        for ax in [ax1, ax2]:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel('Energy (GeV)')
            ax.set_ylabel('NrE (Modulation)')
            ax.legend(loc='center right', fontsize='small')

        ax1.set_title(f'{sp_class} - Varying V star (B={middle_B_field} Gauss)')
        ax2.set_title(f'{sp_class} - Varying B field (V={middle_V_star} km/s)')
        ax.grid(True)
        plt.tight_layout()
        plt.show()    




def plot_extreme_modulations(data_dir, spectral_classes):
    color_map = {
        'F5V': 'blue',
        'G5V': 'green',
        'K5V': 'orange',
        'M5V': 'red',
        'Sol': 'black'
    }
    alpha_map = {
        'F5V': 0.25,
        'G5V': 0.25,
        'K5V': 0.25,
        'M5V': 0.25,
        'Sol': 0.4  # Less transparency for 'Sol' to highlight it
    }
    hatch_map = {
        'F5V': None,
        'G5V': '+',
        'K5V': 'O',
        'M5V': '//',
        'Sol': None
    }

    def plot_class(sp_class, ax):
        max_modulation = np.zeros((500, 3))  # Initialize with zeros
        min_modulation = np.ones((500, 3)) * np.inf  # Initialize with infinity
        
        # Find max and min modulations for the spectral class
        for file in os.listdir(data_dir):
            if file.startswith(sp_class) and file.endswith('.0'):
                r, E, data = np.loadtxt(os.path.join(data_dir, file), unpack=True)
                E = np.reshape(E, (500, 3))
                r = np.reshape(r, (500, 3))
                NrE = np.reshape(data, (500, 3))
                max_modulation = np.maximum(max_modulation, NrE)  # Update max modulation
                min_modulation = np.minimum(min_modulation, NrE)  # Update min modulation

         # Plot the fill between the max and min modulations
        Eplot = E[:, 0] * u.erg.to('GeV')
        Nmin = min_modulation[:, 0]
        Nmax = max_modulation[:, 2]
        alpha_value = alpha_map.get(sp_class, 0.5) 
        hatch_value = hatch_map.get(sp_class)
        ax.fill_between(Eplot, Nmin, Nmax, color=color_map.get(sp_class, 'black'), alpha=alpha_value,hatch=hatch_value, label=sp_class)

    for sp_class in spectral_classes:
        fig, ax = plt.subplots(dpi=300)
        
        # Plot 'Sol' in each plot
        

        # Plot the current spectral class
        plot_class(sp_class, ax)
        plot_class('Sol', ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy (GeV)')
        ax.set_ylabel('NrE (Modulation)')
        ax.set_title(f'Modulation range for {sp_class} including Sol')
        ax.legend(loc='upper right')
        ax.grid(True)
        plt.show()
        
     # Create one figure for all spectral classes
    spectral_classes = ['F5V', 'G5V', 'K5V', 'M5V', 'Sol']
    fig, ax = plt.subplots(dpi=300)

    # Plot each spectral class on the same plot
    for sp_class in spectral_classes:
        plot_class(sp_class, ax)
        ax.set_ylim(1e-3,1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('NrE (Modulation)')
    ax.set_title('Modulation range for each Spectral Class')
    ax.legend(loc='upper right')
    ax.grid(True)
    plt.show()
def plot_solar_modulation_only(data_dir):
    color = 'black'

    fig, ax = plt.subplots(dpi=300)

    sol_files = [file for file in os.listdir(data_dir) if file.startswith('Sol') and file.endswith('.0')]

    for file in sol_files:
        r, E, data = np.loadtxt(os.path.join(data_dir, file), unpack=True)

        NrE = np.reshape(data, (500, 3))
        E = np.reshape(E, (500, 3))

        Eplot = E[:, 0] * u.erg.to('GeV')
        Nmin = NrE[:, 0]
        Nmid = NrE[:, 1]
        Nmax = NrE[:, 2]

        ax.fill_between(Eplot, Nmin, Nmax, color=color, alpha=0.35)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('NrE (Modulation)')
    ax.set_title('Comparison of modulation across known systems')
    ax.legend(loc='upper right')
    ax.grid(True)
    plt.tight_layout()
    plt.show()
# Data directory
data_directory = 'Data/'

# Extract parameters from filenames
spectral_classes, magnetic_fields, stellar_velocities = get_unique_parameters(data_directory)

# Plotting tasks
plot_solar_modulation_only(data_directory)
plot_extreme_modulations(data_directory, spectral_classes)
plot_spectral_class_data(data_directory, spectral_classes, magnetic_fields, stellar_velocities)
plot_varying_parameters(data_directory, spectral_classes, magnetic_fields, stellar_velocities)
