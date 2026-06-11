# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 22:15:33 2024

@author: conor
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
print("Plotting")
plt.rcParams.update({
    "axes.edgecolor": "black",
    "axes.linewidth": 2.0,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "lines.linewidth": 2.5,
    "legend.fontsize": 16
})

plot = 1

if plot == 1:

    df = pd.read_csv("crs-and-exoplanets-main/modulation_results.csv")

    color_mapping = {
        'F': 'blue',
        'G': 'green',
        'K': 'orange',
        'M': 'red',
        'S': 'black'
    }

    df['Base_Spectral_Class'] = df['Spectral_Class'].str[0]
    df['Magnetic_Field'] = df['Spectral_Class'].str.extract(r'_B([0-9.]+)_').astype(float)
    df['Velocity'] = df['Spectral_Class'].str.extract(r'_V([0-9.]+)').astype(float)

    def get_mid_B(df_class):
        B_vals = np.sort(df_class['Magnetic_Field'].unique())
        return B_vals[len(B_vals) // 2]

    s_values = df['S_Value'].unique()

    for s in s_values:
        df_s = df[df['S_Value'] == s]

        fig, ax = plt.subplots(figsize=(10, 6))

        classes = ['F', 'G', 'K', 'M']
        labels = ['F', 'G', 'K', 'M']

        for i, sc in enumerate(classes):
            df_class = df_s[df_s['Base_Spectral_Class'] == sc]

            vals_30 = df_class['Ratio_Above_30GeV']
            vals_all = df_class['Ratio_All_Energy']

            if len(vals_30) > 0:
                vmin_30 = vals_30.min()
                vmax_30 = vals_30.max()
                vcent_30 = 0.5 * (vmin_30 + vmax_30)

                ax.errorbar(
                    i,
                    vcent_30,
                    yerr=[[vcent_30 - vmin_30], [vmax_30 - vcent_30]],
                    fmt='none',
                    color=color_mapping[sc],
                    alpha=1.0,
                    capsize=20,
                    elinewidth=4,
                    capthick=4,
                    zorder=1
                )

            if len(vals_all) > 0:
                vmin_all = vals_all.min()
                vmax_all = vals_all.max()
                vcent_all = 0.5 * (vmin_all + vmax_all)

                ax.errorbar(
                    i,
                    vcent_all,
                    yerr=[[vcent_all - vmin_all], [vmax_all - vcent_all]],
                    fmt='none',
                    color=color_mapping[sc],
                    alpha=0.3,
                    capsize=20,
                    elinewidth=4,
                    capthick=4,
                    linestyle='--',
                    zorder=1
                )

            if len(df_class) > 0:
                B_mid = get_mid_B(df_class)

                df_mid = df_class[
                    (np.isclose(df_class['Magnetic_Field'], B_mid)) &
                    (np.isclose(df_class['Velocity'], 20.0))
                ]

                if len(df_mid) > 0:
                    y_mid_30 = df_mid['Ratio_Above_30GeV'].iloc[0]
                    y_mid_all = df_mid['Ratio_All_Energy'].iloc[0]

                ax.plot(
                    i, y_mid_30,
                    marker='o',
                    markersize=12,
                    markerfacecolor=color_mapping[sc],
                    markeredgecolor='black',
                    markeredgewidth=0.8,
                    alpha=1.0,
                    linestyle='None',
                    zorder=5
                )
                
                ax.plot(
                    i, y_mid_all,
                    marker='o',
                    markersize=12,
                    markerfacecolor=color_mapping[sc],
                    markeredgecolor='black',
                    markeredgewidth=0.8,
                    alpha=0.3,
                    linestyle='None',
                    zorder=5
                )

        ax.set_ylim(0.0, 4.0)
        ax.set_yticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels(labels)

        ax.set_ylabel('Modulation Ratio')
        ax.set_xlabel('Spectral Class')
        label_s = 2.7 if np.isclose(s, 0.0) else s
        ax.set_title(f'Modulation by Spectral Class (s = {label_s})')

        ax.grid(True, axis='y', linestyle='--')
        ax.grid(False, axis='x')

        legend_lines = [
            Line2D([0], [0], color='black', lw=4, alpha=1.0, label='E > 30 GeV'),
            Line2D([0], [0], color='black', lw=4, alpha=0.3, linestyle='--', label='E > 10 MeV'),
            Line2D([0], [0], marker='o', color='black', linestyle='None',
                   markersize=7, markerfacecolor='black',
                   label='Ref (v=20, B_mid)')
        ]
        ax.legend(
            handles=legend_lines,
            fontsize=13,
            loc='upper left',

            handlelength=1.5,
            handletextpad=0.6,
            borderpad=0.6
        )

        plt.show()
