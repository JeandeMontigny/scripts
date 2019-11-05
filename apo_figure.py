#!/usr/bin/env python3
import sys, os, re, warnings
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main():

    dev_days = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    rgc_nb = [6978.94, 5762.4, 4465.67, 4061.95, 3754.3, 3597.26, 3544.54, 3538.93, 3537.47]
    rgc_nb_std = [270.37, 355.08, 144.38, 197.28, 187.25, 199.27, 201.43, 239.88, 132.78]
    rgc_ri = [4.81, 5.12, 5.45, 5.49, 5.17, 5.28, 5.62, 6.03, 6.1]
    rgc_ri_std = [0.57, 0.61, 0.58, 0.35, 0.48, 0.5, 0.69, 0.43, 0.48]

    surface = [11.36, 12.48, 13.08, 13.66, 13.9, 14.4, 14.88, 15.16, 15.05] # +, 15.67] for day 11
    surface_std = [0.73, 2.24, 1.77, 1.83, 1.91, 1.8, 1.65, 2.4, 1.75] # +, 1.47] for day 11

    corrected_rgc_nb = [76905.17, 59836.19, 50788.67, 48172.12, 47281.57, 47635.66, 47927.51, 48095.77, 48527.03]
    corrected_rgc_std = [2705.69, 5261.80, 2497.52, 1684.67, 2145.28, 1812.71, 2644.79, 4136.08, 2696.44]

    death = [0.0, 45.83, 82.58, 92.44693968212697, 99.55, 99.20, 98.49, 100.0]
    death_std = [31.34, 30.81, 14.75, 13.12, 13.58, 15.85, 25.58, 25.74]

    sac_gcl_nb = [np.nan, np.nan, np.nan, 1008.46, 1098.73, 1063.77, 1109.7, 1074.22, 1087.29]
    sac_gcl_nb_std = [np.nan, np.nan, np.nan, 90.96, 121.4, 132.62, 62.03, 108.77, 109.27]
    sac_gcl_ri = [np.nan, np.nan, np.nan, 3.64, 3.93, 3.76, 3.58, 3.63, 3.53]
    sac_gcl_ri_std = [np.nan, np.nan, np.nan, 0.54, 0.26, 0.36, 0.35, 0.36, 0.46]

    sac_inl_nb = [np.nan, np.nan, np.nan, 1216.24, 1334.01, 1357.86, 1268.23, 1359.61, 1285.63]
    sac_inl_nb_std = [np.nan, np.nan, np.nan, 226.22, 167.14, 131.31, 94.4, 157.27, 125.79]
    sac_inl_ri = [np.nan, np.nan, np.nan, 4.12, 4.35, 4.33, 4.51, 4.51, 4.36]
    sac_inl_ri_std = [np.nan, np.nan, np.nan, 0.51, 0.44, 0.37, 0.39, 0.46, 0.66]


#    rgc_thickness = [0.1184, 0.1129, 0.1053, 0.1058, 0.0893, 0.0936, 0.0950]
#    rgc_thickness_std = [0.0070, 0.0078, 0.0070, 0.0056, 0.0066, 0.0125, 0.0105]

#    rgc_diam = [7.36, 8.01, 9.86, 10.98, 11.88, 12.43, 12.43]
#    rgc_diam_std = [0.67, 0.73, 0.53, 0.58, 0.37, 0.20, 0.67]

    createFigure(dev_days, rgc_nb, rgc_nb_std, "Postnatal day", "Cell density", "RGC density during development")
    createFigure(dev_days, rgc_ri, rgc_ri_std, "Postnatal day", "RI", "RGC Regularity Index during development")

    createFigure(dev_days, surface, surface_std, "Postnatal day", "Retinal surface (mm$^2$)", "Retinal surface during development")
    createFigure(dev_days, corrected_rgc_nb, corrected_rgc_std, "Postnatal day", "Estimated cell population", "Estimated RGC population during development")
    createFigure(dev_days[:-1], death, death_std, "Postnatal day", "Cell death (%)", "RGC cumulative death during development")

    createFigure(dev_days, sac_gcl_nb, sac_gcl_nb_std, "Postnatal day", "Cell density", "GCL SAC density during development")
    createFigure(dev_days, sac_gcl_ri, sac_gcl_ri_std, "Postnatal day", "RI", "GCL SAC Regularity Index during development")

    createFigure(dev_days, sac_inl_nb, sac_inl_nb_std, "Postnatal day", "Cell density", "INL SAC density during development")
    createFigure(dev_days, sac_inl_ri, sac_inl_ri_std, "Postnatal day", "RI", "INL SAC Regularity Index during development")

    plt.show()

#--------------------------------------------------------------------------#
def createFigure(x, y, std, x_label, y_label, title):
    """ Create error bar figure
    @params:
        x data list
        y data list
        standard deviation list
        label for x axis (Str)
        label for y axis (Str)
        figure title (Str) """
    plt.figure()

    plt.errorbar(x, y, std)
    plt.xlim(1.5, 10.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

#--------------------------------------------------------------------------#

main()
