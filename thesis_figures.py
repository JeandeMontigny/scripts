#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main():

    dev_days = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    rgc_nb = [6978.94, 5762.4, 4465.67, 4061.95, 3754.3, 3597.26, 3544.54, 3538.93, 3537.47]
    rgc_nb_std = [270.37, 355.08, 144.38, 197.28, 187.25, 199.27, 201.43, 239.88, 132.78]
    rgc_ri = [4.81, 5.12, 5.45, 5.49, 5.17, 5.28, 5.62, 6.03, 6.1]
    rgc_ri_std = [0.57, 0.61, 0.58, 0.35, 0.48, 0.5, 0.69, 0.43, 0.48]

    surface = [11.36, 12.48, 13.08, 13.66, 13.9, 14.4, 14.88, 15.16, 15.05]
    surface_std = [0.73, 2.24, 1.77, 1.83, 1.91, 1.8, 1.65, 2.4, 1.75]

    corrected_rgc_nb = [76905.17, 59836.19, 50788.67, 48172.12, 47281.57, 47635.66, 47927.51, 48095.77, 48527.03]
    corrected_rgc_std = [2705.69, 5261.80, 2497.52, 1684.67, 2145.28, 1812.71, 2644.79, 4136.08, 2696.44]

    death = [0.0, 45.83, 82.58, 92.44, 99.55, 99.20, 98.49, 100.0]
    death_std = [31.34, 30.81, 14.75, 13.12, 13.58, 15.85, 25.58, 25.74]

    sac_gcl_nb = [1114.1, 1078.46, 1098.73, 1063.77, 1094.6, 1074.22, 1087.29]
    sac_gcl_nb_std = [104.56, 110.78, 121.4, 132.62, 110.41, 108.77, 109.27]
    sac_gcl_ri = [3.57, 3.64, 3.78, 3.76, 3.58, 3.63, 3.53]
    sac_gcl_ri_std = [0.44, 0.54, 0.62, 0.36, 0.35, 0.36, 0.46]

    sac_inl_nb = [1345.3, 1292.37, 1334.01, 1357.86, 1297.15, 1359.61, 1285.63]
    sac_inl_nb_std = [141.99, 232.16, 167.14, 131.31, 117.52, 157.27, 125.79]
    sac_inl_ri = [4.18, 4.12, 4.35, 4.33, 4.41, 4.41, 4.36]
    sac_inl_ri_std = [0.48, 0.54, 0.44, 0.37, 0.39, 0.46, 0.66]

    # exclusion diameter of 32um from cell center
    sac_exclusion_factor = [0.772, 0.759, 0.768, 0.765, 0.748, 0.761, 0.759]
    sac_exclusion_factor_std = [0.13, 0.12, 0.09, 0.12, 0.10, 0.08, 0.09]

    real_death = [0.0, 45.83, 82.58, 92.44, 99.55, 99.20, 98.49, 100.0]
    real_death_std = [31.34, 30.81, 14.75, 13.12, 13.58, 15.85, 25.58, 25.74]

    # real vs simulated population
    pop_real = [135000/1350, 76905.17/1350, 59836.19/1350, 50788.67/1350, 48172.12/1350, 47281.57/1350, 47635.66/1350, 47927.51/1350, 48095.77/1350, 48527.03/1350]
    pop_real_std = [2624.75/1350, 2705.69/1350, 5261.80/1350, 2497.52/1350, 1684.67/1350, 2145.28/1350, 1812.71/1350, 2644.79/1350, 4136.08/1350, 2696.44/1350]
    pop_simu = [100, 62.0, 46.1, 39.8, 36.5, 35.6, 35.6, 35.6, 35.6, 35.6]
    pop_simu_std = [0, 4.8, 3.4, 3.2, 2.8, 1.6, 1.4, 0.8, 0.8, 0.8]

    # ri vs death rate high density
    high_dens_death = [1.43, 1.82, 1.83, 2.07, 2.19, 2.32, 2.52, 2.94, 3.1, 3.22, 3.34, 3.53, 4.14, 4.5, 5.1, 5.4, 6.11, 6.36, 6.97, 7.51, 8.11, 7.98, 8.13, 8.43, 8.58, 9.18, 9.25, 9.86, 10.0, 9.54, 10.06, 9.76, 9.86, 10.0, 9.83, 9.94, 9.9, 9.66, 9.94, 9.95]
    high_dens_death_std = [0.35, 0.24, 0.24, 0.31, 0.12, 0.4, 0.4, 0.16, 0.3, 0.29, 0.15, 0.27, 0.38, 0.36, 0.38, 0.27, 0.22, 0.53, 0.21, 0.53, 0.66, 0.64, 0.66, 0.64, 0.84, 0.91, 0.33, 0.38, 0.54, 0.63, 0.41, 0.16, 0.6, 0.74, 0.42, 0.78, 0.45, 0.37, 0.54, 0.42]
    high_dens_death_rate = [98.6, 98.19, 97.2, 95.33, 93.61, 92.06, 89.38, 87.45, 83.86, 80.59, 77.44, 72.24, 68.48, 65.06, 61.62, 59.4, 55.87, 53.39, 50.47, 47.14, 45.45, 42.64, 40.51, 38.79, 36.22, 33.48, 31.17, 29.42, 27.88, 25.42, 24.52, 23.64, 19.12, 16.96, 15.73, 13.78, 12.84, 11.76, 9.6, 6.28]
    # ri vs death rate low density
    low_dens_death = [1.93, 2.13, 2.32, 2.65, 2.88, 3.05, 3.61, 3.34, 3.97, 4.21, 4.39, 4.29, 4.9, 5.09, 5.28, 5.42, 5.29, 5.33, 5.29, 5.17, 5.2, 5.23]
    low_dens_death_std = [0.62, 0.73, 0.63, 0.68, 0.62, 0.65, 0.59, 0.54, 0.27, 0.44, 0.58, 0.51, 0.83, 0.44, 0.82, 0.64, 0.27, 0.63, 0.67, 0.51, 0.98, 0.86]
    low_dens_death_rate = [95.76, 94.88, 90.64, 84.94, 80.99, 71.49, 64.47, 57.6, 50.29, 38.6, 36.99, 34.8, 28.95, 26.61, 21.35, 20.61, 18.71, 17.84, 16.67, 14.77, 8.92, 7.31]

    # mosaic method: delau and ri
    randon_weight = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    delaunay_x = [[2.63, 2.7, 2.76, 2.83, 2.89, 2.96, 3.02, 3.09, 3.15, 3.22, 3.28, 3.35, 3.41, 3.48, 3.54, 3.61, 3.67, 3.74, 3.8, 3.87, 3.93, 4.0, 4.06, 4.13, 4.19, 4.26, 4.32, 4.39, 4.45, 4.52, 4.58, 4.65, 4.71, 4.78, 4.84, 4.91, 4.97, 5.04, 5.1, 5.17, 5.23, 5.3, 5.36, 5.43, 5.49, 5.56, 5.62, 5.69, 5.75, 5.82], [2.19, 2.26, 2.32, 2.39, 2.46, 2.52, 2.59, 2.66, 2.72, 2.79, 2.86, 2.93, 2.99, 3.06, 3.13, 3.2, 3.27, 3.33, 3.4, 3.47, 3.53, 3.6, 3.67, 3.73, 3.8, 3.87, 3.94, 4.0, 4.07, 4.14, 4.2, 4.27, 4.34, 4.4, 4.47, 4.54, 4.61, 4.67, 4.74, 4.81, 4.87, 4.94, 5.01, 5.07, 5.14, 5.21, 5.28, 5.34, 5.41, 5.48], [1.77, 1.85, 1.93, 2.01, 2.09, 2.17, 2.25, 2.32, 2.4, 2.48, 2.56, 2.64, 2.72, 2.8, 2.88, 2.95, 3.03, 3.11, 3.19, 3.27, 3.35, 3.43, 3.51, 3.58, 3.66, 3.74, 3.82, 3.9, 3.98, 4.06, 4.14, 4.21, 4.29, 4.37, 4.45, 4.53, 4.61, 4.69, 4.77, 4.84, 4.92, 5.0, 5.08, 5.16, 5.24, 5.32, 5.4, 5.47, 5.55, 5.63], [1.34, 1.42, 1.5, 1.58, 1.67, 1.75, 1.83, 1.91, 2.0, 2.08, 2.16, 2.24, 2.32, 2.41, 2.49, 2.57, 2.65, 2.73, 2.82, 2.9, 2.98, 3.06, 3.14, 3.23, 3.31, 3.39, 3.47, 3.55, 3.64, 3.72, 3.8, 3.88, 3.97, 4.05, 4.13, 4.21, 4.29, 4.38, 4.46, 4.54, 4.62, 4.7, 4.79, 4.87, 4.95, 5.03, 5.11, 5.2, 5.28, 5.36], [0.69, 0.79, 0.89, 0.99, 1.09, 1.19, 1.29, 1.39, 1.49, 1.59, 1.69, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.91, 3.01, 3.11, 3.21, 3.31, 3.41, 3.51, 3.61, 3.71, 3.81, 3.92, 4.02, 4.12, 4.22, 4.32, 4.42, 4.52, 4.62, 4.72, 4.82, 4.92, 5.03, 5.13, 5.23, 5.33, 5.43, 5.53, 5.63], [0.55, 0.66, 0.76, 0.87, 0.98, 1.08, 1.19, 1.3, 1.4, 1.51, 1.62, 1.72, 1.83, 1.94, 2.04, 2.15, 2.26, 2.36, 2.47, 2.58, 2.68, 2.79, 2.9, 3.0, 3.11, 3.22, 3.32, 3.43, 3.53, 3.64, 3.75, 3.85, 3.96, 4.07, 4.17, 4.28, 4.39, 4.49, 4.6, 4.71, 4.81, 4.92, 5.03, 5.13, 5.24, 5.35, 5.45, 5.56, 5.67, 5.77], [0.06, 0.18, 0.3, 0.42, 0.55, 0.67, 0.79, 0.92, 1.04, 1.16, 1.28, 1.41, 1.53, 1.65, 1.78, 1.9, 2.02, 2.14, 2.27, 2.39, 2.51, 2.64, 2.76, 2.88, 3.0, 3.13, 3.25, 3.37, 3.5, 3.62, 3.74, 3.86, 3.99, 4.11, 4.23, 4.36, 4.48, 4.6, 4.72, 4.85, 4.97, 5.09, 5.21, 5.34, 5.46, 5.58, 5.71, 5.83, 5.95, 6.07], [0.2, 0.33, 0.45, 0.58, 0.7, 0.82, 0.95, 1.07, 1.19, 1.32, 1.44, 1.56, 1.69, 1.81, 1.93, 2.06, 2.18, 2.3, 2.43, 2.55, 2.68, 2.8, 2.92, 3.05, 3.17, 3.29, 3.42, 3.54, 3.66, 3.79, 3.91, 4.03, 4.16, 4.28, 4.4, 4.53, 4.65, 4.78, 4.9, 5.02, 5.15, 5.27, 5.39, 5.52, 5.64, 5.76, 5.89, 6.01, 6.13, 6.26], [0.15, 0.27, 0.4, 0.52, 0.65, 0.77, 0.9, 1.02, 1.15, 1.27, 1.4, 1.53, 1.65, 1.78, 1.9, 2.03, 2.15, 2.28, 2.4, 2.53, 2.65, 2.78, 2.9, 3.03, 3.15, 3.28, 3.41, 3.53, 3.66, 3.78, 3.91, 4.03, 4.16, 4.28, 4.41, 4.53, 4.66, 4.78, 4.91, 5.03, 5.16, 5.29, 5.41, 5.54, 5.66, 5.79, 5.91, 6.04, 6.16, 6.29]]
    delaunay_cumuls = [[0.02, 0.05, 0.1, 0.18, 0.26, 0.37, 0.51, 0.65, 0.79, 0.92, 1.04, 1.14, 1.23, 1.31, 1.37, 1.41, 1.45, 1.48, 1.51, 1.53, 1.56, 1.59, 1.63, 1.67, 1.73, 1.79, 1.84, 1.88, 1.91, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 2.0, 2.03, 2.06, 2.09, 2.13, 2.16, 2.19, 2.22, 2.24, 2.26, 2.28, 2.29, 2.29, 2.3, 2.3], [0.01, 0.02, 0.03, 0.06, 0.09, 0.13, 0.16, 0.21, 0.27, 0.32, 0.37, 0.44, 0.51, 0.6, 0.67, 0.76, 0.83, 0.9, 0.97, 1.03, 1.09, 1.14, 1.2, 1.26, 1.31, 1.37, 1.42, 1.46, 1.51, 1.56, 1.6, 1.65, 1.69, 1.72, 1.76, 1.79, 1.82, 1.85, 1.87, 1.88, 1.9, 1.91, 1.91, 1.92, 1.93, 1.93, 1.94, 1.94, 1.94, 1.94], [0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.09, 0.11, 0.14, 0.17, 0.19, 0.23, 0.26, 0.3, 0.34, 0.39, 0.44, 0.5, 0.55, 0.62, 0.68, 0.74, 0.81, 0.87, 0.93, 0.99, 1.05, 1.11, 1.17, 1.22, 1.27, 1.32, 1.37, 1.41, 1.45, 1.49, 1.52, 1.55, 1.57, 1.59, 1.61, 1.63, 1.64, 1.65, 1.66, 1.66, 1.66, 1.67, 1.67, 1.67], [0.0, 0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.08, 0.1, 0.12, 0.14, 0.17, 0.2, 0.23, 0.26, 0.3, 0.34, 0.39, 0.43, 0.48, 0.53, 0.58, 0.64, 0.69, 0.74, 0.79, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.14, 1.18, 1.22, 1.26, 1.29, 1.33, 1.35, 1.38, 1.4, 1.41, 1.43, 1.44, 1.45, 1.45, 1.46, 1.47, 1.47], [0.0, 0.0, 0.01, 0.01, 0.02, 0.03, 0.03, 0.04, 0.05, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.19, 0.22, 0.25, 0.28, 0.32, 0.36, 0.39, 0.44, 0.48, 0.53, 0.58, 0.62, 0.67, 0.71, 0.76, 0.8, 0.85, 0.89, 0.93, 0.97, 1.01, 1.04, 1.07, 1.1, 1.13, 1.16, 1.19, 1.21, 1.22, 1.24, 1.25, 1.26, 1.28, 1.28, 1.29], [0.0, 0.0, 0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.04, 0.06, 0.07, 0.09, 0.1, 0.12, 0.14, 0.16, 0.19, 0.21, 0.24, 0.27, 0.3, 0.34, 0.37, 0.41, 0.45, 0.49, 0.53, 0.57, 0.6, 0.65, 0.69, 0.73, 0.77, 0.81, 0.84, 0.88, 0.91, 0.94, 0.97, 1.0, 1.03, 1.05, 1.07, 1.09, 1.1, 1.12, 1.13, 1.14, 1.15, 1.16], [0.0, 0.0, 0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.14, 0.16, 0.19, 0.21, 0.24, 0.27, 0.3, 0.33, 0.36, 0.4, 0.43, 0.47, 0.51, 0.55, 0.59, 0.63, 0.66, 0.7, 0.73, 0.77, 0.8, 0.83, 0.86, 0.89, 0.91, 0.93, 0.95, 0.97, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05], [0.0, 0.01, 0.01, 0.01, 0.02, 0.03, 0.03, 0.04, 0.05, 0.06, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.25, 0.28, 0.31, 0.34, 0.37, 0.41, 0.44, 0.48, 0.51, 0.55, 0.58, 0.61, 0.65, 0.68, 0.71, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.97, 0.98], [0.0, 0.01, 0.01, 0.02, 0.02, 0.03, 0.04, 0.05, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.15, 0.17, 0.19, 0.22, 0.25, 0.27, 0.3, 0.33, 0.36, 0.39, 0.43, 0.46, 0.5, 0.53, 0.56, 0.59, 0.63, 0.66, 0.69, 0.72, 0.75, 0.77, 0.8, 0.82, 0.84, 0.86, 0.88, 0.89, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.97]]
    ri_randon_weight = [23.69125, 11.14, 6.71, 4.99, 3.80, 3.12, 2.58, 2.23, 2.09]
    ri_randon_weight_std = [0.98, 0.41, 0.20, 0.17, 0.14, 0.14, 0.08, 0.13, 0.07]

    mig_denstiy = [249.19, 125.55, 205.32, 46.79, 96.96, 35.08, 54.04, 58.98, 21.48, 113.37, 106.94, 43.77, 21.81, 75.87, 132.02, 153.04, 51.15, 50.66, 77.55, 345.2, 52.06, 46.4, 96.23, 44.12, 48.46, 24.61, 56.21, 53.46, 45.97, 45.76, 54.05, 48.16, 48.75, 50.19, 23.34, 19.07, 24.92, 20.05, 48.97, 33.56, 21.29, 19.75, 22.39]
    mig_denstiy_std = [12.66, 10.69, 5.46, 3.74, 5.41, 2.86, 5.46, 4.23, 1.8, 7.68, 6.75, 3.25, 3.46, 3.88, 5.03, 7.5, 6.38, 2.87, 8.32, 11.1, 6.74, 5.11, 2.43, 3.48, 2.48, 3.24, 1.55, 2.55, 2.54, 4.34, 4.74, 3.77, 2.02, 3.02, 4.01, 1.13, 3.85, 2.76, 4.83, 2.81, 4.07, 2.27, 1.11]
    mig_dist = [8.8, 4.97, 6.9, 9.03, 4.48, 8.32, 8.44, 6.99, 7.79, 4.98, 3.99, 7.96, 7.36, 5.73, 5.32, 5.47, 9.4, 8.4, 6.33, 13.39, 8.22, 6.9, 4.82, 6.66, 10.15, 7.7, 9.19, 7.39, 6.62, 8.59, 7.32, 7.52, 9.21, 10.09, 7.95, 7.35, 8.77, 8.87, 8.49, 7.79, 6.91, 6.46, 8.15]
    mig_dist_std = [9.34, 7.17, 8.53, 9.93, 7.19, 9.65, 10.29, 8.95, 9.33, 7.45, 6.02, 9.23, 8.52, 8.45, 7.56, 7.43, 10.38, 9.85, 9.18, 11.23, 9.66, 8.67, 7.86, 9.29, 10.94, 8.24, 11.3, 8.88, 7.99, 9.74, 8.8, 8.63, 10.4, 11.23, 9.18, 7.53, 9.28, 9.87, 9.98, 9.76, 7.8, 7.11, 9.11]


    # cluster position OD to periphery ratio
    centre_peri_p2 = [0.214, 0.28, 0.33, 0.219, 0.25, 0.2, 0.32, 0.189, 0.246, 0.299, 0.378, 0.415, 0.461, 0.338, 0.342, 0.342, 0.365, 0.367, 0.239, 0.357, 0.199, 0.166]
    centre_peri_p3 = [0.284, 0.477, 0.437, 0.353, 0.384, 0.3812, 0.328, 0.442, 0.405, 0.395, 0.41, 0.454, 0.451, 0.315, 0.389, 0.453, 0.357, 0.238, 0.22, 0.273, 0.235, 0.327, 0.328, 0.357, 0.295, 0.333, 0.217, 0.415]
    centre_peri_p4 = [0.552, 0.533, 0.467, 0.384, 0.448, 0.43, 0.44, 0.42, 0.506, 0.661, 0.664, 0.67, 0.683, 0.617, 0.644, 0.62, 0.626, 0.779, 0.585, 0.784, 0.725, 0.728, 0.37, 0.412, 0.454, 0.303, 0.376, 0.493, 0.329, 0.304, 0.526, 0.6, 0.565, 0.569, 0.511, 0.443, 0.473, 0.536, 0.569, 0.555, 0.176, 0.27, 0.547]
    centre_peri_p5 = [0.691, 0.671, 0.56, 0.667, 0.786, 0.666, 0.591, 0.741, 0.757, 0.655, 0.552, 0.665, 0.698, 0.675, 0.636, 0.467, 0.542, 0.437, 0.751, 0.684, 0.728, 0.828, 0.883, 0.559, 0.604, 0.422, 0.516, 0.598, 0.441, 0.576, 0.805, 0.577, 0.659, 0.752]
    centre_peri_p6 = [0.79, 0.657, 0.854, 0.761, 0.738, 0.848, 0.685, 0.737, 0.768, 0.907, 0.876, 0.881, 0.959, 0.877, 0.85, 0.73, 0.814, 0.784, 0.868, 0.819, 0.756, 0.776, 0.78, 0.79, 0.95, 0.884, 0.72, 0.855, 0.927, 0.72, 0.651, 0.77, 0.637, 0.576, 0.56, 0.653, 0.75, 0.76, 0.87]
    centre_peri_p7 = [0.871, 0.908, 0.735, 0.686, 0.815, 0.89, 0.565, 0.91, 0.93, 0.826, 0.81, 0.81, 0.788, 0.819, 0.822, 0.77, 0.788, 0.95, 0.828, 0.853, 0.744, 0.87, 0.9, 0.786, 0.959, 0.96, 0.53, 0.767, 0.77, 0.78, 0.894, 0.797, 0.944, 0.787, 0.752, 0.856, 0.695, 0.598, 0.779]
    centre_peri_p8 = [0.768, 0.945, 0.94, 0.92, 0.737, 0.735, 0.76, 0.877, 0.921, 0.82, 0.86, 0.83, 0.77]
    centre_peri_p9 = [0.83, 0.84, 0.83, 0.92, 0.82, 0.89, 0.93, 0.85, 0.9, 0.85, 0.96, 0.92, 0.96, 0.91, 0.96]

    wave_origin_p2 = [1.29, 1.25, 1.71]
    wave_origin_p3 = [1.05, 1.73, 2.09, 2.17]
    wave_origin_p4 = [1.8, 1.24, 4.6, 3.45, 1.05]
    wave_origin_p5 = [3.9, 2.94, 5.3, 4.7, 3.71, 1.45]
    wave_origin_p6_7 = [2.9, 4, 6, 3.55]
    wave_origin_p8_9 = [3.1, 3.4, 2.8]
    wave_origin_p10_13 = [1.45, 3.35, 0.58, 0.98]

    # -------- front sizes -------- #
    ticks_size = 11; label_size = 12

    # -------- Figure 2 -------- #
    plt.figure(figsize=(10, 5.12))
    plot = plt.subplot(1, 2, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(-7, 1.08, "A", size = 16)
    for i in range(0, len(delaunay_cumuls)):
        plt.plot([val*10 for val in delaunay_x[i]],\
            [val/max(delaunay_cumuls[i]) for val in delaunay_cumuls[i]],\
            label = randon_weight[i])
    plt.xlabel("Segment length", size = label_size)
    plt.ylabel("Cumulative probability", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)
    plt.legend()

    plot = plt.subplot(1, 2, 2)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(0, 26.5, "B", size = 16)
    plt.errorbar(randon_weight, ri_randon_weight, ri_randon_weight_std)
    plt.xlabel("Random weight", size = label_size)
    plt.ylabel("Regularity index", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plt.tight_layout()

    # -------- Figure 3 -------- #
    plt.figure(figsize=(10, 12))

    plot = plt.subplot(2, 2, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(0.8, 7.5e3, "A", size = 16)
    plt.errorbar(dev_days, rgc_nb, rgc_nb_std)
    plt.xlim(1.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Cell density (cells/mm$^2$)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plot = plt.subplot(2, 2, 2)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(0.8, 18, "B", size = 16)
    plt.errorbar(dev_days, surface, surface_std)
    plt.xlim(1.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Retinal surface (mm$^2$)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plot = plt.subplot(2, 2, 3)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(0.8, 8.2e4, "C", size = 16)
    plt.errorbar(dev_days, corrected_rgc_nb, corrected_rgc_std)
    plt.xlim(1.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Estimated cell population", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plot = plt.subplot(2, 2, 4)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(1.8, 135, "D", size = 16)
    plt.errorbar(dev_days[1:], death, death_std)
    plt.xlim(2.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Cumulative apoptosis (%)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plt.tight_layout()

    # -------- Figure 4 -------- #
    plt.figure(figsize=(10, 5.12))
    plot = plt.subplot(1, 2, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(2.9, 1.56e3, "A", size = 16)
    plt.errorbar(dev_days[2:], sac_gcl_nb, sac_gcl_nb_std, label = "GCL population")
    plt.errorbar(dev_days[2:], sac_inl_nb, sac_inl_nb_std, label = "INL population")
    plt.xlim(3.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Cumulative apoptosis (%)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)
    plt.legend()

    plot = plt.subplot(1, 2, 2)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(2.9, 5.15, "B", size = 16)
    plt.errorbar(dev_days[2:], sac_gcl_ri, sac_gcl_ri_std, label = "GCL population")
    plt.errorbar(dev_days[2:], sac_inl_ri, sac_inl_ri_std, label = "INL population")
    plt.xlim(3.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Cumulative apoptosis (%)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)
    plt.legend()

    plt.tight_layout()

    # -------- Figure 5 -------- #
    plt.figure(figsize=(5.12, 5.12))
    plot = plt.subplot(1, 1, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.errorbar(dev_days[2:], sac_exclusion_factor, sac_exclusion_factor_std)
    plt.xlim(3.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Exclusion factor", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plt.tight_layout()

    # -------- Figure 6 -------- #
    plt.figure(figsize=(5.12, 5.12))
    plot = plt.subplot(1, 1, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.errorbar(dev_days, rgc_ri, rgc_ri_std)
    plt.xlim(2.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Regularity index", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plt.tight_layout()

    # -------- Figure 8 -------- #
    plt.figure(figsize=(10, 12))
    plot = plt.subplot(3, 2, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.errorbar(dev_days, rgc_ri, rgc_ri_std)
    plt.xlim(2.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Regularity index", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plt.tight_layout()

    # -------- Figure 9 -------- #
    plt.figure(figsize=(10, 5.12))
    plot = plt.subplot(1, 2, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(-0.3, 106.5, "A", size = 16)
    plt.errorbar([x for x in range(1, 11)], pop_real, yerr = pop_real_std, label = "In-vitro")
    plt.errorbar([x for x in range(1, 11)], pop_simu, yerr = pop_simu_std, label = "In-silico")
    plt.xlim(0.5, 10.5)
    plt.xlabel("Postnatal day", size = label_size)
    plt.ylabel("Cell population (%)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)
    plt.legend()

    plot = plt.subplot(1, 2, 2)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(-6, 11.3, "B", size = 16)
    plt.errorbar(high_dens_death_rate, high_dens_death, high_dens_death_std, label = "High density")
    plt.errorbar(low_dens_death_rate, low_dens_death, low_dens_death_std, label = "Low density")
    plt.axvline(65, color='silver', linestyle='--')
    plt.xlabel("Death rate", size = label_size)
    plt.ylabel("Regularity index", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)
    plt.legend()

    plt.tight_layout()

    # -------- Figure 11 -------- #
    plt.figure(figsize=(10, 5.12))
    plot = plt.subplot(1, 2, 1)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(-5, 1800, "A", size = 16)
    plt.bar([6.42, 12.85, 19.28, 25.71, 32.14, 38.57, 45.0, 51.43, 57.86, 64.29], [1705, 607, 330, 174, 95, 46, 23, 10, 3, 1], width = 6)
    plt.xlabel("Distance (µm)", size = label_size)
    plt.ylabel("Number of cell", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plot = plt.subplot(1, 2, 2)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.text(-20, 26.2, "B", size = 16)
    plt.errorbar(mig_denstiy, mig_dist, xerr = mig_denstiy_std, yerr = mig_dist_std, fmt = 'o', color = "black", ecolor = "gray")
    plt.plot([106.94, 113.36, 125.55, 132.02, 153.03, 205.31, 249.18, 345.19], [4.11, 4.34, 4.79, 5.02, 5.79, 7.68, 9.28, 12.77], color = "red")
    plt.xlabel("Cell type density", size = label_size)
    plt.ylabel("Distance (µm)", size = label_size)
    plt.xticks(size = ticks_size)
    plt.yticks(size = ticks_size)

    plt.tight_layout()


    # -------- Clusters Figure 1 -------- #
    fig, ax1 = plt.subplots()
    ax1.spines['top'].set_visible(False)
    ax1.set_xlabel("Postnatal day", size = label_size)
    ax1.set_ylabel("Cluster D1/D2 ratio", size = label_size)
    ax1.boxplot([centre_peri_p2, centre_peri_p3, centre_peri_p4, centre_peri_p5, centre_peri_p6, centre_peri_p7, centre_peri_p8, centre_peri_p9], labels=dev_days[:len(dev_days)-1], medianprops=dict(color='black'))
    ax1.text(2.5, 0.82, "_______", horizontalalignment='center')
    ax1.text(2.5, 0.82, "***", horizontalalignment='center')
    ax1.text(3.5, 0.92, "_______", horizontalalignment='center')
    ax1.text(3.5, 0.92, "***", horizontalalignment='center')
    ax1.text(4.5, 0.99, "_______", horizontalalignment='center')
    ax1.text(4.5, 0.99, "***", horizontalalignment='center')
    ax2 = ax1.twinx()
    ax2.spines['top'].set_visible(False)
    ax2.set_ylabel("% difference", color='red', size = label_size)
    ax2.plot(dev_days[0:len(dev_days)-2], [33.12, 63.71, 24.59, 25.76, 3.523, 4.863, 6.848], 'r--')
    ax2.set_ylim(0, 100)

    plt.tight_layout()

    # -------- Clusters Figure 2 -------- #
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Postnatal day", size = label_size)
    ax.set_ylabel("Wave origin periphery/centre ratio", size = label_size)
    ax.boxplot([wave_origin_p2, wave_origin_p3, wave_origin_p4, wave_origin_p5, wave_origin_p6_7, wave_origin_p8_9, wave_origin_p10_13], labels=["2", "3", "4", "5", "6-7", "8-9", "10-13"], medianprops=dict(color='black'))

    plt.tight_layout()

#--------------------------------------------------------------------------#

main()
plt.show()
