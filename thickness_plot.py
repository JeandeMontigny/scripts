#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def main():
    thickness = [0.1184, 0.1129, 0.1053, 0.1058, 0.0893, 0.0936, 0.0950]
    std_thickness = [0.0070, 0.0078, 0.0070, 0.0056, 0.0066, 0.0125, 0.0105]

    diam = [7.36, 8.01, 9.86, 10.98, 11.88, 12.43, 12.43]
    std_diam = [0.67, 0.73, 0.53, 0.58, 0.37, 0.20, 0.67]

    density = [6094.0, 4530.5, 4040.7, 3741.4, 3228.9, 3052.1, 3160.9]
    std_density = [105.6, 253.0, 80.3, 145.5, 215.8, 252.3, 119.6]

    ri = [np.average([2.25, 2.63, 2.57, 2.07, 2.47]), np.average([3.0, 2.67, 2.83, 3.33, 3.09]), np.average([3.32, 3.69, 2.84, 3.55, 3.77]), np.average([4.75, 3.27, 3.4, 3.15, 4.34]), np.average([3.34, 3.44, 5.19, 2.79, 3.64]), np.average([3.54, 3.37, 4.26, 4.26, 4.44]), np.average([5.12, 5.74, 4.71, 4.04, 4.82])]
    std_ri = [np.std([2.25, 2.63, 2.57, 2.07, 2.47]), np.std([3.0, 2.67, 2.83, 3.33, 3.09]), np.std([3.32, 3.69, 2.84, 3.55, 3.77]), np.std([4.75, 3.27, 3.4, 3.15, 4.34]), np.std([3.34, 3.44, 5.19, 2.79, 3.64]), np.std([3.54, 3.37, 4.26, 4.26, 4.44]), np.std([5.12, 5.74, 4.71, 4.04, 4.82])]

    label = [3, 4, 5, 6, 7, 8, 9]

    plt.figure(1)
    plt.errorbar(label, thickness, std_thickness)
    plt.xlabel("Postnatal day")
    plt.ylabel("Relative thickness")
    plt.title("RGC layer thickness during development")

    plt.figure(2)
    plt.errorbar(label, diam, std_diam)
    plt.xlabel("Postnatal day")
    plt.ylabel("Diameter (micrometer)")
    plt.title("RGC diameter during development")

    plt.figure(3)
    plt.errorbar(label, density, std_density)
    plt.xlabel("Postnatal day")
    plt.ylabel("Density (cells/$mm^{2}$)")
    plt.title("RGC density during development")

    plt.figure(4)
    plt.errorbar(label, ri, std_ri)
    plt.xlabel("Postnatal day")
    plt.ylabel("Regularity index)")
    plt.title("RGC mosaic regularity during development")


    plt.show()

 
main()
