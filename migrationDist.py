#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import factorial

#--------------------------------------------------------------------------#
def main(distance_file):
    distance=[]
    fichier=open(distance_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        m=re.search( r'(.+)', line, re.M|re.I)
        if m:
            distance.append(float(m.group(1)))

    print("average tangential distance:", round(np.average(distance),2), "micrometers ; with std of", round(np.std(distance),2))

    figure(distance)

#--------------------------------------------------------------------------#
def figure(distance):
    # the bins should be of integer width, because poisson is an integer distribution
    entries, bin_edges, patches = plt.hist(distance, normed=True)

    # calculate binmiddles
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])

    # fit with curve_fit
    parameters, cov_matrix = curve_fit(poisson, bin_middles, entries) 

    # plot poisson-deviation with fitted parameter
    x_plot = np.linspace(0, max(distance))

    plt.plot(x_plot, poisson(x_plot, *parameters), 'r-', lw=2)

    plt.xlabel("distance")
    plt.ylabel("number of cell")
    plt.title("distribution of tangential migration distance")

    plt.show()

#--------------------------------------------------------------------------#
# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [migration distance file]")
