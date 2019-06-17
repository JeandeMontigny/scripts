#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
# parse output_file to construct tabs with average RI and death rate for each simulation parameter
def main(output_file):

    density=0; ri_temps=[]
    # main tabs
    density_list=[]; ri_mean=[]; ri_std=[]
    fichier=open(output_file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # regular expression for line parsing
        # info order: density ri0
        m=re.search( r'(.+) (.+)', line, re.M|re.I)
        if m:
            # if simulation parameter for this result is not the same as before
            if (float(m.group(1)) != density):
                # if temps tabs are not empty, to avoid adding to main tab empty cells
                if not len(ri_temps) == 0:
                    # append average simulation results for this parameter to the main tab
                    ri_mean.append(np.average(ri_temps))
                    ri_std.append(np.std(ri_temps))
                    density_list.append(m.group(1))
                # empty temps tab for next simulation parameter results
                ri_temps=[]
                density = float(m.group(1))
            # append this simulation results to temps tabs
            ri_temps.append(float(m.group(2)))

    # need to append the last param values to main tabs
    ri_mean.append(np.average(ri_temps))
    ri_std.append(np.std(ri_temps))
    density_list.append(m.group(1))

    # call figure construction
    figure(density_list, ri_mean, ri_std)

    return 1
    
#--------------------------------------------------------------------------#
# simple plot of death vs RI (with std)
def figure(density, ri_mean, ri_std):

    plt.figure(1)
    plt.errorbar(density, ri_mean, ri_std, color='black', ecolor='gray')
    plt.ylim(1.5, 11)
    plt.xlabel("density")
    plt.ylabel("regularity index")
    plt.title("cell density impact on RI")

    plt.show()

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    exit("arg error - need 1 arg: [param_RI_study_density.txt file]")

