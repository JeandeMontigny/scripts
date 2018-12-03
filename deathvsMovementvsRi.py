#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
# parse output_file to construct tabs with average RI and death rate for each simulation parameter
def main(output_file):

    moveParam=0; ri_temps=[]; death_temps=[]
    # main tabs
    ri_mean=[]; death=[]; move_param=[]
    fichier=open(output_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        # regular expression for line parsing
        # info order: movementParam deathParam ri0 ri1 ri2 cellDensity cellDeath migrationDistance
        m=re.search( r'(.+) (.+) (.+) (.+) (.+) (.+) (.+) (.+)', line, re.M|re.I)
        if m:
            # if simulation parameter for this result is not the same as before
            if (float(m.group(1)) != moveParam):
                # if temps tabs are not empty, to avoid adding to main tab empty cells
                if not len(ri_temps) == 0:
                    # append average simulation results for this parameter to the main tab
                    ri_mean.append(np.average(ri_temps))
                    death.append(np.average(death_temps))
                    move_param.append(float(m.group(1)))
                # empty temps tab for next simulation parameter results
                ri_temps=[]
                death_temps=[]
                moveParam = float(m.group(1))
            # append this simulation results to temps tabs
            ri_temps.append((float(m.group(3))+float(m.group(4))+float(m.group(5)))/3)
            death_temps.append(float(m.group(7)))

    # need to append the last param values to main tabs
    ri_mean.append(np.average(ri_temps))
    death.append(np.average(death_temps))
    move_param.append(float(m.group(1)))

    # call figure construction
    figure(ri_mean, death, move_param)

#--------------------------------------------------------------------------#
# simple plot of death vs RI (with std)
def figure(ri, deathRate, move_param):

    plt.figure(1)
    # scatter plot with RI value as colour
    plt.scatter(deathRate, move_param, c=ri, marker='s', s=64)
    # plot colour bar and legend
    plt.colorbar().set_label('Regularity index')
    # vertical line for good death rate
    plt.axvline(x=60, color='silver', linestyle='--')
    plt.axvline(x=80, color='silver', linestyle='--')
    # legend
    plt.xlabel("Death rate")
    plt.ylabel("Movement threshold")
    plt.title("Cell death and cell migration impact on RI")

    plt.show()

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [death_movement_RI file]")

