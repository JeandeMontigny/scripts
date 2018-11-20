#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
# parse output_file to construct tabs with average RI and death rate for each simulation parameter
def main(output_file):

    param=0; ri_temps=[]; death_temps=[]
    # main tabs
    ri_mean=[]; ri_std=[]; death=[]
    fichier=open(output_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        # regular expression for line parsing
        # info order: param ri0 ri1 ri2 cellDensity cellDeath migrationDistance
        m=re.search( r'(.+) (.+) (.+) (.+) (.+) (.+) (.+)', line, re.M|re.I)
        if m:
            # if simulation parameter for this result is not the same as before
            if (float(m.group(1)) != param):
                # if temps tabs are not empty, to avoid adding to main tab empty cells
                if not len(ri_temps) == 0:
                    # append average simulation results for this parameter to the main tab
                    ri_mean.append(np.average(ri_temps))
                    ri_std.append(np.std(ri_temps))
                    death.append(np.average(death_temps))
                # empty temps tab for next simulation parameter results
                ri_temps=[]
                death_temps=[]
                param = float(m.group(1))
            # append this simulation results to temps tabs
            ri_temps.append((float(m.group(2))+float(m.group(3))+float(m.group(4)))/3)
            death_temps.append(float(m.group(6)))

    # need to append the last param values to main tabs
    ri_mean.append(np.average(ri_temps))
    ri_std.append(np.std(ri_temps))
    death.append(np.average(death_temps))

    # print corresponding index for specified death value
#    for i, value in enumerate(death):
#        if value > 79 and value < 80:
#            print(i)

    # call figure construction
    figure(ri_mean, ri_std, death)

#--------------------------------------------------------------------------#
# simple plot of death vs RI (with std)
def figure(ri, ri_std, deathRate):

    plt.figure(1)
    plt.errorbar(deathRate, ri, ri_std, color='black', ecolor='gray')
    plt.axhline(3, color='silver')
    plt.xlabel("death rate")
    plt.ylabel("regularity index")
    plt.title("cell death impact on RI")

    plt.show()

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [death_RI file]")

