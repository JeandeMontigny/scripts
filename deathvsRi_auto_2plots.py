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
        m=re.search( r'.+ (.+) (.+) (.+) (.+) (.+) (.+) (.+)', line, re.M|re.I)
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

    # return data for figure construction
    return(ri_mean, ri_std, death)

#--------------------------------------------------------------------------#
# simple plot of death vs RI (with std)
def figure(data1, data2):
    
    ri1 = data1[0]; ri_std1 = data1[1]; deathRate1 = data1[2]
    ri2 = data2[0]; ri_std2 = data2[1]; deathRate2 = data2[2]

    plt.figure(1)
    plt.errorbar(deathRate1, ri1, ri_std1, color='blue', ecolor='lightsteelblue', label="CD")
    plt.errorbar(deathRate2, ri2, ri_std2, color='red', ecolor='lightcoral', label="CF & CD")

    plt.legend(loc=2)

    plt.axhline(3, color='silver')
    plt.axhline(1.8, color='silver', linestyle='--')
    plt.ylim(1.5)
    plt.xlabel("death rate")
    plt.ylabel("regularity index")
    plt.title("cell death impact on RI")

    plt.show()

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==3:
    data1 = main(sys.argv[1])
    data2 = main(sys.argv[2])

    # figure construction
    figure(data1, data2)

    print("done")

else:
    exit("arg error - need 1 arg: [death_RI file] [death_RI file]")

