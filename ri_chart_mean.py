#!/usr/bin/env python3
import sys, os, re, numpy
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

#--------------------------------------------------------------------------#
def main(folder):
    ri_mean = []; death_rate = []; names = []
    for name in os.listdir(folder):
        if name.endswith(".txt"):
            results = read(folder+"/"+name)
            ri_mean.append(results[0])
            death_rate.append(results[1])
            names.append(name)

    simu_nb = len(ri_mean); simu_len = len(ri_mean[0])
    for ri in ri_mean:
        if simu_len != len(ri):
            print("error - simulation length is not the same for every simulation")
            return 0

#    plotIndividualSimus(ri_mean, death_rate, names)

    x_axis = []
    mean_ri_mean = []; mean_ri_std = []
    mean_death_rate = []; mean_cell_density = []; std_death_rate = []
    for i in range(0, simu_len):
        x_axis.append(i)

        mean_ri_mean_temps = []
        mean_death_rate_temps = []
        for j in range(len(ri_mean)):
            mean_ri_mean_temps.append(ri_mean[j][i])
            mean_death_rate_temps.append(death_rate[j][i])

        mean_ri_mean.append( round(np.average(mean_ri_mean_temps), 3) )
        mean_ri_std.append( round(np.std(mean_ri_mean_temps), 3) )

        mean_cell_density.append( 12000 - (12000 * round(np.average(mean_death_rate_temps), 3) /100 )  )
        mean_death_rate.append( round(np.average(mean_death_rate_temps), 3) )
        std_death_rate.append( round(np.std(mean_death_rate_temps), 3) )

    figure(x_axis, mean_ri_mean, mean_ri_std, "regularity index")
    figure(x_axis, mean_death_rate, std_death_rate, "death rate")
    figure(x_axis, mean_cell_density, std_death_rate, "cell density")

    figureRiVsDeath(mean_ri_mean, mean_death_rate)

    plt.show()

    return 1

#--------------------------------------------------------------------------#
def read(output_file):
    tab_ave = []; tab_death = []
    fichier = open(output_file, "r")
    file_lines = fichier.readlines()
    for line in file_lines:
        # on off on-off death
        m = re.search( r'(.+) (.+) (.+) (.+)', line, re.M|re.I)
        if m:
            tab_ave.append(np.average( [float(m.group(1)), float(m.group(2)), float(m.group(3))] ))
            tab_death.append(float(m.group(4)))

    return(tab_ave, tab_death)

#--------------------------------------------------------------------------#
def figure(x, y, std, y_label):
    plt.figure()
    plt.errorbar(x, y, yerr=std, color='black', ecolor='black')
    plt.xlabel("modelling time")
    plt.xticks([])
    plt.ylabel(y_label)
    if y_label == "regularity index":
        plt.ylim(1, 7)
    if y_label == "cell density":
        plt.ylim(3000, 12500)

#--------------------------------------------------------------------------#
def plotIndividualSimus(ri_list, death_list, name_list):
    nb_row = int(np.sqrt(len(ri_list)))
    nb_col = nb_row+1

    plt.figure(98)
    for i in range (0, len(ri_list)):
        ri = ri_list[i]
        plt.subplot(nb_row, nb_col, i+1)
        plt.plot(range(0, len(ri)), ri)
        plt.ylim(1, 7)
        plt.title(name_list[i])
        plt.subplots_adjust(left=0.05, right=0.95, bottom = 0.05, hspace=0.5)

    plt.figure(99)
    for i in range (0, len(ri_list)):
        ri = ri_list[i]
        death = death_list[i]
        diff_ri = []
        for j in range(10, 72):
            diff_ri.append(ri[j]-ri[j-1])
        diff_ri = sig.savgol_filter(diff_ri, 7, 3)
        plt.subplot(nb_row, nb_col, i+1)
        plt.plot(death[10:72], diff_ri, color='black')
        plt.axhline(y=0, color='gray', linestyle='--')
        plt.title(name_list[i])
        plt.subplots_adjust(left=0.05, right=0.95, bottom = 0.05, hspace=0.5)

    plt.show()

#--------------------------------------------------------------------------#
def figureRiVsDeath(mean_ri, mean_death):
    diff_ri = []; diff_death = []

    for i in range(10, 72):
        diff_ri.append(mean_ri[i]-mean_ri[i-1])
    # smooth curve
    diff_ri = sig.savgol_filter(diff_ri, 7, 3)

    plt.figure()
    plt.plot(mean_death[10:72], diff_ri, color='black')
    plt.axhline(y=0, color='gray', linestyle='--')
    plt.xlabel("death rate")
    plt.ylabel("RI evolution")

#--------------------------------------------------------------------------#
if len(sys.argv) == 2:
    if (main(sys.argv[1])):
        print("done")
    else:
        print("erorr")
else:
    exit("arg error - need 1 arg: [RI file folder]")
