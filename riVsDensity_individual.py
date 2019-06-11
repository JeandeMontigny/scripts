#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(fodler):
    list_density=[]; list_ri=[]
    for file in os.listdir(fodler):
        if file.endswith(".txt"):
            output = read_file(fodler, file)
            list_density.append(output[0]); list_ri.append(output[1])

    figure_construction(list_density, list_ri)
    plt.show()

    return 1

#--------------------------------------------------------------------------#
def read_file(fodler, file):
    density=0; line_count=0
    list_ri=[]; list_ri_temps=[]
    file_lines=open(fodler+file, "r").readlines()
    for line in file_lines:
        m=re.search( r'(.+) (.+)', line, re.M|re.I)
        if m:
            line_count+=1
            density = int(m.group(1))
            if np.mod(line_count, 200)==0 and len(list_ri_temps)!=0:
                list_ri_temps.append(float(m.group(2)))
                list_ri.append(list_ri_temps)
                list_ri_temps=[]
            else:
                list_ri_temps.append(float(m.group(2)))

    return (density, merge_lists(list_ri))

#--------------------------------------------------------------------------#
def merge_lists(list_ri):
    list_ave_ri=[]
    for i in range(0, len(list_ri[0])):
        temps=0
        for j in range(0, len(list_ri)):
            temps+=list_ri[j][i]
        list_ave_ri.append(temps/len(list_ri))

    return list_ave_ri

#--------------------------------------------------------------------------#
def figure_construction(list_density, list_ri):
    plt.figure(1)
    for i in range(0, len(list_ri)):
        plt.plot(list(range(0, len(list_ri[i]))), list_ri[i], label=list_density[i])

    plt.ylim(1.5, 11)
    plt.xlabel("simulation time")
    plt.ylabel("regularity index")
#    plt.legend(fontsize="medium", ncol=2, bbox_to_anchor=(1.04,0.5), loc="center left", title="Density")
    plt.legend(loc=0)
    plt.title("layer colapse impact on RI")

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    exit("arg error - need 1 arg: [individual_RI_density_*.txt folder]")
