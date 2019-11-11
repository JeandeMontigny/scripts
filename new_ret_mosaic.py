#!/usr/bin/env python3
import sys, os, re
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(folder, process_position = True, process_ri = True):
    m=re.search( r'.+/results(.+)/', folder, re.M|re.I)
    folder_nb = m.group(1)
    if process_ri:
        processRi(folder+"RI_"+folder_nb+".txt")
    if process_position:
        processPosition(folder+"cells_position/")

#--------------------------------------------------------------------------#
def processRi(file):
    type_list = []; all_steps = []
    current_step = 0; max_step = 0
    type_ri_density = []

    fichier=open(file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # step ri cell_type cell_population total_death
        m = re.search( r'(.+) (.+) (.+) (.+) (.+)', line, re.M|re.I)
        if m:
            step = int(m.group(1)); ri = float(m.group(2)); type = int(m.group(3))
            cell_population = int(m.group(4)); total_death = float(m.group(5))
            if type == -1:
                continue
            if not type in type_list:
                type_list.append(type)
                all_steps.append([type, ri])
                type_ri_density.append([type, [ri, cell_population]])
            else:
                for i in range(0, len(all_steps)):
                    if len(all_steps[i]) > max_step:
                        max_step = len(all_steps[i])
                    if type == all_steps[i][0]:
                        all_steps[i].append(ri)
                    if type == type_ri_density[i][0]:
                        type_ri_density[i].append([ri, cell_population])

    for i in range(0, len(all_steps)):
        if len(all_steps[i]) < max_step:
            for j in range(0, max_step - len(all_steps[i])):
                all_steps[i].insert(1,np.nan)

    riFigure(all_steps, all = True)

    processLastStep(type_ri_density)

    plt.show()

#--------------------------------------------------------------------------#
def riFigure(all_types, all = False):
    plt.figure()
    if all:
        for i in range(0, len(all_types)):
            if all_types[i][0] == 204\
                or all_types[i][0] == 108\
                or all_types[i][0] == 115:
                plt.plot(all_types[i][1:],\
                    label="cell type: " + str(all_types[i][0]))

    ave_ri = []; temps = []
    for i in range(0, len(all_types[0])-1):
        for j in range(0, len(all_types)):
            temps.append(all_types[j][i])
        ave_ri.append(np.average(temps))
        temps = []

    if all:
        plt.plot(ave_ri, linewidth=8, color = 'black')
        plt.ylim(1, 6.5)
        plt.legend()
    else:
        plt.plot(ave_ri, color = 'black')
        plt.ylim(1.5, 3.5)
    plt.xlim(0, len(all_types[0])-1)
    plt.title("Average regularity index")

#--------------------------------------------------------------------------#
def processLastStep(type_ri_density):
    ri = []; density = []
    for i in range(0, len(type_ri_density)):
        ri.append(type_ri_density[i][len(type_ri_density[i])-1][0])
        density.append(type_ri_density[i][len(type_ri_density[i])-1][1])

        last_step.append([ type_ri_density[i][0],
            type_ri_density[i][len(type_ri_density[i])-1] ])

    print("Final step: ri-density correlation magnitude:",\
        round(np.corrcoef(ri, density)[0][1], 2))

    plt.figure()
    plt.plot(density, ri, 'o', color = "black")
    plt.plot(np.unique(density), np.poly1d(np.polyfit(density, ri, 1))(np.unique(density)), color = "red")
    plt.xlabel("Cell density")
    plt.ylabel("RI")

#--------------------------------------------------------------------------#
def processPosition(folder):
    return

#--------------------------------------------------------------------------#
# check arguments
if len(sys.argv) == 2 or len(sys.argv) == 3 or len(sys.argv) == 4:

    if len(sys.argv) and sys.argv[1].endswith(".txt"):
        processRi(sys.argv[1])
        sys.exit()

    if len(sys.argv) == 4:
        if sys.argv[2] == "-p" or sys.argv[2] == "-ri"\
            and sys.argv[3] == "-p" or sys.argv[3] == "-ri":
            main(sys.argv[1], True, True)
            sys.exit()

    if len(sys.argv) == 3:
        if sys.argv[2] == "-p":
            main(sys.argv[1], True, False)
            sys.exit()
        if sys.argv[2] == "-ri":
            main(sys.argv[1], False, True)
            sys.exit()
        else:
            raise SystemExit('Error: invalid option. options: [-p, -ri]')

    else:
        raise SystemExit('Error: invalid option. options: [-p, -ri]')
else:
    raise SystemExit('Error: need at least 1 arg: [result folder]; option: [-p, -ri]')

print("done")
