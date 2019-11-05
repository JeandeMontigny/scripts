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

    fichier=open(file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # step ri cell_type total_death
        m=re.search( r'(.+) (.+) (.+) (.+)', line, re.M|re.I)
        step = int(m.group(1)); ri = float(m.group(2)); type = int(m.group(3))
        if m:
            if type == -1:
                continue
            if not type in type_list:
                type_list.append(type)
                all_steps.append([type, ri])
            else:
                for i in range(0, len(all_steps)):
                    if len(all_steps[i]) > max_step:
                        max_step = len(all_steps[i])
                    if type == all_steps[i][0]:
                        all_steps[i].append(ri)

    for i in range(0, len(all_steps)):
        if len(all_steps[i]) < max_step:
            for j in range(0, max_step - len(all_steps[i])):
                all_steps[i].insert(1,np.nan)

    ri_figure(all_steps, all = True)

#--------------------------------------------------------------------------#
def ri_figure(all_types, all = False):
    plt.figure()
    if all:
        for i in range(0, len(all_types)):
            plt.plot(all_types[i][1:], label=all_types[i][0])

    ave_ri = []; temps = []
    for i in range(0, len(all_types[0])-1):
        for j in range(0, len(all_types)):
            temps.append(all_types[j][i])
        ave_ri.append(np.average(temps))
        temps = []

    if all:
        plt.plot(ave_ri, linewidth=8, color = 'black')
        plt.ylim(1, 6)
    else:
        plt.plot(ave_ri, color = 'black')
        plt.ylim(1.5, 3.5)
    plt.xlim(0, len(all_types[0])-1)
    plt.title("Average regularity index")

    plt.show()

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
