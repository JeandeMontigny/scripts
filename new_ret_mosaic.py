#!/usr/bin/env python3
import sys, os, re
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(folder, process_multiple = True, process_ri = True):
    m=re.search( r'.+/results(.+)/', folder, re.M|re.I)
    folder_nb = m.group(1)
    if process_ri:
        processRi(folder+"RI_"+folder_nb+".txt")
    if process_position:
        processPosition(folder+"cells_position/")

#--------------------------------------------------------------------------#
def processSingle(file, multiple = False):
    type_list = []; type_ri_density = []
    current_step = 0; max_step = 0

    fichier=open(file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # step ri cell_type cell_population total_death
        m = re.search( r'(.+) (.+) (.+) (.+) (.+)', line, re.M|re.I)
        if m:
            step = int(m.group(1)); ri = float(m.group(2))
            type = int(m.group(3)); cell_population = int(m.group(4))
            total_death = float(m.group(5))
            if type == -1:
                continue
            if not type in type_list:
                type_list.append(type)
                type_ri_density.append([type, [ri, cell_population]])
            else:
                for i in range(0, len(type_ri_density)):
                    if len(type_ri_density[i]) > max_step:
                        max_step = len(type_ri_density[i])
                    if type == type_ri_density[i][0]:
                        type_ri_density[i].append([ri, cell_population])

    for i in range(0, len(type_ri_density)):
        if len(type_ri_density[i]) < max_step:
            for j in range(0, max_step - len(type_ri_density[i])):
                type_ri_density[i].insert(1,[np.nan, np.nan])

    if multiple:
        return type_ri_density

    riFigure(type_ri_density, all = True)
    processLastStep(type_ri_density)

#--------------------------------------------------------------------------#
def riFigure(type_ri_density, all = False):
    plt.figure()
    if all:
        for i in range(0, len(type_ri_density)):
            cell_type = type_ri_density[i][0]
            # if cell_type == 204 or cell_type == 108 or cell_type == 115:
            plt.plot([couple[0] for couple in type_ri_density[i][1:]],\
                label="Final cell density: "\
                + str(type_ri_density[i][len(type_ri_density[i])-1][1]))

    ave_ri = []; temps = []
    # for each step
    for i in range(1, len(type_ri_density[0])-1):
        # for each type
        for j in range(0, len(type_ri_density)):
            temps.append(type_ri_density[j][i][0])
        ave_ri.append(np.average(temps))
        temps = []

    if all:
        plt.plot(ave_ri, linewidth=4, color = 'black', label = "average ri")
        plt.ylim(1, 6.5)
        plt.legend()
    else:
        plt.plot(ave_ri, color = 'black')
        plt.ylim(1.5, 3.5)
    plt.xlim(0)
    plt.axvline(12, linestyle='--', color = "gray")
    plt.axvline(67, linestyle='--', color = "gray")
    plt.xlabel("Simulation time")
    plt.ylabel("Regularity index")

#--------------------------------------------------------------------------#
def processLastStep(type_ri_density):
    ri = []; density = []
    for i in range(0, len(type_ri_density)):
        ri.append(type_ri_density[i][len(type_ri_density[i])-1][0])
        density.append(type_ri_density[i][len(type_ri_density[i])-1][1])

    print("Final step: ri-density correlation magnitude:",\
        round(np.corrcoef(ri, density)[0][1], 2))

    plt.figure()
    plt.plot(density, ri, 'o', color = "black")
    plt.plot(np.unique(density), np.poly1d(np.polyfit(density, ri, 1))(np.unique(density)), color = "red")
    plt.xlabel("Cell density")
    plt.ylabel("Regularity index")

#--------------------------------------------------------------------------#
def processMultiple(output_folder):
    sub_folders = [name for name in os.listdir(output_folder) if os.path.isdir(os.path.join(output_folder, name))]
    simus = []
    for folder in sub_folders:
        m=re.search( r'results(.+)', folder, re.M|re.I)
        if m:
            folder_nb = m.group(1)
            simus.append(processSingle(output_folder+"/"+folder+"/RI_"+folder_nb+".txt", multiple = True))

    # [ [type, [ri11, ri12], [ri21, ri22] ]
    #   [type, [ri11, ri12], [ri21, ri22] ] ]
    all_types_all_ri = []; all_types_all_pop = []
    list_type = []
    for simu in simus:
        # for each cell type in this simu
        for i in range(0, len(simu)):
            type = simu[i][0]
            if not type in list_type:
                list_type.append(type)
                all_types_all_ri.append([type])
                all_types_all_pop.append([type])
                # for each step in this cell type of this simu
                for j in range(1, len(simu[i])-1):
                    ri = simu[i][j][0]
                    pop = simu[i][j][1]
                    all_types_all_ri[len(all_types_all_ri)-1].append([ri])
                    all_types_all_pop[len(all_types_all_ri)-1].append([pop])

            else:
                for z in range(0, len(all_types_all_ri)):
                    # get correct type in list
                    if type == all_types_all_ri[z][0]:
                        for j in range(1, len(simu[i])-1):
                            ri = simu[i][j][0]
                            pop = simu[i][j][1]
                            all_types_all_ri[z][j].append(ri)
                            all_types_all_pop[z][j].append(pop)

    ave_ri = []; std_ri = []
    ave_pop = []; std_pop = []
    # for each cell type
    for i in range (0, len(all_types_all_ri)):
        # for each step
        ave_ri.append([all_types_all_ri[i][0]])
        std_ri.append([all_types_all_ri[i][0]])
        ave_pop.append([all_types_all_pop[i][0]])
        std_pop.append([all_types_all_pop[i][0]])
        for j in range(1, len(all_types_all_ri[i])):
            ave_ri[i].append(np.average(all_types_all_ri[i][j]))
            std_ri[i].append(np.std(all_types_all_ri[i][j]))
            ave_pop[i].append(np.average(all_types_all_pop[i][j]))
            std_pop[i].append(np.std(all_types_all_pop[i][j]))

    processLastStepMultiple(all_types_all_ri, all_types_all_pop)

    plt.figure()
    plt.errorbar([x for x in range(0, len(ave_ri[35])-1)],\
                 ave_ri[35][1:], yerr=std_ri[35][1:],\
                 label = str(ave_ri[35][0]))
    plt.legend()
    plt.show()

    return

#--------------------------------------------------------------------------#
def processLastStepMultiple(all_types_all_ri, all_types_all_pop):
    # ri vs density plot with x and y error bar

    return

#--------------------------------------------------------------------------#
# check arguments
if len(sys.argv) == 2:

    if len(sys.argv) and sys.argv[1].endswith(".txt"):
        processSingle(sys.argv[1])
    else:
        processMultiple(sys.argv[1])

else:
    raise SystemExit('Error: need at least 1 arg: [output folder, RI file]')

plt.show()
print("done")
