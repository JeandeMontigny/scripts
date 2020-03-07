#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(distance_file):
    cells_gcl = []; cells_inl = []
    cells_gcl_x = []; cells_gcl_y = []
    cells_inl_x = []; cells_inl_y = []
    fichier = open(distance_file, "r")
    file_lines = fichier.readlines()
    for  line in file_lines:
        m = re.search( r'(.+) (.+) (.+) (.+)', line, re.M|re.I)
        if m:
            if int(m.group(1)) == 0:
                cells_gcl.append([float(m.group(2)), float(m.group(3))])
                cells_gcl_x.append(float(m.group(2)))
                cells_gcl_y.append(float(m.group(3)))
            elif int(m.group(1)) == 1:
                cells_inl.append([float(m.group(2)), float(m.group(3))])
                cells_inl_x.append(float(m.group(2)))
                cells_inl_y.append(float(m.group(3)))

    ri_gcl = get_ri(cells_gcl)
    ri_inl = get_ri(cells_inl)
    ri_combined = get_ri(cells_gcl + cells_inl)
    print("ri gcl:", ri_gcl, "; ri inl:", ri_inl, "; ri combined:", ri_combined)
    exclusion_factor = get_exclusion_fact(cells_gcl, cells_inl)
    print("gcl-inl exclusion factor:", exclusion_factor, "(1 = complete exclusion)")

    figure(cells_gcl_x, cells_gcl_y, cells_inl_x, cells_inl_y)

#--------------------------------------------------------------------------#
def get_ri(cells):
    shortest_dist_list = get_shortest_dist_list(cells)

    return round((np.average(shortest_dist_list) / np.std(shortest_dist_list)) , 2)

#--------------------------------------------------------------------------#
def get_shortest_dist_list(cells):
    shortest_dist_list = []
    for i in range(0, len(cells)):
        distance_list = []
        cell_a = cells[i]
        for j in range(0, len(cells)):
            cell_b = cells[j]
            temps_distance = get_distance(cell_a, cell_b)
            if temps_distance != 0:
                distance_list.append(temps_distance)
        shortest_dist_list.append(min(distance_list))

    return shortest_dist_list

#--------------------------------------------------------------------------#
def get_distance(cell_a, cell_b):

    return round(np.sqrt( \
        np.square(cell_a[0] - cell_b[0]) + \
        np.square(cell_a[1] - cell_b[1]) ), 2)

#--------------------------------------------------------------------------#
def get_exclusion_fact(pop_a, pop_b):
    redondancy = 0
    r = 16
    for cell_a in pop_a:
        for cell_b in pop_b:
            if (np.square(cell_a[0] - cell_b[0]) + np.square(cell_a[1] - cell_b[1]) <= np.square(r)):
                redondancy+=1
                # exit loop to avoid multiple positif for one cell
                break
    return round(1 - redondancy/len(pop_a+pop_b)*2, 2)

#--------------------------------------------------------------------------#
def figure(cells_gcl_x, cells_gcl_y, cells_inl_x, cells_inl_y):
    plt.figure()
    sub_figure(cells_gcl_x, cells_gcl_y, c='b')
    sub_figure(cells_inl_x, cells_inl_y, c='r')
    # plt.legend(["gcl", "inl"])
    # sub_figure(cells_gcl_x+cells_inl_x, cells_gcl_y+cells_inl_y)

    plt.show()

#--------------------------------------------------------------------------#
def sub_figure(cells_x, cells_y, c='b'):
    plt.scatter(cells_x, cells_y, color=c, s=4)
    plt.xlabel("x position")
    plt.ylabel("y position")

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [SAC position export file]")
