#!/usr/bin/env python3
import sys, os, re
import numpy as np
import random as rd
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(output_folder):
    # only one monolayer mosaic is considered here
    nb_cells = 400 # 200 per nb_types
    space = 500
    nb_types = 2
    #TODO: test optimal death rate
    death_rate = 0.65

    density = nb_cells/nb_types * 1e6 / (space * space)

    cells = createCells(nb_cells, space, nb_types, dif=True)
    initial_ri = getRi(cells, nb_types)

    # plotMosaic(cells, nb_types)

    print("initial density (per cell type):", int(density), "- initial ri:", initial_ri)

#    fate_mosaic = fate(cells)
#    death_mosaic = death(fate_mosaic, death_rate)
    death_mosaic = death(cells, death_rate)
    migration_mosaic = migration(cells, repetition=25, exclusion=30)

#    fate_ri = getRi(fate_mosaic, nb_types)
    death_ri = getRi(death_mosaic, nb_types)
    migration_ri = getRi(migration_mosaic, nb_types)

    print()
#    print("non AB fate:", fate_ri)
    print("non AB death:", death_ri)
    print("non AB migration:", migration_ri)

#    plotMosaic(fate_mosaic, nb_types, "Mosaic created with non AB fate mechanism")
    plotMosaic(death_mosaic, nb_types, "Mosaic created with non AB death mechanism")
    plotMosaic(migration_mosaic, nb_types, "Mosaic created with non AB migration mechanism")

    plt.show()

    return 1
#--------------------------------------------------------------------------#
def plotMosaic(cells, nb_types, title):
    plt.figure()
    for i in range(0, len(cells)):
        plt.plot(cells[i][0], cells[i][1], 'o', color=rgb(0, nb_types, cells[i][2]))
    plt.title(title)

#--------------------------------------------------------------------------#
def rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r

    return r/255, g/255, b/255

#--------------------------------------------------------------------------#
def createCells(nb_cells, space, nb_types, dif=True):
    cells = []
    for i in range(0, nb_cells):
        new_cell = [round(rd.uniform(0, space), 2), round(rd.uniform(0, space), 2), int(rd.uniform(0, nb_types))] if dif else [round(rd.uniform(0, space), 2), round(rd.uniform(0, space), 2), -1]
        while positionConflict(cells, new_cell):
            new_cell = [round(rd.uniform(0, space), 2), round(rd.uniform(0, space), 2), int(rd.uniform(0, nb_types))] if dif else [round(rd.uniform(0, space), 2), round(rd.uniform(0, space), 2), -1]
        cells.append(new_cell)

    return cells

#--------------------------------------------------------------------------#
def positionConflict(cells, new_cell):
    for cell in cells:
        if twodDistance(cell, new_cell) < 7:
            return True
    return False

#--------------------------------------------------------------------------#
def writeRi(positions_list, nb_types):
    #TODO
    getRi(positions_list, nb_types)

#--------------------------------------------------------------------------#
def getRi(positions_list, nb_types):
    ri = []
    for type in range(0, nb_types):
        shortest_dist_list = getDistLists(positions_list, type)
        ri.append(round((np.average(shortest_dist_list)/np.std(shortest_dist_list)), 2))

    return round(np.average(ri), 2), round(np.std(ri), 2)

#--------------------------------------------------------------------------#
def getDistLists(coord_list, type):
    shortest_dist_list = []
    for i in range(0, len(coord_list)):
        distance_list = []
        cell_coord = coord_list[i]
        cell_type = coord_list[i][2]
        if cell_type != type:
            continue
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            other_cell_type = coord_list[j][2]
            if other_cell_type != type:
                continue
            tempsDistance = twodDistance(cell_coord, other_cell_coord)
            # if cell is not itself
            if tempsDistance != 0:
                distance_list.append(tempsDistance)
        # add shortest distance
        shortest_dist_list.append(min(distance_list))

    return shortest_dist_list

#--------------------------------------------------------------------------#
def twodDistance(a, b):
    return np.sqrt(np.square(a[0] - b[0]) + np.square(a[1] - b[1]))

#--------------------------------------------------------------------------#
def newCoord(cell_ref, cell_to_move, dist, factor=4):
    x = cell_to_move[0] + (factor*(cell_to_move[0]-cell_ref[0])/dist)
    y = cell_to_move[1] + (factor*(cell_to_move[1]-cell_ref[1])/dist)
    return round(x, 2), round(y, 2), cell_to_move[2]

#--------------------------------------------------------------------------#
def fate(cells, write=False):
    for i in range(0, len(cells)):
        cell_a = cells[i]
        if cell_a[2] == -1:
            closest_cell = 1e6; closest_cell_type = -1
            for j in range(0, len(cells)):
                cell_b = cells[j]
                if cell_b[2] != -1:
                    distance = twodDistance(cell_a, cell_b)
                    if distance < closest_cell:
                        closest_cell = distance
                        closest_cell_type = cell_b[2]
            # if no differentiated cells yet
            if closest_cell_type == -1:
                new_type = int(rd.uniform(0, 2))
            elif closest_cell_type == 1:
                new_type = 0
            elif closest_cell_type == 0:
                new_type =1
            cells[i] = [cell_a[0], cell_a[1], new_type]
            if write:
                writeRi(cells, 2)

    return cells

#--------------------------------------------------------------------------#
def death(cells, death_rate):
    repetition = int(len(cells)*death_rate)
    for iteration in range(0, repetition):
        print("death mechanism iteration", iteration, "out of", repetition, end="\r", flush=True)
        # find closest cells couple
        closest_cells = 1e6; closest_couple = []
        for i in range(0, len(cells)):
            cell_a = cells[i]
            for j in range(0, len(cells)):
                cell_b = cells[j]
                # if homotypic cells
                if cell_b[2] == cell_a[2]:
                    distance = twodDistance(cell_a, cell_b)
                    if distance > 0 and distance < closest_cells:
                        closest_cells = distance
                        closest_couple = [cell_a, cell_b]
        # decide if cell_a or cell_b has to be removed
        second_closest_a = 1e6
        cell_a = closest_couple[0]
        for i in range(0, len(cells)):
            cell_c = cells[i]
            distance = twodDistance(cell_a, cell_c)
            if distance < second_closest_a and distance > closest_cells:
                second_closest_a = distance

        second_closest_b = 1e6
        cell_b = closest_couple[1]
        for i in range(0, len(cells)):
            cell_c = cells[i]
            distance = twodDistance(cell_b, cell_c)
            if distance < second_closest_b and distance > closest_cells:
                second_closest_b = distance
        # select cell to remove between a and b
        cell_to_remove = cell_a if second_closest_a < second_closest_b else cell_b
        # remove selected cell
        cells.remove(cell_to_remove)

    return cells

#--------------------------------------------------------------------------#
def migration(cells, repetition=25, exclusion=25):
    for iteration in range(0, repetition):
        print("migration mechanism iteration", iteration, "out of", repetition, end="\r", flush=True)
        # for each cell, repulse all cells that are too close
        for i in range(0, len(cells)):
            cell_a = cells[i]
            for j in range(0, len(cells)):
                cell_b = cells[j]
                # if homotypic cells
                if cell_b[2] == cell_a[2]:
                    distance = twodDistance(cell_a, cell_b)
                    if distance > 0 and distance < exclusion:
                        cells[j] = newCoord(cell_a, cell_b, distance)

    return cells

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [output folder]')
