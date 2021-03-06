#!/usr/bin/env python3
import sys, os, re
import numpy as np

#--------------------------------------------------------------------------#
def main(gcl_coord_file, inl_coord_file):
    """Extract x,y coordinate from ImageJ cell position extraction of amacrine cells in GCl and INL layer.
    Compute the mosaic Regularity Index (RI) corresponding to each layer as well as the combined mosaic.
    Compute the exclusion factor that represents how mosaics exclude one each other (ie if cells dont overlap)
    Argument: files containing X,Y coordinates of every cells"""

    # extract cells position
    gcl_cells_position = getCellsPosition(gcl_coord_file)
    inl_cells_position = getCellsPosition(inl_coord_file)

    # independant mosaics
    gcl_shortest_dist_list = getShortestDistList(gcl_cells_position)
    inl_shortest_dist_list = getShortestDistList(inl_cells_position)
    gcl_ri = getRi(gcl_shortest_dist_list)
    inl_ri = getRi(inl_shortest_dist_list)
    print("gcl ri:", gcl_ri)
    print("inl ri:", inl_ri)
    print("GCL/INL ratio:", round(len(gcl_cells_position)/len(inl_cells_position), 2))

    # combined mosaic
    both_cells_position = mergeTab(gcl_cells_position, inl_cells_position)
    both_shortest_dist_list = getShortestDistList(both_cells_position)
    both_ri = getRi(both_shortest_dist_list)
    print("combined ri:", both_ri)

    # exclusion
    exclusionFactor = getExclusionFactor(gcl_cells_position, inl_cells_position)
    print("exclusion factor:", exclusionFactor)

#--------------------------------------------------------------------------#
def getExclusionFactor(tab_1, tab_2):
    """ return a double between [0,1]: 1 if mosaics completely excluse one each other, 0 if they are totaly superimposed """
    redondancy = 0
    r = 10
    for cell_1 in tab_1:
        for cell_2 in tab_2:
            if (np.square(cell_1[0] - cell_2[0]) + np.square(cell_1[1] - cell_2[1]) <= np.square(r)):
                redondancy+=1
                # exit loop to avoid multiple positif for one cell
                break
    return round(1 - redondancy/len(tab_1+tab_2)*2, 2)

#--------------------------------------------------------------------------#
def mergeTab(tab_1, tab_2):
    return (tab_1 + tab_2)

#--------------------------------------------------------------------------#
def getRi(shortest_dist_list):
    return round((np.average(shortest_dist_list)/np.std(shortest_dist_list)), 2)

#--------------------------------------------------------------------------#
def getCellsPosition(coord_file):
    cells_position=[]
    fichier=open(coord_file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # ...,x,y
        m=re.search( r'[0-9]+,[0-9.]+,[0-9.]+,[0-9.]+,[0-9.]+,([0-9.]+),([0-9.]+)', line, re.M|re.I)
        if m :
            cells_position.append([float(m.group(1)), float(m.group(2))])
    return cells_position

#--------------------------------------------------------------------------#
def getShortestDistList(coord_list):
    shortest_dist_list = []
    for i in range(0, len(coord_list)):
        shortestDist = 100000
        cell_coord = coord_list[i]
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            tempsDistance = np.sqrt(np.square(cell_coord[0] - other_cell_coord[0]) + np.square(cell_coord[1] - other_cell_coord[1]))
            # if cell is not itself and is closer
            if tempsDistance != 0 and tempsDistance < shortestDist:
                shortestDist = tempsDistance
        shortest_dist_list.append(shortestDist)
    return shortest_dist_list

#--------------------------------------------------------------------------#
if len(sys.argv)!=3:
    exit("arg error - need 2 arg: [GCL_cells_coord], [INL_cells_coord]")
else:
    main(sys.argv[1], sys.argv[2])
    print("done")

