#!/usr/bin/env python3
import sys, os, re
import numpy as np

# compute RI from imagej cell position extraction
#--------------------------------------------------------------------------#
def main(coord_file):
    """Extract x,y coordinate from ImageJ cell position extraction then compute the Regularity Index (RI) corresponding to this cell population.
    Argument: file containing X,Y coordinates of every cells, sorted as "index,x,y" format """
    cells_position=[]
    fichier=open(coord_file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # index,x,y
        # m=re.search( r'[0-9]+,([0-9]+),([0-9]+)', line, re.M|re.I)
        # with "point selection" method: index,Area,Mean,Min,Max,X,Y
        m=re.search( r'[0-9]+,[0-9.]+,[0-9.]+,[0-9.]+,[0-9.]+,([0-9.]+),([0-9.]+)', line, re.M|re.I)
        if m :
            cells_position.append([float(m.group(1)), float(m.group(2))])

    shortest_dist_list=getShortestDistList(cells_position)
    ri = np.average(shortest_dist_list)/np.std(shortest_dist_list)
    print("ri for this mosaic:", round(ri, 2))

#--------------------------------------------------------------------------#
def getShortestDistList(coord_list):
    shortest_dist_list = []
    for i in range(0, len(coord_list)):
        distance_list = []
        cell_coord = coord_list[i]
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            tempsDistance = np.sqrt(np.square(cell_coord[0] - other_cell_coord[0]) + np.square(cell_coord[1] - other_cell_coord[1]))
            # if cell is not itself
            if tempsDistance != 0:
                distance_list.append(tempsDistance)
        # add shortest distance
        shortest_dist_list.append(min(distance_list))
    return shortest_dist_list

#--------------------------------------------------------------------------#
if len(sys.argv)!=2:
    exit("arg error - need 1 arg: [cells_position_file]")
else:
    main(sys.argv[1])
    print("done")

