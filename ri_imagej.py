#!/usr/bin/env python3
import sys, os, re
import numpy as np
#import matplotlib.pyplot as plt

# compute RI from imagej cell position extraction
#--------------------------------------------------------------------------#
def main(coord_file):
    cells_position=[]
    fichier=open(coord_file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # index,x,y
#        m=re.search( r'[0-9]+,([0-9]+),([0-9]+)', line, re.M|re.I)
        # with "point selection": index,Area,Mean,Min,Max,X,Y
        m=re.search( r'[0-9]+,[0-9.]+,[0-9.]+,[0-9.]+,[0-9.]+,([0-9.]+),([0-9.]+)', line, re.M|re.I)
        if m :
            cells_position.append([float(m.group(1)), float(m.group(2))])

    shortest_dist_list=getShortestDistList(cells_position)
    ri = np.average(shortest_dist_list)/np.std(shortest_dist_list)
    print("ri for that mosaic:", round(ri, 2))


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
if len(sys.argv)!=2:
    exit("arg error - need 2 arg: [cells_position_file]")
else:
    main(sys.argv[1])
    print("done")

