#!/usr/bin/env python3
import sys, os, re
import numpy as np

#--------------------------------------------------------------------------#
def main(folder, output_path):
    # create and open output file
    output=open(output_path+"ri_measures_output.out", "a")
    # walk across folders
    for dirpath, dirnames, filenames in os.walk(folder):
        # for each files in directory
        for name in filenames:
            # if file is a .txt
            if name.endswith(".txt"):
                # call read function and get ri value for this file. round value to only 2 decimals
                ri = round(read(dirpath+"/"+name), 2)
                # create string for output
                output_string=dirpath[35:]+" "+name+" "+str(ri)
                print(output_string)
                # write output inside the output file
                output.write(output_string+"\n")
    # close output file
    output.close()

#--------------------------------------------------------------------------#
def read(coord_file):
    """Return RI value. To do that, extract x,y coordinate from ImageJ cell position extraction then compute the Regularity Index (RI) corresponding to this cell population.
    Argument: file containing X,Y coordinates of every cells, sorted as 'index Area Mean Min Max X Y' """
    cells_position=[]
    fichier=open(coord_file, "r")
    file_lines=fichier.readlines()
    # first, extract X Y position for each cell
    # for each line of this file
    for line in file_lines:
        # if line doesn't start with space character, ie, if it's not the file legend
        if line[0]!= " ":
            # file strucure: index Area Mean Min Max X Y
            m=re.search( r'[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t([0-9.]+)\t([0-9.]+)', line, re.M|re.I)
            if m :
                cells_position.append([round(float(m.group(1)), 2), round(float(m.group(2)), 2)])

    # then call getShortestDistList function for this whole file
    shortest_dist_list=getShortestDistList(cells_position)
    # calculate ri, which is: mean of distances / std of distances
    ri = np.average(shortest_dist_list)/np.std(shortest_dist_list)

    return ri

#--------------------------------------------------------------------------#
def getShortestDistList(coord_list):
    """Return a list containing shortest neighbour distance for each cell
    Argument: tab of cells coordinates"""
    shortest_dist_list = []
    # for each cell coordinate in this list
    for i in range(0, len(coord_list)):
        distance_list = []
        cell_coord = coord_list[i]
        # for each other cell in this list
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            # calculate Euclidean distance between these two cells
            tempsDistance = np.sqrt(np.square(cell_coord[0] - other_cell_coord[0]) + np.square(cell_coord[1] - other_cell_coord[1]))
            # if cell is not itself
            if tempsDistance != 0:
                distance_list.append(tempsDistance)
        # add shortest distance
        shortest_dist_list.append(min(distance_list))

    return shortest_dist_list

#--------------------------------------------------------------------------#
if len(sys.argv)!=3:
    exit("arg error - need 2 arg: [folder with position files] [absolute path of output file]")
else:
    main(sys.argv[1], sys.argv[2])
    print("done")
