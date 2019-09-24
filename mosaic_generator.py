#!/usr/bin/env python3
import sys, os, re
import numpy as np
import random as rd

#--------------------------------------------------------------------------#
def main(output_folder, weight):
    if float(weight) > 1 or float(weight) < 0:
        SystemExit('random weigh should be between 0 and 1')
    cell_per_dim = 20
    cells_space = 25
    rand = float(weight) * cells_space

    if os.path.isfile(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+weight+".txt"):
        os.remove(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+weight+".txt")
    output_file=open(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+weight+".txt", "a")

    for i in range(0, cell_per_dim):
        for j in range(0, cell_per_dim):
            position = str( round(rd.uniform(-rand, rand), 2)+(i*cells_space) ) + " " + str( round(rd.uniform(-rand, rand), 2)+(j*cells_space) ) + "\n"
            output_file.write(position)
    output_file.close()

    return 1

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==3:
    if main(sys.argv[1], sys.argv[2]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [mosaic ouput folder] [random weigh: 0-1]')
