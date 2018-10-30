#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(distance_file):
    distance=[]
    fichier=open(distance_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        m=re.search( r'(.+)', line, re.M|re.I)
        if m:
            distance.append(float(m.group(1)))

    print("average tangential distance:", round(np.average(distance),2), "micrometers ; with std of", round(np.std(distance),2))
    figure(distance)

#--------------------------------------------------------------------------#
def figure(distance):
    plt.figure(1)
    plt.hist(distance)
    plt.xlabel("distance")
    plt.ylabel("number of cell")
    plt.title("distribution of tangential migration distance")

    plt.show()

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [migration distance file]")
