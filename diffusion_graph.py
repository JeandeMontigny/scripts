#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(diffusion_file):
    x=[]; y=[];
    fichier=open(diffusion_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        m=re.search( r'(.+) (.+)', line, re.M|re.I)
        if m:
            x.append(float(m.group(1)))
            y.append(float(m.group(2)))

    figure(x, y)

#--------------------------------------------------------------------------#
def figure(x, y):
    plt.figure(1)

    # draw cell boundaries
    plt.axvline(x=45, color='silver')
    plt.axvline(x=57, color='silver')

    plt.plot(x,y)
    plt.xlabel("distance")
    plt.ylabel("concentration")

    plt.show()

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [diffusion file]")
