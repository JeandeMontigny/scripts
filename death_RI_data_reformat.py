#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(output_file):

    param=0
    fichier=open(output_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:

        if line[0] != 'a':
            if line[0] == 'D':
                param = float(line[2:])
            continue

        m=re.search( r'average ri: (.+) with std: (.+) ; cell death: (.+)', line, re.M|re.I)
        if m:
            ri = m.group(1)
            death = m.group(3)
            print(param, ri, ri, ri, 1700, death, 0)
#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [death_RI file]")

