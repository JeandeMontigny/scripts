#!/usr/bin/env python3
import sys, os, re, numpy
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(output_file):
    tab_on=[]; tab_off=[]; tab_onoff=[]; tab_ave=[]
    fichier=open(output_file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # on off on-off
        m=re.search( r'(.+) (.+) (.+) .+', line, re.M|re.I)
        if m:
            tab_on.append(float(m.group(1)))
            tab_off.append(float(m.group(2)))
            tab_onoff.append(float(m.group(3)))
            tab_ave.append((float(m.group(1))+float(m.group(2))+float(m.group(3)))/3)

    figure(tab_on, tab_off, tab_onoff, tab_ave)

#--------------------------------------------------------------------------#
def figure(tab1, tab2, tab3, tab4):
    plt.figure(1)
    plt.plot(tab1, label="On cells")
    plt.plot(tab2, label="Off cells")
    plt.plot(tab3, label="On-Off cells")
    plt.plot(tab4, label="Mean")
    plt.xlabel("modelling time")
    plt.ylabel("regularity index")
    plt.axis([0, len(tab1), 0, 10]) # 15, 1, 5.5
    plt.legend()

    plt.show()

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [RI file]")

