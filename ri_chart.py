#!/usr/bin/env python3
import sys, os, re, numpy
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(output_file):
    tab_on=[]; tab_off=[]; tab_onoff=[]; tab_ave=[]
    fichier=open(output_file, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        m=re.search( r'(.+) (.+) (.+)', line, re.M|re.I)
#        m=re.search( r'(.+) (.+)', line, re.M|re.I) # read RI for on (1) and off (2) cells
        if m:
            tab_on.append(float(m.group(1)))
            tab_off.append(float(m.group(2)))
            tab_onoff.append(float(m.group(3)))
            tab_ave.append((float(m.group(1))+float(m.group(2))+float(m.group(3)))/3)
#            tab_ave.append((float(m.group(1))+float(m.group(2)))/2)

#    for index in range(10):
#        tab_on[index]=None
#        tab_off[index]=None
#        tab_ave[index]=None

    figure(tab_on, tab_off, tab_onoff, tab_ave)
#    figure(tab_on, tab_off, tab_ave)

#--------------------------------------------------------------------------#
def figure(tab1, tab2, tab3, tab4):
#def figure(tab1, tab2, tab4):
#    fig_ave=plt.figure(1)
    plt.figure(1)
    plt.plot(tab1, label="On cells")
    plt.plot(tab2, label="Off cells")
    plt.plot(tab3, label="On-Off cells")
    plt.plot(tab4, label="All cells")
    plt.xlabel("modelling time")
    plt.ylabel("regularity index")
#    plt.axis([0, len(tab1), 1, max(tab1)+0.2])
    plt.axis([0, len(tab1), 0, 10]) # 15, 1, 5.5
    plt.legend()

#    sub1=fig_ave.add_subplot(2,2,1); plt.plot(tab1); sub1.set_title("on cells")
#    sub2=fig_ave.add_subplot(2,2,2); plt.plot(tab2); sub2.set_title("off cells")
#    sub3=fig_ave.add_subplot(2,2,3); plt.plot(tab3); sub3.set_title("on-off cells")
#    sub4=fig_ave.add_subplot(2,2,4); plt.plot(tab4); sub4.set_title("all cells")

    plt.show()

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    main(sys.argv[1])
    print("done")
else:
    exit("arg error - need 1 arg: [RI file]")

