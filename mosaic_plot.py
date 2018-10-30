#!/usr/bin/env python3
import sys, os, re
import matplotlib.pyplot as plt

#TODO: what is this scipt for?
#--------------------------------------------------------------------------#
def main(folder):
    for coord_file in os.listdir(folder):
        x=[]; y=[]
        if coord_file.endswith(".txt"):
            with open(folder+"/"+coord_file, 'r') as fichier:
                for line in fichier:
                    # x y
                    m=re.search( r'(.+) (.+)', line, re.M|re.I)
                    if m :
                        x.append(float(m.group(1)))
                        y.append(float(m.group(2)))
            myPlot(x, y, coord_file)
    if sys.argv[2]=="all":
        plt.show()

#--------------------------------------------------------------------------#
def myPlot(x, y, file_name):
    plt.scatter(x, y, 4)
    if sys.argv[2]=="-":
        plt.title(file_name)
        plt.show()

#--------------------------------------------------------------------------#
if len(sys.argv)!=3:
    exit("arg error - need 2 arg: [files directory] [all,-]")
else:
    if sys.argv[2]=="all" or sys.argv[2]=="-" :
        main(sys.argv[1])
        print("done")
    else:
        exit("arg error - need 2 arg: [files directory] [all,-]")

