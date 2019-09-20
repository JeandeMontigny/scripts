#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(fodlers):
    label = []; mean_area = []; std_area = []

    list_dir = [name for name in os.listdir(fodlers) if os.path.isdir(os.path.join(fodlers, name))]
    for folder in list_dir:
        temps = extractFolderData(fodlers+folder)
        if temps:
            label.append(temps[0])
            mean_area.append(temps[1])
            std_area.append(temps[2])

    plot(mean_area, std_area, label)

    return 1

#--------------------------------------------------------------------------#
def extractFolderData(folder):
    list_file = [name for name in os.listdir(folder) if not os.path.isdir(os.path.join(folder, name))]
    areas = []
    if len(list_file) == 0:
        return 0
    # extract folder corresponding developmental day. assume it won't be higher that P99
    label = folder[-2:] if folder[-2:-1] == 'P' else folder[-3:]
    for file in list_file:
        areas.append(readFile(folder+"/"+file))

    return label, round(np.average(areas), 2), round(np.std(areas), 2)

#--------------------------------------------------------------------------#
def readFile(file):
    file_open = open(file, "r")
    file_lines = file_open.readlines()
    file_lines[1]
    m = re.search(r'1\t([0-9.]+).+', file_lines[1], re.M|re.I)

    return float(m.group(1))/1e6 if m else 0

#--------------------------------------------------------------------------#
def plot(mean, std, label):
#    x = [value for value in range(0, len(mean))]
    plt.errorbar(label, mean, std)
    plt.xlabel("Postnatal day")
    plt.ylabel("Retinal surface")

    plt.show()

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [surface meta folder]')
