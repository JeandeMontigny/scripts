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
            label.append(temps[0]); mean_area.append(temps[1]); std_area.append(temps[2])

    plot(sortData(label, mean_area, std_area))

    return 1

#--------------------------------------------------------------------------#
def extractFolderData(folder):
    list_file = [name for name in os.listdir(folder) if not os.path.isdir(os.path.join(folder, name))]
    areas = []
    if len(list_file) == 0:
        return 0
    # extract folder corresponding developmental day. assume it won't be higher that P99
    label = folder[-1:] if folder[-2:-1] == 'P' else folder[-2:]
    for file in list_file:
        areas.append(readFile(folder+"/"+file))

    return int(label), round(np.average(areas), 2), round(np.std(areas), 2)

#--------------------------------------------------------------------------#
def readFile(file):
    file_lines = open(file, "r").readlines()
    m = re.search(r'1\t([0-9.]+).+', file_lines[1], re.M|re.I)

    return float(m.group(1))/1e6 if m else 0

#--------------------------------------------------------------------------#
def sortData(label, mean, std):
    # insertion sort
    i = 1
    while i < len(label):
        j = i
        while j > 0 and label[j-1] > label[j]:
            label[j], label[j-1] = label[j-1], label[j]
            mean[j], mean[j-1] = mean[j-1], mean[j]
            std[j], std[j-1] = std[j-1], std[j]
            j = j-1
        i = i+1

    return label, mean, std

#--------------------------------------------------------------------------#
def plot(data):
    plt.errorbar(data[0], data[1], data[2])
    plt.xlabel("Postnatal day")
    plt.ylabel("Retinal surface (mm$^2$)")

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
