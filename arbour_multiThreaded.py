#!/usr/bin/env python3
import sys, os, re
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

import threading
from queue import Queue
import time

#--------------------------------------------------------------------------#
def main(fodler):
    number_of_file = get_nb_files(fodler)
    print(str(number_of_file) + " files to process")

    output_tab = []

    def process_queue():
        while True:
            current_file = file_queue.get()
            output_tab.append(process_file(read_file(current_file)))
            file_queue.task_done()

    file_queue = Queue()
    threads_nb = os.cpu_count()

    for i in range(32):
        t = threading.Thread(target=process_queue)
        t.daemon = True
        t.start()

    start = time.time()

    for file in os.listdir(fodler):
        if file.endswith(".swc"):
            file_queue.put(fodler+file)

    file_queue.join()

    # print(threading.enumerate())
    print("Execution time: " + str(time.time() - start))

    return 1

#--------------------------------------------------------------------------#
def read_file(file):
    nb_dend_root = 0; nb_dend_seg = 0; line_count = 0; line_non_point = 0; prev_parent = 0
    coord_soma = [0, 0, 0]; coord_point = [0, 0, 0]
    coord_tab = []; distance_branching_tab = []; distance_terminal_tab = []
    branching_index = []; terminasion_index = []

    for line in open(file, "r").readlines():
        m = re.search( r' ?([0-9]+) ([0-9]+) (.+) (.+) (.+) (.+) ([-0-9]+).?', line, re.M|re.I)

        if m and (line[0] != "#" or line[0] != "\n"):
            line_count+=1
            point_id = int(m.group(1)); point_type = int(m.group(2));
            point_x = float(m.group(3)); point_y = float(m.group(4)); point_z = float(m.group(5))
            point_radius = float(m.group(6)); point_parent = int(m.group(7))

            # if this is the soma
            if point_parent == -1:
                coord_soma = [point_x, point_y, point_z]
            # if parent is the soma
            if point_parent == 1:
                nb_dend_root+=1
            # if point is a dendrite or a branching of a terminal
            if point_type == 3 or point_type == 5 or point_type == 6:
                nb_dend_seg+=1
                coord_point = [point_x, point_y, point_z]
                coord_tab.append([point_x - coord_soma[0], point_y - coord_soma[1], point_z - coord_soma[2]])
                # if point is a branching
                if point_type == 5:
                    distance_branching_tab.append(distance_3d(coord_point, coord_soma))
                # if point is a terminal
                if point_type == 6:
                    distance_terminal_tab.append(distance_3d(coord_point, coord_soma))

                # if the point is not a son of the previous point and don't belong to the soma
            if point_parent < prev_parent and point_parent > 1:
                # it means that the parent of that point is a branching point
                branching_index.append(point_parent)
                # and that the previous point is a terminaison point
                terminasion_index.append(prev_id)

            prev_parent = point_parent
            prev_id = point_id

        elif line[0] == "#" or line[0] == "\n":
            line_non_point+=1

    if nb_dend_seg == 0 or len(distance_branching_tab) == 0:
        warn("file only contains a cell body or no branching points")
        return None

    return distance_terminal_tab, distance_branching_tab, coord_tab

#--------------------------------------------------------------------------#
def process_file(output):
    # output: distance_terminal_tab, distance_branching_tab, coord_tab
    average_terminal_distance=round(np.average(output[0]), 3)
    average_terminal_distance_std=round(np.std(output[0]), 3)
    disc_diam_95=disc_span_95(output[0])
    average_branch_distance=round(np.average(output[1]), 3)
    average_branch_distance_std=round(np.std(output[1]), 3)
    aniso_score = anisometry(output[2])
    z_coord_distrib = get_z_distrib(output[2])
    peaks = peak_detector(z_coord_distrib)
    cell_type = type_finder(peaks)

    return

#--------------------------------------------------------------------------#
def distance_3d(point1, point2):
    if len(point1) != 3 or len(point2) != 3:
        raise SystemExit('Error: can\'t calculate distance if point are not in 3D')

    return round(np.sqrt(np.square(point1[0] - point2[0]) + np.square(point1[1] - point2[1]) + np.square(point1[2] - point2[2])), 3)

#--------------------------------------------------------------------------#
def disc_span_95(tab):
    """ Return the distance from soma containing 95% of terminal dendrites (mean + 2 std).
    @params:
        list of terminal point distance from soma """
    tab_length=len(tab); tab_95=[]
    tab.sort()
    for dist in tab:
        if len(tab_95)<0.95*tab_length:
            tab_95.append(dist)
        else:
            return max(tab_95)

#--------------------------------------------------------------------------#
def anisometry(coord_tab):
    tab = []
    for coord in coord_tab:
        tab.append(iso_distributor(coord))

    theo_val = float(len(tab)/8)
    if len(tab) < 2 or theo_val == 0:
        return 0
    return round(((abs(tab.count(1)-theo_val)/theo_val)+(abs(tab.count(2)-theo_val)/theo_val)+(abs(tab.count(3)-theo_val)/theo_val)+(abs(tab.count(4)-theo_val)/theo_val)+(abs(tab.count(5)-theo_val)/theo_val)+(abs(tab.count(6)-theo_val)/theo_val)+(abs(tab.count(7)-theo_val)/theo_val)+(abs(tab.count(8)-theo_val)/theo_val))/14, 3)

#--------------------------------------------------------------------------#
def iso_distributor(coord):
    x = coord[0]; y = coord[1]
    if x < 0:
        if y < 0:
            return 1 if abs(x) > abs(y) else 2
        else:
            return 3 if abs(x) > y else 4
    else:
        if y < 0:
            return 5 if x > abs(y) else 6
        else:
            return 7 if x > y else 8

#--------------------------------------------------------------------------#
def count(tab, inf_boundary, sup_boundary):
    count = 0
    for i in tab:
        if abs(i) >= inf_boundary and abs(i) < sup_boundary:
            count+=1
    return count

#--------------------------------------------------------------------------#
def count_z(multi_tab, inf_boundary, sup_boundary):
    count = 0
    for tab in multi_tab:
        if abs(tab[2]) >= inf_boundary and abs(tab[2]) < sup_boundary:
            count+=1
    return count

#--------------------------------------------------------------------------#
def get_z_distrib(coord_tab):
    tab_z_count = []; interval = []
    for i in frange(0,60, 0.5):
        tab_z_count.append(count_z(coord_tab, i, i+0.5))
        interval.append(i)
    tab_z_count = np.asarray(tab_z_count)

    return tab_z_count

#--------------------------------------------------------------------------#
def peak_detector(z_distrib, width_value=3):
    # if width is too low, return 0
    if width_value < 1.5:
        return 0
    peaks, _ = find_peaks(z_distrib, width=width_value)
    # if no peaks have been found, lower the width
    if len(peaks) == 0:
        return peak_detector(z_distrib, width_value-0.5)
    return peaks

#--------------------------------------------------------------------------#
def type_finder(peaks):
    if np.amax(peaks) < 50:
        return "on"
    if np.amin(peaks) > 50:
        return "off"
    if np.amin(peaks) < 50 and np.amax(peaks) > 50:
        return "on-off"
    else:
        return "other"

#--------------------------------------------------------------------------#
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump

#--------------------------------------------------------------------------#
def figure_construction():
    plt.figure(1)

    plt.ylim(1.5, 11)
    plt.xlabel("simulation time")
    plt.ylabel("regularity index")
#    plt.legend(fontsize="medium", ncol=2, bbox_to_anchor=(1.04,0.5), loc="center left", title="Density")
    plt.legend(loc=0)
    plt.title("layer colapse impact on RI")

#--------------------------------------------------------------------------#
def get_nb_files(folder):
    """ Return the number of files to process, including files in sub-folders.
    @params:
        folder containing sub-folders and files """
    nb_of_files = 0
    for name in os.listdir(folder):
        if name.endswith(".swc"):
            nb_of_files += 1

    return nb_of_files

#--------------------------------------------------------------------------#
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 0, length = 100, fill = "#"):
    """ Create progress bar
    @params:
        iteration   - Required: current iteration (Int)
        total       - Required: total iterations (Int)
        decimals    - Optional: positive number of decimals in percent complete (Int)
        length      - Optional: character length of bar (Int)
        fill        - Optional: bar fill character (Str)"""
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s [%s] %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print new line on complete
    if iteration == total:
        print("\n")

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [swc folder]')
