#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(surface_folders, density_folders):
    label = []; global_cell_nb_ave = []; global_cell_nb_std = []; cell_pop_list = []
    # get folder corresponding to the same day for surface and density
    for surface_day in os.listdir(surface_folders):
        for density_day in os.listdir(density_folders):
            if density_day == surface_day:
                day_global_cell_nb = []
                # extract surface and name for each retina of this day
                temps = getSurface(surface_folders+surface_day); ret_surface_list = temps[0]; retina_list = temps[1]
                # for each whole mount retina of this dev day
                for i in range(0, len(retina_list)-1):
                    ret_density = getDensity(density_folders+density_day, retina_list[i])
                    # if this whole mount retina had density measures
                    if ret_density > 0:
                        day_global_cell_nb.append(ret_density*ret_surface_list[i])

                label.append(int(density_day[1:]))
                cell_pop_list.append(day_global_cell_nb)
                global_cell_nb_ave.append(np.average(day_global_cell_nb))
                global_cell_nb_std.append(np.std(day_global_cell_nb))

    # hard copy to avoid referenced copy
    unsorted_label = label.copy()
    plot(sortData(label, global_cell_nb_ave, global_cell_nb_std), "Postnatal day", "Estimated cell population")
    apoptosePlot(sortBis(unsorted_label, cell_pop_list))

    plt.show()

    return 1

#--------------------------------------------------------------------------#
def getSurface(folder):
    list_file = [name for name in os.listdir(folder) if not os.path.isdir(os.path.join(folder, name))]
    areas = []; ret_name = []
    if len(list_file) == 0:
        return 0
    for file in list_file:
        # extract retina name
        ret_name.append(file[:12])
        # open and read file
        file_lines = open(folder+"/"+file, "r").readlines()
        m = re.search(r'1\t([0-9.]+).+', file_lines[1], re.M|re.I)
        areas.append(float(m.group(1))/1e6 if m else 0)

    return areas, ret_name

#--------------------------------------------------------------------------#
def getDensity(folder, retina):
    nb_of_cells = 0
    # get all RGC cell position files, corresponding only to this precise whole mount retina 
    files = [f for f in os.listdir(folder) if f.endswith("RGC.txt") and f.startswith(retina)]
    # if no cell position corresponsing to this retina, return
    if len(files) == 0:
        return 0
    # for each cell position file of files list
    for coord_file in files:
        file_lines = open(folder+"/"+coord_file, "r").readlines()
        for line in file_lines:
            if line[0] != " ":
                # file strucure: index Area Mean Min Max X Y
                m = re.search( r'[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t([0-9.]+)\t([0-9.]+)', line, re.M|re.I)
                if m :
                    nb_of_cells += 1
    # return average number of cell for this retina
    return ((nb_of_cells / (335.4*335.4)) * 10e5)/len(files)

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
def sortBis(label, pop_list):
    i = 1
    while i < len(label):
        j = i
        while j > 0 and label[j-1] > label[j]:
            label[j], label[j-1] = label[j-1], label[j]
            pop_list[j], pop_list[j-1] = pop_list[j-1], pop_list[j]
            j = j-1
        i = i+1

    return label, pop_list

#--------------------------------------------------------------------------#
def plot(data, x_label, y_label):
    plt.figure()
    plt.errorbar(data[0], data[1], data[2])
    plt.xlabel(x_label)
    plt.ylabel(y_label)

#--------------------------------------------------------------------------#
def apoptosePlot(data):
    cell_pop_list = data[1]
    estimated_apoptose = []; estimated_apoptose_std = []
    # for each development day
    for i in range(1, len(cell_pop_list)):
        day_temps = []
        # for each population measure of day x+1
        for j in range(0, len(cell_pop_list[i])):
            # for each population measure of day x
            for k in range(0, len(cell_pop_list[i-1])):
                # measure j day x+1 - measure k day x
                estimated_apo = cell_pop_list[i][j] - cell_pop_list[i-1][k]
                day_temps.append(estimated_apo)

        estimated_apoptose.append(np.average(day_temps))
        estimated_apoptose_std.append(np.std(day_temps))


    cumulative_death = []; cumulative_death_std = []
    max_val = max(estimated_apoptose); min_val = min(estimated_apoptose)
    for i in range(0, len(estimated_apoptose)):
        cumulative_death.append(100*(estimated_apoptose[i] - min_val) / (max_val - min_val))
        cumulative_death_std.append(100*(estimated_apoptose_std[i] - max_val) / (min_val - max_val))

    plot([data[0][:-1], estimated_apoptose, estimated_apoptose_std], "Postnatal day", "Estimated apoptosis")
    plot([data[0][:-1], cumulative_death, cumulative_death_std], "Postnatal day", "Cumulative apoptosis (%)")

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==3:
    if main(sys.argv[1], sys.argv[2]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [surface meta folder] [density meta folder]')
