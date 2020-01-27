#!/usr/bin/env python3
import sys, os, re, time
from multiprocessing import Pool
import numpy as np
from scipy.signal import find_peaks
import scipy.cluster.vq as vq
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D
import scipy.signal as sig

#--------------------------------------------------------------------------#
def main(fodler):
    start = time.time()
    threads_nb = int(os.cpu_count()/2)
    file_queue = [fodler+file for file in os.listdir(fodler) if file.endswith(".swc")]

    print("Processing", len(file_queue), "files, using", threads_nb, "CPUs")
    # threads_nb processors will run read_file method for each file in file_queue
    output = regroup_data(Pool(threads_nb).map(read_file, enumerate(file_queue)))

    analyse(output, figures = True, clustering = False, charac = False)

    print("Execution time:", round(time.time() - start, 2), "sec")
    plt.show()

    return 1

#--------------------------------------------------------------------------#
def read_file(enumerate_obj):
    nb_dend_root = 0; nb_dend_seg = 0; line_count = 0; prev_parent = 0
    coord_soma = [0, 0, 0]; coord_point = [0, 0, 0]
    coord_tab = []; distance_branching_tab = []; distance_terminal_tab = []
    z_terminal_tab = []
    branching_index = []; terminasion_index = []

    for line in open(enumerate_obj[1], "r").readlines():
        m = re.search( r' ?([0-9]+) ([0-9]+) (.+) (.+) (.+) (.+) ([-0-9]+).?', line, re.M|re.I)

        if m and (line[0] != "#" or line[0] != "\n"):
            line_count+=1
            point_id = int(m.group(1)); point_type = int(m.group(2));
            point_x = float(m.group(3)); point_y = float(m.group(4)); point_z = abs(float(m.group(5)))
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
                    distance_branching_tab.append(distance_xd(coord_point, coord_soma))
                # if point is a terminal
                if point_type == 6:
                    z_terminal_tab.append(point_z - coord_soma[2])
                    distance_terminal_tab.append(distance_xd(coord_point, coord_soma))

            # if the point is not a son of the previous point or of the soma
            if point_parent < prev_parent and point_parent > 1:
                # it means that the parent of that point is a branching point
                branching_index.append(point_parent)
                # and that the previous point is a terminaison point
                terminasion_index.append(prev_id)

            prev_parent = point_parent
            prev_id = point_id

    if nb_dend_seg == 0 or len(distance_branching_tab) == 0:
        print("file contains only a cell body or no branching points")
        return None

    return process_file([distance_terminal_tab, distance_branching_tab, coord_tab, z_terminal_tab, nb_dend_root])

#--------------------------------------------------------------------------#
def process_file(output):
    average_terminal_distance = round(np.average(output[0]), 3)
    # average_terminal_distance_std = round(np.std(output[0]), 3)
    disc_diam_95 = disc_span_95(output[0])
    nb_branching = len(output[1])
    average_branch_distance = round(np.average(output[1]), 3)
    # average_branch_distance_std = round(np.std(output[1]), 3)
    aniso_score = anisometry(output[2])
    mass_center = get_mass_center(output[2])
    z_coord_distrib = get_z_distrib(output[2], smooth = True)
    peaks = peak_detector(z_coord_distrib)
    cell_type = type_finder(peaks)
    ave_term_lam_depth = np.average(output[3])

    return average_terminal_distance, disc_diam_95, average_branch_distance, \
        aniso_score, z_coord_distrib, peaks, ave_term_lam_depth, \
        nb_branching, output[4], mass_center

#--------------------------------------------------------------------------#
def analyse(output, figures = True, clustering = False, pca = False, charac = False, print_results = True):

    if print_results:
        root = output[8]; diam = output[1]
        branch_nb = output[7]; aniso = output[3]
        print(len(output[0]), "cells", \
        "\nnb root:", round(np.average(root), 2), \
            "- std:", round(np.std(root), 2), \
        "\ndiameter:", round(np.average(diam), 2), \
            "- std:", round(np.std(diam), 2), \
        "\nbranchin nb:", round(np.average(branch_nb), 2), \
            "- std:", round(np.std(branch_nb), 2), \
        "\naniso:", round(np.average(aniso), 2), \
            "- std:", round(np.std(aniso), 2)
        )

    if figures:
        figure_construction(output)

    if charac:
        subtype_charac(output, figures = True)

    if clustering:
        clusters = k_clustering([output[1], output[3], output[7], output[2]], \
            4, emfo = False, figure = True, x_label = "arbour diameter", \
            y_label = "aniso score", z_label = "branching number", title = "all")

        # clustering separated by types (on, off, on-off)
        data = split_by_type(output)
        clusters_1 = k_clustering([data[0][1], data[0][3], data[0][5], \
            data[0][2]], 4, emfo = False, figure = True, \
            x_label = "arbour diameter", y_label = "aniso score", \
            z_label = "branching number", title = "on")

        clusters_2 = k_clustering([data[1][1], data[1][3], data[1][5], \
            data[1][2]], 4, emfo = False, figure = True, \
            x_label = "arbour diameter", y_label = "aniso score", \
            z_label = "branching number", title = "off")

        clusters_3 = k_clustering([data[2][1], data[2][3], data[2][5], \
            data[2][2]], 4, emfo = False, figure = True, \
            x_label = "arbour diameter", y_label = "aniso score", \
            z_label = "branching number", title = "on-off")

#TODO: find new PCA module
#    if pca:
#        pca_analysis([output[0], output[1], output[2], output[3]])

#--------------------------------------------------------------------------#
def subtype_charac(output, figures = False):
    on, off, on_off = split_by_type(output)

    # average on z terminal depth: 14.06 - z center of mass 13.78
    # average off z terminal depth: 35.28 - z center of mass 33.37
    # ave on-off terminal depth: 28.79 - z center of mass 31.32d

    diam = []; root = []; aniso = []; branch_nb = []
    # on-off DSGC (all 4 types in 1 category)
    # for i in range(0, len(on_off[1])):
    #     if on_off[4][i] > 26.1 and on_off[1][i] is not None :

    for i in range(0, len(output[1])):
        if output[1][i] is not None:
        # if output[1][i] is not None \
        #     and output[1][i] < 110  \
        #     and output[9][i][2] > 28 and output[3][i] > 0.5:
            diam.append(output[1][i])
            aniso.append(output[3][i])
            branch_nb.append(output[7][i])
            root.append(output[8][i])

    # for i in range(0, len(on_off[1])):
    #     if on_off[1][i] is not None:
    #     # if output[1][i] is not None \
    #     #     and output[1][i] < 110  \
    #     #     and output[9][i][2] > 28 and output[3][i] > 0.5:
    #         diam.append(on_off[1][i])
    #         aniso.append(on_off[3][i])
    #         branch_nb.append(on_off[5][i])
    #         root.append(on_off[6][i])

    print(len(diam), "mini-led", \
    "\nnb root:", round(np.average(root), 2), \
        "- std:", round(np.std(root), 2), \
    "\ndiameter:", round(np.average(diam), 2), \
        "- std:", round(np.std(diam), 2), \
    "\nbranchin nb:", round(np.average(branch_nb), 2), \
        "- std:", round(np.std(branch_nb), 2), \
    "\naniso:", round(np.average(aniso), 2), \
        "- std:", round(np.std(aniso), 2)
    )

    if figures:
        fig_violin(root, title = "root number distribution")
        fig_violin(diam, title = "diameter distribution")
        fig_violin(branch_nb, title = "branchin number distribution")
        fig_violin(aniso, title = "aniso score distribution")

#--------------------------------------------------------------------------#
def distance_xd(point1, point2):
    if len(np.shape(point1)) == 0 and len(np.shape(point2)) == 0:
        return round(np.sqrt(np.square(point1 - point2)), 3)

    if len(point1) != len(point2):
        raise SystemExit('Error: can\'t calculate distance if points don\'t have the same number of dimensions')

    sum_square = 0
    for dim in range(0, len(point1)):
        sum_square += np.square(point1[dim] - point2[dim])

    return round(np.sqrt(sum_square), 3)

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
            break
    return max(tab_95)

#--------------------------------------------------------------------------#
def anisometry(coord_tab):
    tab = [iso_distributor(coord) for coord in coord_tab]

    theo_val = float(len(tab)/8)
    if len(tab) < 2 or theo_val == 0:
        return 0
    return round( ( \
        (abs(tab.count(1) - theo_val)/theo_val) + \
        (abs(tab.count(2) - theo_val)/theo_val) + \
        (abs(tab.count(3) - theo_val)/theo_val) + \
        (abs(tab.count(4) - theo_val)/theo_val) + \
        (abs(tab.count(5) - theo_val)/theo_val) + \
        (abs(tab.count(6) - theo_val)/theo_val) + \
        (abs(tab.count(7) - theo_val)/theo_val) + \
        (abs(tab.count(8) - theo_val)/theo_val) \
        ) /14, 3)

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
def get_mass_center(coord_tab):
    x_center = []; y_center = []; z_center = []
    for point_coord in coord_tab:
        x_center.append(point_coord[0])
        y_center.append(point_coord[1])
        z_center.append(point_coord[2])

    return [round(np.average(x_center), 2), \
        round(np.average(y_center), 2), \
        round(np.average(z_center), 2)]

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
def get_z_distrib(coord_tab, smooth = False, figures = False):
    tab_z_count = []; interval = []
    for i in frange(0, 75, 0.5):
        tab_z_count.append(count_z(coord_tab, i, i+0.5))
        interval.append(i)
    tab_z_count = np.asarray(tab_z_count)
    if smooth:
        tab_z_count = sig.savgol_filter(tab_z_count, 5, 3)

    if figures:
        plt.figure()
        plt.plot(tab_z_count)
        plt.show()

    return tab_z_count

#--------------------------------------------------------------------------#
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump

#--------------------------------------------------------------------------#
def peak_detector(z_distrib, width_value = 3, plot = False):
    # if width is too low for result to be meaningful, return 0
    if width_value < 1.5:
        return 0
    peaks, _ = find_peaks(z_distrib, height = max(z_distrib)/4, width = width_value)
    # if no peaks have been found, lower the width
    if len(peaks) == 0:
        return peak_detector(z_distrib, width_value-0.5)
    if plot:
        plt.figure()
        plt.plot(z_distrib)
        for point in peaks:
            plt.plot(point, max(z_distrib)*1.05, 'ro')
        plt.title(type_finder(peaks))
        plt.show()

    return peaks

#--------------------------------------------------------------------------#
def type_finder(peaks):
    # if just one lamination peak
    if len(peaks) == 1:
        return "on" if peaks < 42 else "off" # 50
    # if distance between peaks is higher to be two lamination levels
    elif np.amax(peaks) - np.amin(peaks) >= 20:
        return "on-off"
    elif np.amin(peaks) < 42: # 50
        return "on"
    elif np.amax(peaks) >= 42: # 50
        return "off"
    else:
        return "other"

#--------------------------------------------------------------------------#
def split_by_type(output, print_types_nb = False):
    mean_z_terminal_on = []; aniso_sub_on = []; diam_sub_on = []
    term_dist_sub_on = []; branching_dist_sub_on = []
    z_coord_distrib_on = []
    branching_nb_on = []; root_nb_on = []; mass_center_on = []

    mean_z_terminal_off = []; aniso_sub_off = []; diam_sub_off = []
    term_dist_sub_off = []; branching_dist_sub_off = []
    z_coord_distrib_off = []
    branching_nb_off = []; root_nb_off = []; mass_center_off = []

    mean_z_terminal_on_off = []; aniso_sub_on_off = []; diam_sub_on_off = []
    term_dist_sub_on_off = []; branching_dist_sub_on_off = []
    z_coord_distrib_on_off = []
    branching_nb_on_off = []; root_nb_on_off = []; mass_center_on_off = []

    for i in range(0, len(output[5])):
        type = type_finder(output[5][i])
        if type == "on":
            term_dist_sub_on.append(output[0][i])
            diam_sub_on.append(output[1][i])
            branching_dist_sub_on.append(output[2][i])
            aniso_sub_on.append(output[3][i])
            z_coord_distrib_on.append(output[4][i])
            mean_z_terminal_on.append(output[6][i])
            branching_nb_on.append(output[7][i])
            root_nb_on.append(output[8][i])
            mass_center_on.append(output[9][i])
        elif type == "off":
            term_dist_sub_off.append(output[0][i])
            diam_sub_off.append(output[1][i])
            branching_dist_sub_off.append(output[2][i])
            aniso_sub_off.append(output[3][i])
            z_coord_distrib_off.append(output[4][i])
            mean_z_terminal_off.append(output[6][i])
            branching_nb_off.append(output[7][i])
            root_nb_off.append(output[8][i])
            mass_center_off.append(output[9][i])
        elif type == "on-off":
            term_dist_sub_on_off.append(output[0][i])
            diam_sub_on_off.append(output[1][i])
            branching_dist_sub_on_off.append(output[2][i])
            aniso_sub_on_off.append(output[3][i])
            z_coord_distrib_on_off.append(output[4][i])
            mean_z_terminal_on_off.append(output[6][i])
            branching_nb_on_off.append(output[7][i])
            root_nb_on_off.append(output[8][i])
            mass_center_on_off.append(output[9][i])

    if print_types_nb:
        print(len(term_dist_sub_on), "on cells,", len(term_dist_sub_off),\
            "off cells,", len(term_dist_sub_on_off), "on-off cells" )

    on = [term_dist_sub_on, diam_sub_on, branching_dist_sub_on, \
        aniso_sub_on, mean_z_terminal_on, branching_nb_on, root_nb_on, \
        z_coord_distrib_on, mass_center_on]
    off = [term_dist_sub_off, diam_sub_off, branching_dist_sub_off, \
        aniso_sub_off, mean_z_terminal_off, branching_nb_off, root_nb_off, \
        z_coord_distrib_off, mass_center_off]
    on_off = [term_dist_sub_on_off, diam_sub_on_off, \
        branching_dist_sub_on_off, aniso_sub_on_off, mean_z_terminal_on_off, \
        branching_nb_on_off, root_nb_on_off, z_coord_distrib_on_off, \
        mass_center_on_off]

    return on, off, on_off

#--------------------------------------------------------------------------#
def regroup_data(tab):
    terminal_dist = []; disc_diam = []; branch_dist = []; aniso = []
    z_coord_distrib = []; peaks = []; z_terminal_coord_distrib = []
    nb_branching = []; nb_root = []; mass_center = []
    for cell in tab:
        if cell is not None:
            terminal_dist.append(cell[0]); disc_diam.append(cell[1])
            branch_dist.append(cell[2]); aniso.append(cell[3])
            z_coord_distrib.append(cell[4]); peaks.append(cell[5])
            z_terminal_coord_distrib.append(cell[6])
            nb_branching.append(cell[7]); nb_root.append(cell[8])
            mass_center.append(cell[9])

    return terminal_dist, disc_diam, branch_dist, aniso, z_coord_distrib, \
        peaks, z_terminal_coord_distrib, nb_branching, nb_root, mass_center

#--------------------------------------------------------------------------#
def figure_construction(tab):
    # fig_violin(tab[1], title = "diameter distribution")
    # fig_violin(tab[2], title = "branching distance distribution")
    # fig_violin(tab[3], title = "anisometry score distribution")

    data = split_by_type(tab, print_types_nb = True)
    # all_cells = [tab[0], tab[1], tab[2], tab[3], tab[6], tab[7]]
    # on_cells = data[0]; off_cells = data[1]; on_off_cells = data[2]
    # diameter, aniso, branching_nb
    all_cells = [tab[1], tab[3], tab[7]]
    on_cells = [data[0][1], data[0][3], data[0][5]]
    off_cells = [data[1][1], data[1][3], data[1][5]]
    on_off_cells = [data[2][1], data[2][3], data[2][5]]

    # print(len(on_cells[0]), "on cells", \
    # "\nnb root:", round(np.average(on_cells[6]), 2), \
    #     "- std:", round(np.std(on_cells[6]), 2), \
    # "\ndiameter:", round(np.average(on_cells[1]), 2), \
    #     "- std:", round(np.std(on_cells[1]), 2), \
    # "\nbranchin nb:", round(np.average(on_cells[5]), 2), \
    #     "- std:", round(np.std(on_cells[5]), 2), \
    # "\naniso:", round(np.average(on_cells[3]), 2), \
    #     "- std:", round(np.std(on_cells[3]), 2)
    # )

    # print(len(off_cells[0]), "off cells", \
    # "\nnb root:", round(np.average(off_cells[6]), 2), \
    #     "- std:", round(np.std(off_cells[6]), 2), \
    # "\ndiameter:", round(np.average(off_cells[1]), 2), \
    #     "- std:", round(np.std(off_cells[1]), 2), \
    # "\nbranchin nb:", round(np.average(off_cells[5]), 2), \
    #     "- std:", round(np.std(off_cells[5]), 2), \
    # "\naniso:", round(np.average(off_cells[3]), 2), \
    #     "- std:", round(np.std(off_cells[3]), 2)
    # )

    # fig_violin(on_cells[0], title = "terminal distance")
    # fig_violin(all_cells[1], title = "arbour diameter")
    # # fig_violin(on_cells[2], title = "branching distance")
    # fig_violin(all_cells[3], title = "anisometry score")
    # # fig_violin(on_cells[4], title = "mean terminal depth")
    # fig_violin(all_cells[5], title = "branching number")

    # print("all_cells")
    # print(all_cells)
    #
    # print("on_cells")
    # print(on_cells)
    #
    # print("off_cells")
    # print(off_cells)
    #
    # print("on_off_cells")
    # print(on_off_cells)

    multiple_fig_violin([all_cells, on_cells, off_cells, on_off_cells], \
        labels = ["", "all", "", "on", "", "off", "", "on-off"], \
        titles = ["arbour diameter", "anisometry score", "branching number"])

#--------------------------------------------------------------------------#
def fig_violin(tab, x_label = "", y_label = "", title = ""):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    violin_parts = ax.violinplot(tab, showmeans = True, showmedians = False)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('gray')
        pc.set_edgecolor('black')
    violin_parts['cbars'].set_color('black')
    violin_parts['cmins'].set_color('black')
    violin_parts['cmaxes'].set_color('black')
    violin_parts['cmeans'].set_color('black')
    violin_parts['cmeans'].set_linestyle('--')
    # violin_parts['cmedians'].set_color('dimgray')

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

#--------------------------------------------------------------------------#
def multiple_fig_violin(tab, labels, titles):
    if len(tab[0]) != len(tab[1]) or len(tab[1]) != len(tab[2]):
        print("error in multiple_fig_violin: lists should be of same size")
        return

    import dendrites_real_data
    real_all_cells, real_on_cells, real_off_cells, real_on_off_cells =\
     dendrites_real_data.get_real_data()
    real_cells =\
     [real_all_cells, real_on_cells, real_off_cells, real_on_off_cells]

    colours_labels = []; initialise_label = True
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        colours_labels.append((mpatches.Patch(color=color), label))

    for i in range(0, len(tab[0])):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # real data
        violin_parts = ax.violinplot([cell_type[i] for cell_type in real_cells], showmeans = True, showextrema = False)
        # keep just left half
        for b in violin_parts['bodies']:
            m = np.mean(b.get_paths()[0].vertices[:, 0])
            b.get_paths()[0].vertices[:, 0] =\
             np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
            # b.set_color('royalblue')
        b = violin_parts['cmeans']
        for j in range(0, len(b.get_paths())):
            m = np.mean(b.get_paths()[j].vertices[:, 0])
            b.get_paths()[j].vertices[:, 0] =\
             np.clip(b.get_paths()[j].vertices[:, 0], -np.inf, m)
            b.set_color('blue')
        if initialise_label:
            add_label(violin_parts, "In-vitro")

        # these cells
        violin_parts = ax.violinplot([cell_type[i] for cell_type in tab], showmeans = True, showextrema = False)
        # keep just right half
        for b in violin_parts['bodies']:
            m = np.mean(b.get_paths()[0].vertices[:, 0])
            b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
            # b.set_color('orange')
        b = violin_parts['cmeans']
        for j in range(0, len(b.get_paths())):
            m = np.mean(b.get_paths()[j].vertices[:, 0])
            b.get_paths()[j].vertices[:, 0] = np.clip(b.get_paths()[j].vertices[:, 0], m, np.inf)
            b.set_color('darkorange')
        if initialise_label:
            add_label(violin_parts, "In-silico")
        initialise_label = False

        plt.legend(*zip(*colours_labels))
        ax.set_xticklabels(labels)
        ax.set_title(titles[i])

#--------------------------------------------------------------------------#
def fig_plot(multi_tab, x_label = "", y_label = "", title = ""):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for tab in multi_tab:
        ax.plot(tab)

#--------------------------------------------------------------------------#
def figure_cloud_3d(tab, colours = None, x_label = "x-axis", y_label = "y-axis", z_label = "z-axis", title = ""):
    fig=plt.figure()
    ax=Axes3D(fig)
    ax.scatter(tab[0], tab[1], tab[2], c = colours, s = 100)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_title(title)

#--------------------------------------------------------------------------#
def get_colours(list):
    min_idx = min(list); max_idx = max(list)
    list_colours = [rgb(min_idx, max_idx, index) for index in list]

    return list_colours

#--------------------------------------------------------------------------#
def rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r

    return r/255, g/255, b/255

#--------------------------------------------------------------------------#
def k_clustering(tab, k_value = 4, emfo = False, figure = False, \
    x_label = "x-axis", y_label = "y-axis", z_label = "z-axis", title = ""):
    tab_ndarray =  np.array(tab).T
    whitened = vq.whiten(tab_ndarray, check_finite = False)
    centroids, labels = vq.kmeans2(data = whitened, k = k_value)

    if emfo:
        emfo_k(tab)
    if figure:
        colours = get_colours(labels)
        figure_cloud_3d([tab[0], tab[1], tab[2]], colours, x_label = x_label, y_label = y_label, z_label = z_label, title = title)

    return centroids, labels

#--------------------------------------------------------------------------#
def emfo_k(tab, min_k = 1, max_k = 15):
    """ Elbow Method For Optimal k """
    distortion = []
    data = [measure for measure in tab]
    for k in range(min_k, max_k):
        clusters = k_clustering(data, k)
        distortion.append(get_distortion(data, clusters))
    plt.figure()
    plt.plot([k for k in range(min_k, max_k)], distortion)
    plt.xlabel("Number of clusters (k)")
    plt.ylabel("Distortion score")
    plt.title("Elbow method for k-means clustering optimal cluster number")

#--------------------------------------------------------------------------#
def get_distortion(data, clusters):
    """ return the sum of squared distances from each point to its assigned centroid
    """
    list_dist = []
    if len(np.shape(data)) == 1:
        for i in range(0, len(data)):
            centroid_coord = clusters[0][clusters[1][i]]
            list_dist.append(np.square(distance_xd(data[i], centroid_coord)))
    else:
        for i in range(0, len(data[0])):
            point_coord = [measure[i] for measure in data]
            centroid_coord = clusters[0][clusters[1][i]]
            list_dist.append(np.square(distance_xd(point_coord, centroid_coord)))

    return sum(list_dist)

#--------------------------------------------------------------------------#
def pca_analysis(tab):
    param1 = np.array(tab[0]); param2 = np.array(tab[1])
    param3 = np.array(tab[2]); param4 = np.array(tab[3])
    matrix = []
    for i in range (0, len(param1)):
        matrix.append([float(param1[i]), float(param2[i]), float(param3[i]), float(param4[i])])
    matrix = np.array(matrix)
    results = PCA(matrix)
    print("PCA - variance percentages for each component:", results.fracs)

    plot3d_pca(results)

#--------------------------------------------------------------------------#
def plot3d_pca(results):
    x = []; y = []; z = []
    for item in results.Y:
        x.append(item[0])
        y.append(item[1])
        z.append(item[2])
    figure_cloud_3d([x, y, z])

#--------------------------------------------------------------------------#
def print_progress_bar (iteration, total, prefix = '', suffix = '', decimals = 0, length = 100, fill = "#"):
    """ Create progress bar
    @params:
        iteration - Required: current iteration (Int)
        total     - Required: total iterations (Int)
        decimals  - Optional: positive number of decimals in percent complete (Int)
        length    - Optional: character length of bar (Int)
        fill      - Optional: bar fill character (Str)"""
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s [%s] %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print new line on complete
    if iteration == total:
        print("\n")

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2 or len(sys.argv)==3:
    path = sys.argv[1]
    if len(sys.argv)==3:
        if sys.argv[2] == "-r":
            import correction
            correction.main(path)
            path = path+"rewrite/"
        else:
            raise SystemExit('Error: invalid option. available option: [-r]')
    if main(path):
        print("done")
    else:
        print("error during execution, is", sys.argv[1], "a valid directory?")
else:
    raise SystemExit('Error: need at least 1 arg: [swc folder]; option: [-r]')
