#!/usr/bin/env python3
import sys, os, re, time
from multiprocessing import Pool
import numpy as np
from scipy.signal import find_peaks
import scipy.cluster.vq as vq
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D
import scipy.signal as sig

#--------------------------------------------------------------------------#
def main(fodler):
    start = time.time()
    threads_nb = int(os.cpu_count()/2)
    file_queue = [fodler+file for file in os.listdir(fodler) if file.endswith(".swc")]

    print(len(file_queue), "files to process, using", threads_nb, "CPUs")
    # threads_nb processors will run read_file method for each file in file_queue
    output = regroup_data(Pool(threads_nb).map(read_file, enumerate(file_queue)))

    analyse(output, clustering = False)

    print("Execution time:", round(time.time() - start, 2), "sec")
    plt.show()

    return 1

#--------------------------------------------------------------------------#
def read_file(enumerate_obj):
    nb_dend_root = 0; nb_dend_seg = 0; line_count = 0; prev_parent = 0
    coord_soma = [0, 0, 0]; coord_point = [0, 0, 0]
    coord_tab = []; distance_branching_tab = []; distance_terminal_tab = []
    branching_index = []; terminasion_index = []

    for line in open(enumerate_obj[1], "r").readlines():
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
                    distance_branching_tab.append(distance_xd(coord_point, coord_soma))
                # if point is a terminal
                if point_type == 6:
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
        warn("file only contains a cell body or no branching points")
        return None

    return process_file([distance_terminal_tab, distance_branching_tab, coord_tab])

#--------------------------------------------------------------------------#
def process_file(output):
    # output: distance_terminal_tab, distance_branching_tab, coord_tab
    average_terminal_distance=round(np.average(output[0]), 3)
    average_terminal_distance_std=round(np.std(output[0]), 3)
    disc_diam_95=disc_span_95(output[0])
    average_branch_distance=round(np.average(output[1]), 3)
    average_branch_distance_std=round(np.std(output[1]), 3)
    aniso_score = anisometry(output[2])
    z_coord_distrib = get_z_distrib(output[2], smooth = True)
    peaks = peak_detector(z_coord_distrib)
    cell_type = type_finder(peaks)

    return average_terminal_distance, disc_diam_95, average_branch_distance, \
        aniso_score, z_coord_distrib, peaks

#--------------------------------------------------------------------------#
def analyse(output, figures = True, clustering = False, pca = False):
    if figures:
        figure_construction(output)

    if clustering:
        emfo_k(output)
        clusters = k_clustering([output[0], output[1], output[2], output[3]], 5)
        colours = [str(clust/5) for clust in clusters[1]]
        figure_cloud_3d([output[1], output[2], output[3]], colours)

    if pca:
        pca_analysis([output[0], output[1], output[2], output[3]])

    z_distrib = output[4]
    # TODO:analysis

#--------------------------------------------------------------------------#
def distance_xd(point1, point2):
    if len(point1) != len(point2):
        raise SystemExit('Error: can\'t calculate distance if points haven\'t got the same number of dimensions')

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
def get_z_distrib(coord_tab, smooth = False):
    tab_z_count = []; interval = []
    for i in frange(0,75, 0.5):
        tab_z_count.append(count_z(coord_tab, i, i+0.5))
        interval.append(i)
    tab_z_count = np.asarray(tab_z_count)
    if smooth:
        tab_z_count = sig.savgol_filter(tab_z_count, 5, 3)

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
    peaks, _ = find_peaks(z_distrib, height=max(z_distrib)/4, width=width_value)
    # if no peaks have been found, lower the width
    if len(peaks) == 0:
        return peak_detector(z_distrib, width_value-0.5)
    if plot:
        plt.figure()
        plt.plot(z_distrib)
        for point in peaks:
            plt.plot(point, max(z_distrib)*1.05, 'ro')
        plt.show()

    return peaks

#--------------------------------------------------------------------------#
def type_finder(peaks):
    if np.amax(peaks) < 59:
        return "on"
    if np.amin(peaks) > 59:
        return "off"
    if np.amin(peaks) < 59 and np.amax(peaks) > 59:
        return "on-off"
    else:
        return "other"

#--------------------------------------------------------------------------#
def regroup_data(tab):
    terminal_dist = []; disc_diam = []; branch_dist = []; aniso = []
    z_coord_distrib = []; peaks = []
    for cell in tab:
        if cell[1] is not None:
            terminal_dist.append(cell[0]); disc_diam.append(cell[1])
            branch_dist.append(cell[2]); aniso.append(cell[3])
            z_coord_distrib.append(cell[4]); peaks.append(cell[5])

    return terminal_dist, disc_diam, branch_dist, aniso, z_coord_distrib, peaks

#--------------------------------------------------------------------------#
def figure_construction(tab):
    # diam = tab[1]; branch = tab[2]; aniso = tab[3]
    # fig_violin(aniso, title = "anisometry score distribution")
    # fig_violin(diam, title = "diameter distribution")
    # fig_violin(branch, title = "branching distance distribution")

    z_distrib = tab[4]; peaks = tab[5]
    z_on = []; z_off = []; z_on_off = []
    for i in range(0, len(peaks)):
        if type_finder(peaks[i]) == "on":
            z_on.append(z_distrib[i])
        elif type_finder(peaks[i]) == "off":
            z_off.append(z_distrib[i])
        elif type_finder(peaks[i]) == "on-off":
            z_on_off.append(z_distrib[i])

    print(len(z_on), "on cells,", len(z_off), "off cells,", len(z_on_off), "on-off cells")

    plt.figure(97)
    for distrib in z_on:
        plt.plot(distrib)
        plt.xlim(0, 120)
        plt.title("on")
    plt.figure(98)
    for distrib in z_off:
        plt.plot(distrib)
        plt.xlim(0, 120)
        plt.title("off")
    plt.figure(99)
    for distrib in z_on_off:
        plt.plot(distrib)
        plt.xlim(0, 120)
        plt.title("on-off")

    #TODO: other measures plot
    # terminal_dist, disc_diam, branch_dist, aniso, z_coord_distrib, peaks

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
def fig_plot(multi_tab, x_label = "", y_label = "", title = ""):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for tab in multi_tab:
        ax.plot(tab)

#--------------------------------------------------------------------------#
def figure_cloud_3d(tab, colours = None, x_label = "x-axis", y_label = "y-axis", z_label = "z-axis", title = ""):
    fig=plt.figure()
    ax=Axes3D(fig)
    ax.scatter(tab[0], tab[1], tab[2], color = colours, s = 100)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_title(title)

#--------------------------------------------------------------------------#
def k_clustering(tab, k_value = 5):
    tab_ndarray =  np.array(tab).T
    whitened = vq.whiten(tab_ndarray, check_finite = False)
    centroids, labels = vq.kmeans2(data = whitened, k = k_value)

    return centroids, labels

#--------------------------------------------------------------------------#
def emfo_k(output, min_k = 1, max_k = 15):
    """ Elbow Method For Optimal k """
    distortion = []
    data = [output[0], output[1], output[2], output[3]]
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
    for i in range(0, len(data[0])):
        point_coord = [data[0][i], data[1][i], data[2][i], data[3][i]]
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
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [swc folder]')
