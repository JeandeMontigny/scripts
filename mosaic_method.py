#!/usr/bin/env python3
import sys, os, re
import numpy as np
import random as rd
import scipy.spatial as ss
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
#TODO: several repetition for stats
def main(folder):
    create_mosaics = False
    figure_creation = True
    stat_analysis = True

    if create_mosaics:
        for rand in range(1, 10):
            mosaic_creation(folder, rand/10)

    delau_list = []; voro_list = []; ri_list = []; dist_list = []; random_weight_list = []

    for coord_file in [file for file in os.listdir(folder) if file.endswith(".txt")]:
        random_weight = coord_file[len(coord_file)-7:len(coord_file)-4] if coord_file[len(coord_file)-6] == '.' else coord_file[len(coord_file)-5:len(coord_file)-4]
        random_weight_list.append(float(random_weight))

        # read file
        cells_position = read(folder+coord_file)

        # data analyse
        delau_list.append(delaunay(folder, cells_position, random_weight))
        voro_list.append(voronoi(folder, cells_position, random_weight))
        ri_val, dist = ri(cells_position)
        ri_list.append(ri_val)
        dist_list.append(dist)

    random_weight_list, delau_list, voro_list, ri_list = sortData(random_weight_list, delau_list, voro_list, ri_list)

    if figure_creation:
        delauPlot(folder, random_weight_list, delau_list)
        voroPlot(folder, random_weight_list, voro_list)
        riPlot(folder, random_weight_list, ri_list)
        distPlot(folder, random_weight_list, dist_list)

    if stat_analysis:
    #TODO: cumulative distribution stat analysis.
    #      Kolmogorov-Smirnov test: statistic, pvalue = ks_2samp(distrib1, distrib2)
    #      analyse for method sensitivity and specificity?

        delayStats(random_weight_list, delau_list)
        voroStats(random_weight_list, voro_list)
        riStats(random_weight_list, ri_list)

    return 1

#--------------------------------------------------------------------------#
def mosaic_creation(output_folder, weight):
    if float(weight) > 1 or float(weight) < 0:
        SystemExit('random weigh should be between 0 and 1')
    cell_per_dim = 20
    cells_space = 25
    rand = float(weight) * cells_space

    if os.path.isfile(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+str(weight)+".txt"):
        os.remove(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+str(weight)+".txt")
    output_file=open(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+str(weight)+".txt", "a")

    for i in range(0, cell_per_dim):
        for j in range(0, cell_per_dim):
            position = str( round(rd.uniform(-rand, rand), 2)+(i*cells_space) ) + " " + str( round(rd.uniform(-rand, rand), 2)+(j*cells_space) ) + "\n"
            output_file.write(position)
    output_file.close()

#--------------------------------------------------------------------------#
def read(coord_file):
    cells_position = []
    file_lines = open(coord_file, "r").readlines()
    for line in file_lines:
        m = re.search(r'([0-9.]+)\ ([0-9.]+)\ *[0-9.]*\n', line, re.M|re.I)
        if m:
            cells_position.append([float(m.group(1)), float(m.group(2))])

    return cells_position

#--------------------------------------------------------------------------#
def delaunay(output_folder, positions_list, random_weight):
    tri = ss.Delaunay(positions_list, qhull_options="QJ")

    return tri

    plt.figure()
    plt.triplot(tri.points[:,0], tri.points[:,1], tri.simplices)
    plt.plot(tri.points[:,0], tri.points[:,1], 'o', color='black')
    plt.title("Delaunay triangulation")

    seg_length = []; seg_done = []
    # simplices: index of tri.points points forming triangles
    for triangle in tri.simplices:
        segments = [[tri.points[triangle[0]], tri.points[triangle[1]]], [tri.points[triangle[0]], tri.points[triangle[2]]], [tri.points[triangle[1]], tri.points[triangle[2]]]]
        for seg in segments:
            if not segDone(seg, seg_done):
                seg_length.append(dist(seg))
                seg_done.append(seg)

    # delaunay segment length density distribution
    plt.figure()
    n = plt.hist(np.sort(seg_length)[:int(len(seg_length)-len(seg_length)*0.05)], bins=int(len(seg_length)/40), density=True, cumulative=False, histtype='bar')
    plt.title("Delaunay triangulation segment length density distribution")
    plt.savefig(output_folder+"delau_seg_distrib_"+random_weight+".png")

    # delaunay segment length cumulative density
    density = []
    for i in range(0, len(n[0])):
        density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
    plt.figure()
    plt.plot(n[1][:len(n[1])-1], density/sum(n[0]))
    plt.title("Delaunay triangulation segment length cumulative density")
    plt.savefig(output_folder+"delau_seg_cumul_"+random_weight+".png")

#--------------------------------------------------------------------------#
def segDone(new_seg, seg_list):
    for seg in seg_list:
        if seg[0][0] == new_seg[0][0] and seg[0][1] == new_seg[0][1] and seg[1][0] == new_seg[1][0] and seg[1][1] == new_seg[1][1]\
            or seg[0][0] == new_seg[1][0] and seg[0][1] == new_seg[1][1] and seg[1][0] == new_seg[0][0] and seg[1][1] == new_seg[0][1]:
            return True
    return False

#--------------------------------------------------------------------------#
def dist(points):
    return round(np.sqrt(np.square(points[0][0] - points[1][0]) + np.square(points[0][1] - points[1][1])), 2)

#--------------------------------------------------------------------------#
def voronoi(output_folder, positions_list, random_weight):
    voro = ss.Voronoi(positions_list, qhull_options="Qc")

    return voro

    # voro.vertices: delaunays circles centers
    # voro.ridge_points: index of voro.vertices points couple forming lines
    # voro.ridge_vertices: index of voro.vertices points couple forming lines, including outside of region as -1
    # voro.regions: index of vertices points composing regions. closed region if no -1 index
    ss.voronoi_plot_2d(voro)

    # calculate areas and angles of every internal domains
    areas_list = []; angles_list = []
    for region in voro.regions:
        # if region doesn't have point outside. ignore non closed domains (border domains), avoid border effect
        if -1 not in region:
            domain_points = []
            for index in region:
                domain_points.append(list(voro.vertices[index]))
            angles_list.append(polygoneAngles(domain_points))
            areas_list.append(polygoneArea(domain_points))

    # voronoi area density distribution
    plt.figure()
    n = plt.hist(np.sort(areas_list)[:int(len(areas_list)-len(areas_list)*0.1)], bins=int(len(areas_list)/10), density=True, cumulative=False, histtype='bar')
    plt.title("Voronoi domains area density distribution")
    plt.savefig(output_folder+"voro_area_distrib_"+random_weight+".png")

    # voronoi area cumulative density
    density = []
    for i in range(0, len(n[0])):
        density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
    plt.figure()
    plt.plot(n[1][:len(n[1])-1], density/sum(n[0]))
    plt.title("Voronoi domains area cumulative density")
    plt.savefig(output_folder+"voro_area_cumul_"+random_weight+".png")

    # voronoi angles distribution
    angles = []
    for angle_sub_list in angles_list:
        for angle in angle_sub_list:
            angles.append(np.degrees(angle))
    plt.figure()
    n = plt.hist(angles, bins=int(len(angles)/100), density=True, cumulative=False, histtype='bar')
    plt.title("Voronoi domains angle density distribution")
    plt.savefig(output_folder+"voro_angle_distrib_"+random_weight+".png")

    # voronoi angles cumulative density
    density = []
    for i in range(0, len(n[0])):
        density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
    plt.figure()
    plt.plot(n[1][:len(n[1])-1], density/sum(n[0]))
    plt.title("Voronoi domains angle cumulative density")
    plt.savefig(output_folder+"voro_angle_cumul_"+random_weight+".png")

#--------------------------------------------------------------------------#
def polygoneAngles(points_coord):
    angle_list = []
    for i in range(0, len(points_coord)):
        if i == 0:
            angle_list.append(getAngle(points_coord[i+1], points_coord[i], points_coord[len(points_coord)-1]))
        elif i == len(points_coord)-1:
            angle_list.append(getAngle(points_coord[0], points_coord[i], points_coord[i-1]))
        else:
             angle_list.append(getAngle(points_coord[i+1], points_coord[i], points_coord[i-1]))

    return angle_list

#--------------------------------------------------------------------------#
def getAngle(b, a, c):
    angle = (np.square(dist([a, b])) + np.square(dist([a, c])) - np.square(dist([b, c]))) / (2 * dist([a, b]) * dist([a, c]))
    if angle < -1:
        angle = -1
    if angle > 1:
        angle = 1

    return np.arccos(angle)

#--------------------------------------------------------------------------#
def polygoneArea(points_coord):
    x = []; y = []
    for coord in points_coord:
        x.append(coord[0])
        y.append(coord[1])

    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

#--------------------------------------------------------------------------#
def ri(positions_list):
    shortest_dist_list, dist_list = getDistLists(positions_list)

    return round((np.average(shortest_dist_list)/np.std(shortest_dist_list)), 2), dist_list

#--------------------------------------------------------------------------#
def getDistLists(coord_list):
    shortest_dist_list = []; dist_list = []
    for i in range(0, len(coord_list)):
        distance_list = []
        cell_coord = coord_list[i]
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            tempsDistance = np.sqrt(np.square(cell_coord[0] - other_cell_coord[0]) + np.square(cell_coord[1] - other_cell_coord[1]))
            # if cell is not itself
            if tempsDistance != 0:
                distance_list.append(tempsDistance)
                dist_list.append(tempsDistance)
        # add shortest distance
        shortest_dist_list.append(min(distance_list))

    return shortest_dist_list, dist_list

#--------------------------------------------------------------------------#
def sortData(label, a, b, c):
    # insertion sort
    i = 1
    while i < len(label):
        j = i
        while j > 0 and label[j-1] > label[j]:
            label[j], label[j-1] = label[j-1], label[j]
            a[j], a[j-1] = a[j-1], a[j]
            b[j], b[j-1] = b[j-1], b[j]
            c[j], c[j-1] = c[j-1], c[j]
            j = j-1
        i = i+1

    return label, a, b, c

#--------------------------------------------------------------------------#
def delauPlot(folder, random_weight_list, delau_list):
    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    fig2 = plt.figure(); ax2=fig2.add_subplot(111)
    for tri_index in range(0, len(delau_list)):
        tri = delau_list[tri_index]
        seg_length = []; seg_done = []
        # simplices: index of tri.points points forming triangles
        for triangle in tri.simplices:
            segments = [[tri.points[triangle[0]], tri.points[triangle[1]]], [tri.points[triangle[0]], tri.points[triangle[2]]], [tri.points[triangle[1]], tri.points[triangle[2]]]]
            for seg in segments:
                if not segDone(seg, seg_done):
                    seg_length.append(dist(seg))
                    seg_done.append(seg)

        # delaunay segment length density distribution
        n = ax1.hist(np.sort(seg_length)[:int(len(seg_length)-len(seg_length)*0.05)], bins=int(len(seg_length)/40), density=True, cumulative=False, histtype='bar')

        # delaunay segment length cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
        ax2.plot(n[1][:len(n[1])-1], density/sum(n[0]), label=str(random_weight_list[tri_index]))

    ax1.set_title("Delaunay triangulation segment length density distribution")
    ax1.get_figure().savefig(folder+"delau_seg_distrib_vs_rand.png")

    ax2.set_title("Delaunay triangulation segment length cumulative density")
    ax2.legend()
    ax2.get_figure().savefig(folder+"delau_seg_cumul_vs_rand.png")

#--------------------------------------------------------------------------#
def voroPlot(folder, random_weight_list, voro_list):
    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    fig2 = plt.figure(); ax2=fig2.add_subplot(111)
    fig3 = plt.figure(); ax3=fig3.add_subplot(111)
    fig4 = plt.figure(); ax4=fig4.add_subplot(111)
    for voro_index in range(0, len(voro_list)):
        areas_list = []; angles_list = []
        voro = voro_list[voro_index]
        for region in voro.regions:
            # if region doesn't have point outside. ignore non closed domains (border domains), avoid border effect
            if -1 not in region:
                domain_points = []
                for index in region:
                    domain_points.append(list(voro.vertices[index]))
                angles_list.append(polygoneAngles(domain_points))
                areas_list.append(polygoneArea(domain_points))
        # voronoi area density distribution
        n = ax1.hist(np.sort(areas_list)[:int(len(areas_list)-len(areas_list)*0.1)], bins=int(len(areas_list)/10), density=True, cumulative=False, histtype='bar')
        # voronoi area cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
        ax2.plot(n[1][:len(n[1])-1], density/sum(n[0]), label=str(random_weight_list[voro_index]))
        # voronoi angles distribution
        angles = []
        for angle_sub_list in angles_list:
            for angle in angle_sub_list:
                angles.append(np.degrees(angle))
        n = ax3.hist(angles, bins=int(len(angles)/100), density=True, cumulative=False, histtype='bar')
        # voronoi angles cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
        ax4.plot(n[1][:len(n[1])-1], density/sum(n[0]), label=str(random_weight_list[voro_index]))

    ax1.set_title("Voronoi domains area density distribution")
    ax1.get_figure().savefig(folder+"voro_area_distrib_vs_rand.png")

    ax2.set_title("Voronoi domains area cumulative density")
    ax2.legend()
    ax2.get_figure().savefig(folder+"voro_area_cumul_vs_rand.png")

    ax3.set_title("Voronoi domains angles density distribution")
    ax3.get_figure().savefig(folder+"voro_angle_distrib_vs_rand.png")

    ax4.set_title("Voronoi domains angles cumulative density")
    ax4.legend()
    ax4.get_figure().savefig(folder+"voro_angle_cumul_vs_rand.png")

#--------------------------------------------------------------------------#
def riPlot(folder, random_weight_list, ri):
    plt.figure()
    plt.plot(random_weight_list, ri)
    plt.title("RI measure depending on randomness weight for mosaic creation")
    plt.savefig(folder+"ri_vs_random.png")

#--------------------------------------------------------------------------#
def distPlot(folder, random_weight_list, voro_list):
    #TODO
    return

#--------------------------------------------------------------------------#
def delayStats(random_weight_list, delau_list):

    return

#--------------------------------------------------------------------------#
def voroStats(random_weight_list, voro_list):

    return

#--------------------------------------------------------------------------#
def riStats(random_weight_list, ri_list):

    return

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 2 arg: [input-output folder]')
