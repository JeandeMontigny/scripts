#!/usr/bin/env python3
import sys, os, re
import numpy as np
import scipy.spatial as ss
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(coord_file):

    cells_position = read(coord_file)

    delaunay(cells_position)
    voronoi(cells_position)

    plt.show()

    return 1

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
def delaunay(positions_list):
    tri = ss.Delaunay(positions_list, qhull_options="QJ")

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

    plt.figure()
    n = plt.hist(seg_length, bins=int(len(seg_length)/10), density=False, cumulative=False, histtype='bar')
    plt.title("Delaunay triangulation segment length distribution")

    density = []
    for i in range(0, len(n[0])-1):
        density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
    plt.figure()
    plt.plot(density/sum(n[0]))
    plt.title("Delaunay triangulation segment length cumulative density")

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
def voronoi(positions_list):
    voro = ss.Voronoi(positions_list)

    ss.voronoi_plot_2d(voro)

    # voronoi domain size distribution

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [cells coordinates file]')
