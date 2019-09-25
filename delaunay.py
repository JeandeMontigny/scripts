#!/usr/bin/env python3
import sys, os, re
import numpy as np
import random as rd
import scipy.spatial as ss
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(coord_file):
    cells_position = read(coord_file)

    delaunay(cells_position)
    voronoi(cells_position)
    print(ri(cells_position))

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

    # delaunay segment length density distribution
    plt.figure()
    n = plt.hist(np.sort(seg_length)[:int(len(seg_length)-len(seg_length)*0.05)], bins=int(len(seg_length)/40), density=True, cumulative=False, histtype='bar')
    plt.title("Delaunay triangulation segment length density distribution")

    # delaunay segment length cumulative density
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
    voro = ss.Voronoi(positions_list, qhull_options="Qc")

    # voro.vertices: delaunays circles centers
    # voro.ridge_points: index of voro.vertices points couple forming lines
    # voro.ridge_vertices: index of voro.vertices points couple forming lines, including outside of region as -1
    # voro.regions: index of vertices points composing regions. closed region if no -1 index
    ss.voronoi_plot_2d(voro, show_vertices=False)

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

    # voronoi area cumulative density
    density = []
    for i in range(0, len(n[0])-1):
        density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
    plt.figure()
    plt.plot(density/sum(n[0]))
    plt.title("Voronoi domains area cumulative density")

    # voronoi angles distribution
    angles = []
    for angle_sub_list in angles_list:
        for angle in angle_sub_list:
            angles.append(np.degrees(angle))
    plt.figure()
    n = plt.hist(angles, bins=int(len(angles)/100), density=True, cumulative=False, histtype='bar')
    plt.title("Voronoi domains angle density distribution")

    # voronoi angles cumulative density
    density = []
    for i in range(0, len(n[0])-1):
        density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
    plt.figure()
    plt.plot(density/sum(n[0]))
    plt.title("Voronoi domains angle cumulative density")

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
    shortest_dist_list = getShortestDistList(positions_list)

    return round((np.average(shortest_dist_list)/np.std(shortest_dist_list)), 2)

#--------------------------------------------------------------------------#
def getShortestDistList(coord_list):
    shortest_dist_list = []
    for i in range(0, len(coord_list)):
        distance_list = []
        cell_coord = coord_list[i]
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            tempsDistance = np.sqrt(np.square(cell_coord[0] - other_cell_coord[0]) + np.square(cell_coord[1] - other_cell_coord[1]))
            # if cell is not itself
            if tempsDistance != 0:
                distance_list.append(tempsDistance)
        # add shortest distance
        shortest_dist_list.append(min(distance_list))

    return shortest_dist_list

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [cells coordinates file]')
