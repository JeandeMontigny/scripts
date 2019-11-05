#!/usr/bin/env python3
import sys, os, re
import numpy as np
import scipy as sp
import random as rd
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(folder):
    create_mosaics = True
    mosaic_repetitions = 8
    figure_creation = True
    stat_analysis = True
    show_plots = False

    if create_mosaics:
        for rand in range(1, 10):
            mosaic_creation(folder, rand/10, mosaic_repetitions)

    delau_list = []; voro_list = []; ri_list = []; short_dist_list =[]; dist_list = []; random_weight_list = []

    for coord_file in [file for file in os.listdir(folder) if file.startswith("mosaic") and file.endswith(".txt")]:
        random_weight = coord_file[len(coord_file)-7:len(coord_file)-4] if coord_file[len(coord_file)-6] == '.' else coord_file[len(coord_file)-5:len(coord_file)-4]
        random_weight_list.append(float(random_weight))

        # read file
        cells_position = read(folder+coord_file)

        # data analyse
        delau_list.append(delaunay(folder, cells_position, random_weight))
        voro_list.append(voronoi(folder, cells_position, random_weight))
        ri_val, short_dist, dist = ri(cells_position)
        ri_list.append(ri_val); short_dist_list.append(short_dist); dist_list.append(dist)

    random_weight_list, delau_list, voro_list, ri_list, short_dist_list, dist_list = sortData(random_weight_list, delau_list, voro_list, ri_list, short_dist_list, dist_list)

    # average data
    weight_single_list = getWeightSingleList(random_weight_list, mosaic_repetitions)
    delau_x_cumul, delau_ave_cumul = delauAveCumul(weight_single_list, mosaic_repetitions, delau_list)
    area_x_cumul, area_ave_cumul, angle_x_cumul, angle_ave_cumul = voroAveCumul(random_weight_list, mosaic_repetitions, voro_list)
    ri_ave, ri_std, shortest_x_list, shortest_ave_list, dist_x_list, dist_ave_list = riAveCumul(random_weight_list, mosaic_repetitions, ri_list, short_dist_list, dist_list)

    if figure_creation:
        CumulDensityPlot(weight_single_list, delau_x_cumul, delau_ave_cumul, "Averaged delaunay triangulation segment length cumulative density", folder+"delau_seg_cumul_vs_rand.png")

        CumulDensityPlot(weight_single_list, area_x_cumul, area_ave_cumul, "Averaged voronoi domains area cumulative density", folder+"voro_area_cumul_vs_rand.png")

        CumulDensityPlot(weight_single_list, angle_x_cumul, angle_ave_cumul, "Averaged voronoi domains angles cumulative density", folder+"voro_angle_cumul_vs_rand.png")

        riPlot(weight_single_list, ri_ave, ri_std, "RI score depending on randomness weight for mosaic creation", folder+"ri_vs_random.png")

        CumulDensityPlot(weight_single_list, shortest_x_list, shortest_ave_list, "Averaged closest neighbour cumulative density", folder+"closest_neighbour_cumul_vs_rand.png")

        CumulDensityPlot(weight_single_list, dist_x_list, dist_ave_list, "Averaged neighbour distance cumulative density", folder+"neighbour_cumul_vs_rand.png")

    if stat_analysis:
    # methods sensitivity analysis
        significanceTable(folder, weight_single_list,\
         delauStats(weight_single_list, delau_ave_cumul),\
         voroStats(weight_single_list, area_ave_cumul, angle_ave_cumul),\
         riStats(mosaic_repetitions, weight_single_list, ri_list, shortest_ave_list, dist_ave_list))
    #TODO: false positives (specificity)

    if show_plots:
        plt.show()

    return 1

#--------------------------------------------------------------------------#
def mosaic_creation(output_folder, weight, n):
    if float(weight) > 1 or float(weight) < 0:
        SystemExit('random weigh should be between 0 and 1')
    cell_per_dim = 20
    cells_space = 25
    rand = float(weight) * 20

    for repetition in range(0, n):
        if os.path.isfile(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+str(repetition)+"_"+str(weight)+".txt"):
            os.remove(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+str(repetition)+"_"+str(weight)+".txt")
        output_file=open(output_folder+"mosaic_"+str(cell_per_dim*cell_per_dim)+"cells_"+str(repetition)+"_"+str(weight)+".txt", "a")
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
def delaunay(output_folder, positions_list, random_weight, plot=False):
    tri = sp.spatial.Delaunay(positions_list)#, qhull_options="QJ")

    if plot:
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

    return tri

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
def voronoi(output_folder, positions_list, random_weight, plot=False):
    voro = sp.spatial.Voronoi(positions_list)#, qhull_options="Qc") # Qc: sensitive parameter

    if plot:
        # voro.vertices: delaunays circles centers
        # voro.ridge_points: index of voro.vertices points couple forming lines
        # voro.ridge_vertices: index of voro.vertices points couple forming lines, including outside of region as -1
        # voro.regions: index of vertices points composing regions. closed region if no -1 index
        sp.spatial.voronoi_plot_2d(voro)

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

    return voro

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

    return round((np.average(shortest_dist_list)/np.std(shortest_dist_list)), 2), shortest_dist_list, dist_list

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
def sortData(label, a, b, c, d, e):
    # insertion sort
    i = 1
    while i < len(label):
        j = i
        while j > 0 and label[j-1] > label[j]:
            label[j], label[j-1] = label[j-1], label[j]
            a[j], a[j-1] = a[j-1], a[j]
            b[j], b[j-1] = b[j-1], b[j]
            c[j], c[j-1] = c[j-1], c[j]
            d[j], d[j-1] = d[j-1], d[j]
            e[j], e[j-1] = e[j-1], e[j]
            j = j-1
        i = i+1

    return label, a, b, c, d, e

#--------------------------------------------------------------------------#
def CumulDensityPlot(weight_single_list, x_cumul, ave_cumul, title, save_name):
    fig2 = plt.figure(); ax2=fig2.add_subplot(111)
    for i in range(0, len(weight_single_list)):
        ax2.plot(x_cumul[i], ave_cumul[i]/max(ave_cumul[i]), label=str(weight_single_list[i]))
    ax2.set_title(title)
    ax2.legend()
    ax2.get_figure().savefig(save_name)

#--------------------------------------------------------------------------#
def old_delauPlot(folder, random_weight_list, delau_list):
    # repetitions
    nb_of_repetitions = 1; i = 0
    while random_weight_list[i] == random_weight_list[i+1]:
        nb_of_repetitions += 1; i += 1
    weight_single_list = []
    for i in range(0, int(len(random_weight_list)/nb_of_repetitions)):
        weight_single_list.append(random_weight_list[i*nb_of_repetitions])

    list_seg_length = []
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
        list_seg_length.append(seg_length)

    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    fig2 = plt.figure(); ax2=fig2.add_subplot(111)

    ave_cumul = []; temps_ave_cumul = []; x_cumul = []; temps_x_cumul = []
    for index in range(0, len(list_seg_length)):
        n = ax1.hist(np.sort(list_seg_length[index])[:int(len(list_seg_length[index])-len(list_seg_length[index])*0.05)], bins=int(len(list_seg_length[index])/40))

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1]) # n[1][:len(n[1])-1]

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0
                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            ave_cumul.append(ave_cumul_point)
            x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1])

    for i in range(0, len(weight_single_list)):
        ax2.plot(x_cumul[i], ave_cumul[i]/max(ave_cumul[i]), label=str(weight_single_list[i]))
    ax2.set_title("Averaged delaunay triangulation segment length cumulative density")
    ax2.legend()
    ax2.get_figure().savefig(folder+"delau_seg_cumul_vs_rand.png")

#--------------------------------------------------------------------------#
def singleDelauPlot(folder, random_weight_list, delau_list):
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
def old_voroPlot(folder, random_weight_list, voro_list):
    # repetitions
    nb_of_repetitions = 1; i = 0
    while random_weight_list[i] == random_weight_list[i+1]:
        nb_of_repetitions += 1; i += 1
    weight_single_list = []
    for i in range(0, int(len(random_weight_list)/nb_of_repetitions)):
        weight_single_list.append(random_weight_list[i*nb_of_repetitions])

    list_areas_list = []; list_angles_list = []
    for voro_index in range(0, len(voro_list)):
        areas_list = []; angles_list = []
        voro = voro_list[voro_index]
        for region in voro.regions:
            # if region doesn't have point outside. ignore non closed domains (border domains), avoid border effect
            if len(region) > 0 and -1 not in region:
                domain_points = []
                for index in region:
                    domain_points.append(list(voro.vertices[index]))
                areas_list.append(polygoneArea(domain_points))
                angles_list.append(polygoneAngles(domain_points))
        list_areas_list.append(areas_list)
        list_angles_list.append(angles_list)

    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    fig2 = plt.figure(); ax2=fig2.add_subplot(111)
    ave_cumul = []; temps_ave_cumul = []; x_cumul = []; temps_x_cumul = []
    for index in range(0, len(list_areas_list)):
        n = ax1.hist(np.sort(list_areas_list[index])[:int(len(list_areas_list[index])-len(list_areas_list[index])*0.1)], bins=int(len(list_areas_list[index])/10), density=True)

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0

                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            ave_cumul.append(ave_cumul_point)
            x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

    for i in range(0, len(weight_single_list)):
        ax2.plot(x_cumul[i], ave_cumul[i]/max(ave_cumul[i]), label=str(weight_single_list[i]))
    ax2.set_title("Averaged voronoi domains area cumulative density")
    ax2.legend()
    ax2.get_figure().savefig(folder+"voro_area_cumul_vs_rand.png")

    fig3 = plt.figure(); ax3=fig3.add_subplot(111)
    fig4 = plt.figure(); ax4=fig4.add_subplot(111)
    ave_cumul = []; temps_ave_cumul = []; x_cumul = []; temps_x_cumul = []
    for index in range(0, len(list_angles_list)):
        angles = []
        for angle_sub_list in list_angles_list[index]:
            for angle in angle_sub_list:
                angles.append(np.degrees(angle))
        n = ax3.hist(angles, bins=int(len(angles)/100), density=True)

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0
                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            ave_cumul.append(ave_cumul_point)
            x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

    for i in range(0, len(weight_single_list)):
        ax4.plot(x_cumul[i], ave_cumul[i]/max(ave_cumul[i]), label=str(weight_single_list[i]))
    ax4.set_title("Averaged voronoi domains angles cumulative density")
    ax4.legend()
    ax4.get_figure().savefig(folder+"voro_angle_cumul_vs_rand.png")

#--------------------------------------------------------------------------#
def singleVoroPlot(folder, random_weight_list, voro_list):
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
def riPlot(weight_single_list, ri_ave, ri_std, title, save_name):
    plt.figure()
    plt.errorbar(weight_single_list, ri_ave, ri_std)
    plt.title(title)
    plt.savefig(save_name)

#--------------------------------------------------------------------------#
def riPlots(folder, random_weight_list, ri, short_dist_list, dist_list):
    # repetitions
    nb_of_repetitions = 1; i = 0
    while random_weight_list[i] == random_weight_list[i+1]:
        nb_of_repetitions += 1; i += 1

    # group ri values by random weight used
    ri_grouped = []
    for i in range(0, int(len(ri)/nb_of_repetitions)):
        ri_grouped.append(ri[i*nb_of_repetitions:i*nb_of_repetitions+nb_of_repetitions])

    # plot ri value vs random with error bar
    ri_ave = []; ri_std = []
    for group_vals in ri_grouped:
        ri_ave.append(np.average(group_vals))
        ri_std.append(np.std(group_vals))
    plt.figure()
    plt.errorbar(weight_single_list, ri_ave, ri_std)
    plt.title("RI measure depending on randomness weight for mosaic creation")
    plt.savefig(folder+"ri_vs_random.png")

    fig2 = plt.figure(); ax2=fig2.add_subplot(111)
    x_cumul, ave_cumul = getAveCumulative(short_dist_list, nb_of_repetitions)
    for i in range(0, len(weight_single_list)):
        ax2.plot(x_cumul[i], ave_cumul[i]/max(ave_cumul[i]), label=str(weight_single_list[i]))
    ax2.set_title("Averaged closest neighbour cumulative density")
    ax2.legend()
    ax2.get_figure().savefig(folder+"closest_neighbour_cumul_vs_rand.png")

    fig4 = plt.figure(); ax4=fig4.add_subplot(111)
    x_cumul, ave_cumul = getAveCumulative(dist_list, nb_of_repetitions)
    for i in range(0, len(weight_single_list)):
        ax4.plot(x_cumul[i], ave_cumul[i]/max(ave_cumul[i]), label=str(weight_single_list[i]))
    ax4.set_title("Averaged neighbour distance cumulative density")
    ax4.legend()
    ax4.get_figure().savefig(folder+"neighbour_cumul_vs_rand.png")

#--------------------------------------------------------------------------#
def singleRiPlots(folder, random_weight_list, ri, short_dist_list, dist_list):
    # ri value vs random
    plt.figure()
    plt.plot(random_weight_list, ri)
    plt.title("RI measure depending on randomness weight for mosaic creation")
    plt.savefig(folder+"ri_vs_random.png")

    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    fig2 = plt.figure(); ax2=fig2.add_subplot(111)
    fig3 = plt.figure(); ax3=fig3.add_subplot(111)
    fig4 = plt.figure(); ax4=fig4.add_subplot(111)

    for dist_index in range(0, len(short_dist_list)):
        # closest neighbour distance distribution
        n = ax1.hist(short_dist_list[dist_index], density=True)
        # closest neighbour distance cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
        ax2.plot(n[1][:len(n[1])-1], density/sum(n[0]), label=str(random_weight_list[dist_index]))

    for dist_index in range(0, len(dist_list)):
        # neighbours distance distribution
        n = ax3.hist(dist_list[dist_index], density=True)
        # neighbours distance culumative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))
        ax4.plot(n[1][:len(n[1])-1], density/sum(n[0]), label=str(random_weight_list[dist_index]))

    ax1.set_title("Closest neighbour distance distribution")
    ax1.get_figure().savefig(folder+"closest_neighbour_distrib_vs_rand.png")

    ax2.set_title("Closest neighbour cumulative density")
    ax2.legend()
    ax2.get_figure().savefig(folder+"closest_neighbour_cumul_vs_rand.png")

    ax3.set_title("Neighbour distance distribution")
    ax3.get_figure().savefig(folder+"neighbour_distrib_vs_rand.png")

    ax4.set_title("Neighbour distance cumulative density")
    ax4.legend()
    ax4.get_figure().savefig(folder+"neighbour_cumul_vs_rand.png")

#--------------------------------------------------------------------------#
def getAveCumulative(data, nb_of_repetitions):
    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    ave_cumul = []; temps_ave_cumul = []; x_cumul = []; temps_x_cumul = []
    for index in range(0, len(data)):
        n = ax1.hist(data[index], bins=50, density=True)
        plt.close(fig1)

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1]) # n[1][:len(n[1])-1]

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0
                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            ave_cumul.append(ave_cumul_point)
            x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1])

    return x_cumul, ave_cumul

#--------------------------------------------------------------------------#
def getWeightSingleList(random_weight_list, nb_of_repetitions):
    weight_single_list = []
    for i in range(0, int(len(random_weight_list)/nb_of_repetitions)):
        weight_single_list.append(random_weight_list[i*nb_of_repetitions])

    return weight_single_list

#--------------------------------------------------------------------------#
def delauAveCumul(weight_single_list, nb_of_repetitions, delau_list):
    list_seg_length = []
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
        list_seg_length.append(seg_length)

    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    ave_cumul = []; temps_ave_cumul = []; x_cumul = []; temps_x_cumul = []
    for index in range(0, len(list_seg_length)):
        n = ax1.hist(np.sort(list_seg_length[index])[:int(len(list_seg_length[index])-len(list_seg_length[index])*0.05)], bins=50, density=True)

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1]) # n[1][:len(n[1])-1]

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0

                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            ave_cumul.append(ave_cumul_point)
            x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1])

    return x_cumul, ave_cumul

#--------------------------------------------------------------------------#
def voroAveCumul(weight_single_list, nb_of_repetitions, voro_list):
    list_areas_list = []; list_angles_list = []
    for voro_index in range(0, len(voro_list)):
        areas_list = []; angles_list = []
        voro = voro_list[voro_index]
        for region in voro.regions:
            # if region doesn't have point outside. ignore non closed domains (border domains), avoid border effect
            if len(region) > 0 and -1 not in region:
                domain_points = []
                for index in region:
                    domain_points.append(list(voro.vertices[index]))
                areas_list.append(polygoneArea(domain_points))
                angles_list.append(polygoneAngles(domain_points))
        list_areas_list.append(areas_list)
        list_angles_list.append(angles_list)

    fig1 = plt.figure(); ax1=fig1.add_subplot(111)
    area_ave_cumul = []; temps_ave_cumul = []; area_x_cumul = []; temps_x_cumul = []
    for index in range(0, len(list_areas_list)):
        n = ax1.hist(np.sort(list_areas_list[index])[:int(len(list_areas_list[index])-len(list_areas_list[index])*0.1)], bins=50, density=True)
        # for nice figure, increase *0.1 to ~0.25
        plt.close(fig1)

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0

                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            area_ave_cumul.append(ave_cumul_point)
            area_x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

    fig3 = plt.figure(); ax3=fig3.add_subplot(111)
    angle_ave_cumul = []; temps_ave_cumul = []; angle_x_cumul = []; temps_x_cumul = []
    for index in range(0, len(list_angles_list)):
        angles = []
        for angle_sub_list in list_angles_list[index]:
            for angle in angle_sub_list:
                angles.append(np.degrees(angle))
        n = ax3.hist(angles, bins=25, density=True)
        plt.close(fig3)

        # cumulative density
        density = []
        for i in range(0, len(n[0])):
            density.append((n[0][i] + density[i-1]) if i > 0 else (n[0][i]))

        # if end of this random value group
        if index % nb_of_repetitions == nb_of_repetitions-1:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

            ave_cumul_point = []; ave_cumul_x = []
            for i in range(0, len(temps_ave_cumul[0])):
                temps_density = 0; temps_x = 0
                for repetition in range(0, len(temps_ave_cumul)):
                    temps_density += temps_ave_cumul[repetition][i]
                    temps_x = temps_x_cumul[repetition][i]
                # average cumulative density value for this point
                ave_cumul_point.append(temps_density/len(temps_ave_cumul))
                ave_cumul_x.append(temps_x/len(temps_x_cumul))

            # construct average density distribution for this random value
            angle_ave_cumul.append(ave_cumul_point)
            angle_x_cumul.append(ave_cumul_x)
            temps_ave_cumul = []
            temps_x_cumul = []
        else:
            temps_ave_cumul.append(density)
            temps_x_cumul.append(n[1][:len(n[1])-1])

    return area_x_cumul, area_ave_cumul, angle_x_cumul, angle_ave_cumul

#--------------------------------------------------------------------------#
def riAveCumul(weight_single_list, nb_of_repetitions, ri_list, short_dist_list, dist_list):
    ri_grouped = []
    for i in range(0, int(len(ri_list)/nb_of_repetitions)):
        ri_grouped.append(ri_list[i*nb_of_repetitions:i*nb_of_repetitions+nb_of_repetitions])
    ri_ave = []; ri_std = []
    for group_vals in ri_grouped:
        ri_ave.append(np.average(group_vals))
        ri_std.append(np.std(group_vals))

    short_x_cumul, short_ave_cumul = getAveCumulative(short_dist_list, nb_of_repetitions)

    dist_x_cumul, dist_ave_cumul = getAveCumulative(dist_list, nb_of_repetitions)

    return ri_ave, ri_std, short_x_cumul, short_ave_cumul, dist_x_cumul, dist_ave_cumul

#--------------------------------------------------------------------------#
def significanceTable(output_folder, weight_single_list, delau, voro, dists):

    columns = []; rows =["Delaunay segment length", "Voronoi areas", " Voronoi angles", "Regularity index", "Closest neighbour", "Neighbours distances"]
    for i in range(0, len(weight_single_list)-1):
        columns.append(str(weight_single_list[i])+" - "+str(weight_single_list[i+1]))

    cells_data = [delau, voro[0], voro[1], dists[0], dists[1], dists[2]]

    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off'); ax.axis('tight')
    # colWidths: list of float. decrease colWidths[0]
    ax.table(cellText=cells_data, cellLoc='center', rowLabels=rows, rowLoc='center', colLabels=columns, loc='center')
    plt.subplots_adjust(left=0.4, right=0.94)
    plt.savefig(output_folder+"methods_significance_table.png", dpi = 256)

#--------------------------------------------------------------------------#
def delauStats(weight_single_list, delau_cumul):
    # significance
    seg_lenght_signi = []
    for i in range(0, len(weight_single_list)-1):
        seg_lenght_signi.append(signi_star(sp.stats.ks_2samp(delau_cumul[i], delau_cumul[i+1])))

    return seg_lenght_signi

#--------------------------------------------------------------------------#
def voroStats(weight_single_list, area_cumul, angle_cumul):
    # significance
    area_signi = []; angle_signi = []
    for i in range(0, len(weight_single_list)-1):
        area_signi.append(signi_star(sp.stats.ks_2samp(area_cumul[i], area_cumul[i+1])))
        angle_signi.append(signi_star(sp.stats.ks_2samp(angle_cumul[i], angle_cumul[i+1])))

    return area_signi, angle_signi

#--------------------------------------------------------------------------#
def riStats(nb_of_repetitions, random_weight_list, ri_list, short_dist_list, dist_list):
    ri_signi = []; ri_grouped = []; short_dist_signi = []; dist_signi = []
    for i in range(0, int(len(ri_list)/nb_of_repetitions)):
        ri_grouped.append(ri_list[i*nb_of_repetitions:i*nb_of_repetitions+nb_of_repetitions])
    # significance
    for i in range(0, len(random_weight_list)-1):
        ri_signi.append(signi_star(sp.stats.ttest_ind(ri_grouped[i], ri_grouped[i+1])))
        short_dist_signi.append(signi_star(sp.stats.ks_2samp(short_dist_list[i], short_dist_list[i+1])))
        dist_signi.append(signi_star(sp.stats.ks_2samp(dist_list[i][:len(dist_list)-1], dist_list[i+1][:len(dist_list)-1])))

    return ri_signi, short_dist_signi, dist_signi

#--------------------------------------------------------------------------#
def signi_star(stat):
    if stat[1] < 0.05:
        if stat[1] < 0.001:
            return "***"
        elif stat[1] < 0.01:
            return "**"
        else:
            return "*"
    else:
        return "."

#--------------------------------------------------------------------------#
# check number of arguments
if len(sys.argv)==2:
    if main(sys.argv[1]):
        print("done")
    else:
        print("error during execution")
else:
    raise SystemExit('Error: need 1 arg: [input-output folder]')
