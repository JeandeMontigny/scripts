#!/usr/bin/env python3
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(folder):
    temps_branching = []; temps_length = []
    for name in os.listdir(folder):
        if name.endswith(".swc"):
            results = read(folder+"/"+name)
            temps_branching.append(results[0])
            temps_length.append(results[1])

    branching_ave = round(np.average(temps_branching), 2); branching_std = round(np.std(temps_branching), 2)
    length_ave = round(np.average(temps_length), 2); length_std = round(np.std(temps_length), 2)
    print("average of", branching_ave, "branching point, with std of", branching_std, "\naverage dendritic length of", length_ave, " with std of", length_std)

    plot(branching_ave, branching_std, length_ave, length_std)

#--------------------------------------------------------------------------#
def read(file_name):
    swc_file = open(file_name, "r")

    previous_point_coord = [0, 0, 0]; previous_branching_point_coord = [0, 0, 0]
    branching_nb = 0; length = 0

    for line in swc_file:
        m=re.search( r' ?([0-9]+) ([0-9]+) (.+) (.+) (.+) (.+) ([-0-9]+).?', line, re.M|re.I)
        if m and (line[0] != "#" or line[0] != "\n"):
            point_coord = [float(m.group(3)), float(m.group(4)), float(m.group(5))]

            # if the point is the soma
            if int(m.group(2) == -1):
                # previous point become coordinate of the soma
                previous_point_coord = point_coord

            # if the point is a segment, branching point or a terminal point
            if int(m.group(2)) == 3 or int(m.group(2)) == 4 or int(m.group(2)) == 5 or int(m.group(2)) == 6:
                # add this segment length to total length
                length += dist(point_coord, previous_point_coord)
                previous_point_coord = point_coord
                # if it is a branching point
                if int(m.group(2)) == 5:
                    branching_nb += 1
                    # previous branching point become this point
                    previous_branching_point_coord = point_coord
                # if it is a terminal point
                if int(m.group(2)) == 6:
                    # then the previous point needs to be the previous branching point
                    previous_point_coord = previous_branching_point_coord

    return branching_nb, length

#--------------------------------------------------------------------------#
def dist(current, previous):
    return round(np.sqrt(np.square(current[0]-previous[0])+np.square(current[1]-previous[1])+np.square(current[2]-previous[2])), 3)

#--------------------------------------------------------------------------#
def plot(branching_ave, branching_std, length_ave, length_std):
    figure([branching_ave, branching_std], [37.19, 10.04], "Average number of branching points")
    figure([length_ave, length_std], [5942.89, 1590.46], "Average length of dendritic tree")

    plt.show()

#--------------------------------------------------------------------------#
def figure(simu, real, title):
    plt.figure()
    y_pos = np.arange(2)
    plt.bar(y_pos, [simu[0], real[0]], 0.75, color=['dimgrey', 'darkgray'], yerr=[simu[1], real[1]], error_kw=dict(ecolor='black'), align='center')
    plt.xticks(y_pos, ["simulated neurons", "real neurons"])
    plt.title(title)

#--------------------------------------------------------------------------#
if len(sys.argv) != 3:
    exit("arg error - need 2 arg: [folder swc files] [rewrite,-]")
else:
    if sys.argv[2]=="rewrite":
        import correction
        correction.main(sys.argv[1])
        arg=sys.argv[1]+"/rewrite"
        main(arg)
    else:
        main(sys.argv[1])
    print("done")
