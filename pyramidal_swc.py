#!/usr/bin/env python3
import sys, os, re
import numpy as np

#--------------------------------------------------------------------------#
def main(folder):
    for name in os.listdir(folder):
        if name.endswith(".swc"):
            print(read(folder+"/"+name))

#--------------------------------------------------------------------------#
def read(file_name):
    swc_file = open(file_name, "r")

    point_coord = [0, 0, 0]; previous_point_coord = [0, 0, 0]; previous_branching_point_coord = [0, 0, 0]
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
            if int(m.group(2)) == 3 or int(m.group(2)) == 5 or int(m.group(2)) == 6:
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
