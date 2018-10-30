#!/usr/bin/env python3
import sys, os, re

# corrects swc file export, adding branching point label to branching points
# corrected files are placed in the rewrite folder, inside the folder containing original swc files
#--------------------------------------------------------------------------#
def main(path_dir):
    os.chdir(path_dir)
    if not os.path.isdir("./rewrite"):
        os.mkdir("./rewrite")
    print("correcting files..")
    number_of_file = len([name for name in os.listdir('.') if os.path.isfile(name)])
    nb_file = 0
    for file in os.listdir("./"):
        if file.endswith(".swc"):
            nb_file+=1
            # dynamic print of corrected file number
            sys.stdout.write("\r"+str(nb_file)+"/"+str(number_of_file))
            sys.stdout.flush()
            
            lecture(file)

    print("done")
#--------------------------------------------------------------------------#
def lecture(file_name):
    fichier=open("./"+file_name, "r")
    file_lines=fichier.readlines()
    branching_index=[]; prev_parent=0; terminasion_index=[]; prev_id=0; line_non_cell=0
    for line in file_lines:
            # (label) (parent label) 
        m = re.search( r' ?([0-9]+) [0-9]+ .+ .+ .+ .+ ([-0-9]+).?', line, re.M|re.I)
        # if line is not a comment or empty
        if m and (line[0]!="#" or line[0]!="\n"):
            # if the parent of this point is a branching point
            if (int(m.group(2))<prev_parent) and (int(m.group(2))>1):
                branching_index.append(int(m.group(2)))
                terminasion_index.append(prev_id)
            prev_parent=int(m.group(2))
            prev_id=int(m.group(1))
        elif line[0]=="#" or line[0]=="\n":
            line_non_cell+=1

    terminasion_index.append(len(file_lines)-line_non_cell)

    write(file_name, branching_index, terminasion_index)

#--------------------------------------------------------------------------#
def write(file_name, branching_index, terminasion_index):
    output=open("rewrite/re_"+file_name, "w")
    fichier=open("./"+file_name, "r")
    file_lines=fichier.readlines()
    for  line in file_lines:
        m = re.search( r' ?([0-9]+) ([0-9]+) (.+) (.+) (.+) (.+) ([-0-9]+).?', line, re.M|re.I)
        if m and (line[0]!="#" or line[0]!="\n"):
            if int(m.group(1)) in branching_index and int(m.group(2))==3:
                output.write(m.group(1)+" 5 "+m.group(3)+" "+m.group(4)+" "+m.group(5)+" "+m.group(6)+" "+m.group(7)+"\n")
            elif int(m.group(1)) in terminasion_index and int(m.group(2))==3:
                output.write(m.group(1)+" 6 "+m.group(3)+" "+m.group(4)+" "+m.group(5)+" "+m.group(6)+" "+m.group(7)+"\n")
            else:
                output.write(m.group(1)+" "+m.group(2)+" "+m.group(3)+" "+m.group(4)+" "+m.group(5)+" "+m.group(6)+" "+m.group(7)+"\n")

