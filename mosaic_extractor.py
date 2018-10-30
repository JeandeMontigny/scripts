#!/usr/bin/env python3
import sys, os, re

#--------------------------------------------------------------------------#
def main(input_file, folderName):
    RGC_dico={}
    if not os.path.exists(folderName):
        os.makedirs(folderName)
    else:
        for fileName in os.listdir(folderName):
            os.remove(folderName+"/"+fileName)

    with open(input_file, 'r') as fichier:
        for line in fichier:
            #TODO: just use the bias index (towards 1 means ON response, towards -1 OFF and around 0 ON-OFF. We set normally the boundaries to +/- 0.33)
            m=re.search( r'(.+) .+ .+ (.+) (.+) .+ .+ .+ .+ (.+)', line, re.M|re.I) # name FF_RQ channel x y ... Bias_Index
            if m and str(m.group(1))!="RGC_type":
                RGC_type=str(m.group(1))
                x_pos=round(float(m.group(2)), 2)
                y_pos=round(float(m.group(3)), 2)
                insertIntoDico(RGC_dico, RGC_type, x_pos, y_pos)

    for each_type in RGC_dico.keys():
        writeToFile(folderName, each_type, RGC_dico[each_type])

#--------------------------------------------------------------------------#
def writeToFile(folder, cell_type, positions):
    sortie=open(folder+"/"+cell_type+".txt", "a")
    for position in positions:
        sortie.write(str(position).replace("(", "").replace(")", "").replace(",", "")+"\n")

#--------------------------------------------------------------------------#
def insertIntoDico(dico, RGC_type, x, y):
    if not RGC_type in dico:
        dico[RGC_type] = [(x,y)]
    else:
        dico[RGC_type].append((x,y))

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
    m=re.search( r'.+/(.+)\.txt', sys.argv[1], re.M|re.I)
    if m:
        main(sys.argv[1], m.group(1))
    else:
        exit("file name parsing error. Should be a .txt file")
    print("done")
else:
    exit("arg error - need 1 arg: [gerrit file]")

