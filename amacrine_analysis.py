#!/usr/bin/env python3
import sys, os, re
import numpy as np

#--------------------------------------------------------------------------#
def main(folder):
    # extract info from folder files
    folder_data = getFolderData(folder)
    # create relevant pairs of GCL/INL files
    files_couple = createPairs(folder_data)
    # process each GCL/INL pair
    results = analyse(files_couple, folder)

    # compute means and std
    tab_gcl_ri=[]; tab_inl_ri=[]; tab_gcl_inl_ratio=[]; tab_both_ri=[]; tab_exclusionFactor=[]
    for gcl_ri, inl_ri, gcl_inl_ratio, both_ri, exclusionFactor in results:
        tab_gcl_ri.append(gcl_ri); tab_inl_ri.append(inl_ri); tab_gcl_inl_ratio.append(gcl_inl_ratio); tab_both_ri.append(both_ri); tab_exclusionFactor.append(exclusionFactor)

    print("gcl ri:", np.average(tab_gcl_ri), np.std(tab_gcl_ri))
    print("inl ri:", np.average(tab_inl_ri), np.std(tab_inl_ri))
    print("gcl/inl ratio:", np.average(tab_gcl_inl_ratio), np.std(tab_gcl_inl_ratio))
    print("both ri:", np.average(tab_both_ri), np.std(tab_both_ri))
    print("exclusion factor:", np.average(tab_exclusionFactor), np.std(tab_exclusionFactor))

#--------------------------------------------------------------------------#
def analyse(files_couple, folder):
    results = []
    for couple in files_couple:
        results.append(processCouple(folder+"/"+couple[0], folder+"/"+couple[1]))
    return results

#--------------------------------------------------------------------------#
def processCouple(gcl_coord_file, inl_coord_file):

    # extract cells position
    gcl_cells_position = getCellsPosition(gcl_coord_file)
    inl_cells_position = getCellsPosition(inl_coord_file)
    # independant mosaics
    gcl_ri = getRi(getShortestDistList(gcl_cells_position))
    inl_ri = getRi(getShortestDistList(inl_cells_position))
    # cell layers ratio
    gcl_inl_ratio = len(gcl_cells_position)/len(inl_cells_position)
    # combined mosaic
    both_ri = getRi(getShortestDistList(mergeTab(gcl_cells_position, inl_cells_position)))
    # exclusion
    exclusionFactor = getExclusionFactor(gcl_cells_position, inl_cells_position)

    return [gcl_ri, inl_ri, gcl_inl_ratio, both_ri, exclusionFactor]

#--------------------------------------------------------------------------#
def getExclusionFactor(tab_1, tab_2):
    """ return a double between [0,1]: 1 if mosaics completely excluse one each other, 0 if they are totaly superimposed """
    redondancy = 0
    r = 10
    for cell_1 in tab_1:
        for cell_2 in tab_2:
            if (np.square(cell_1[0] - cell_2[0]) + np.square(cell_1[1] - cell_2[1]) <= np.square(r)):
                redondancy+=1
                # exit loop to avoid multiple positif for one cell
                break
    return round(1 - redondancy/len(tab_1+tab_2)*2, 2)

#--------------------------------------------------------------------------#
def mergeTab(tab_1, tab_2):
    return (tab_1 + tab_2)

#--------------------------------------------------------------------------#
def getRi(shortest_dist_list):
    return round((np.average(shortest_dist_list)/np.std(shortest_dist_list)), 2)

#--------------------------------------------------------------------------#
def getCellsPosition(coord_file):
    cells_position=[]
    fichier=open(coord_file, "r")
    file_lines=fichier.readlines()
    for line in file_lines:
        # ...,x,y
        m=re.search( r'[0-9]+,[0-9.]+,[0-9.]+,[0-9.]+,[0-9.]+,([0-9.]+),([0-9.]+)', line, re.M|re.I)
        if m :
            cells_position.append([float(m.group(1)), float(m.group(2))])
    return cells_position

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
def createPairs(folder_data):
    files_couple = []
    for data1 in folder_data:
        for data2 in folder_data:
            if (data1[1] == data2[1]) and (data1[2] != data2[2]):
                if data1[2] == 'GCL':
                    files_couple.append([data1[3], data2[3]])
                else:
                    files_couple.append([data2[3], data1[3]])
                folder_data.remove(data2)
                break
    return files_couple
#--------------------------------------------------------------------------#
def getFolderData(folder):
    folder_data = []
    for file in os.listdir(folder):
        m=re.search( r'p(.)_peri_(.)_(...)\.csv', file, re.M|re.I)
        if m:
            folder_data.append([m.group(1), m.group(2), m.group(3), file])
    return folder_data

#--------------------------------------------------------------------------#
if len(sys.argv)!=2:
    exit("arg error - need 1 arg: [1dev_day_amacrine_position_folder]")
else:
    main(sys.argv[1])

