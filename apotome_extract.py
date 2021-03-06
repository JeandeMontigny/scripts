#!/usr/bin/env python3
import sys, os, re, warnings
import numpy as np
import matplotlib.pyplot as plt

# TODO: improve: SAC_INL are parsed twice
#--------------------------------------------------------------------------#
def main(folder, output_path):
    """ main: create output files containing number of cells per image and corresponding RI values.
    Walk through folders to find coordinate files
    @params:
        folder containing sub-folders and coordinate files"""
    # create and open output file
    output_RGC = open(output_path+"ri_measures_output_RGC.out", "a")
    output_SAC_GCL = open(output_path+"ri_measures_output_SAC_GCL.out", "a")
    output_SAC_INL = open(output_path+"ri_measures_output_SAC_INL.out", "a")
    # tab for figures
    figures_tab = []
    # tab for exclusion factor figure
    exclusion_tab = []
    # prepare progression bar
    nb_of_files = getNbOfFiles(folder)
    current_file_nb = 0
    # walk across folders
    for dirpath, dirnames, filenames in os.walk(folder):
        # for each files in directory
        for name in filenames:
            # if file is a .txt
            if name.endswith(".txt"):
                # print progress bar
                current_file_nb += 1
                printProgressBar(current_file_nb, nb_of_files)
                # call read function to get the number of cells and ri value for this file
                values = read(dirpath, name)
                nb_of_cells = values[0]
                density = (nb_of_cells / (335.4*335.4)) * 10e5
                ri = values[1]
                # create string for output - format: dev_day, nb_of_cells, ri
                m = re.search(r'.+/P([0-9]+)', dirpath, re.M|re.I)
                string_output = str(m.group(1))+" "+str(density)+" "+str(ri)+"\n"
                tab_cell = [m.group(1), density, ri]
                if len(values) == 3:
                    exclusion = values[2]
                    exclusion_tab.append([int(m.group(1)), float(exclusion)])
                # write output into corresponding output file
                if name.endswith("RGC.txt"):
                    tab_cell.append("RGC"); figures_tab.append(tab_cell)
                    output_RGC.write(string_output)
                if name.endswith("GCL.txt"):
                    tab_cell.append("SAC_GCL"); figures_tab.append(tab_cell)
                    output_SAC_GCL.write(string_output)
                if name.endswith("INL.txt"):
                    tab_cell.append("SAC_INL"); figures_tab.append(tab_cell)
                    output_SAC_INL.write(string_output)

    # close output files
    output_RGC.close(); output_SAC_GCL.close(); output_SAC_INL.close()
    # create figures
    figureExclusion(exclusion_tab)
    figures(figures_tab)
    # plot figures
    plt.show()

#--------------------------------------------------------------------------#
def getNbOfFiles(folder):
    """ Return the number of files to process, including files in sub-folders.
    @params:
        folder containing sub-folders and files """
    nb_of_files = 0
    for dirpath, dirnames, filenames in os.walk(folder):
        for name in filenames:
            if name.endswith(".txt"):
                nb_of_files += 1

    return nb_of_files

#--------------------------------------------------------------------------#
def read(folder, name):
    """ Return number of cells and rounded RI value.
    To do that, extract x,y coordinate from ImageJ cell position extraction then compute the Regularity Index (RI) corresponding to this cell population.
    @params:
        folder containing coordinate files
        coordinate file with cells info, sorted as 'index Area Mean Min Max X Y' """
    cells_position = []
    nb_of_cells = 0
    coord_file = folder+"/"+name
    fichier = open(coord_file, "r")
    file_lines = fichier.readlines()
    # first, extract X Y position for each cell
    # for each line of this file
    for line in file_lines:
        # if line doesn't start with space character, ie, if it's not the file legend
        if line[0] != " ":
            # file strucure: index Area Mean Min Max X Y
            m = re.search( r'[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t([0-9.]+)\t([0-9.]+)', line, re.M|re.I)
            if m :
                cells_position.append([round(float(m.group(1)), 2), round(float(m.group(2)), 2)])
                nb_of_cells += 1
    # check if data has been extracted 
    if len(cells_position) == 0 or nb_of_cells == 0:
        warnings.warn("error after parsing file "+coord_file+" - check file format")
        return 0, 0
    # then call getShortestDistList function for this whole file
    shortest_dist_list = getShortestDistList(cells_position)
    # calculate ri, which is: mean of distances / std of distances
    ri = np.average(shortest_dist_list) / np.std(shortest_dist_list)
    # if file is gcl amacrine cells, we have to process the pair GCL/INL
    if name.endswith("GCL.txt"):
        # get corresponding inl pair cells position
        inl_cells_position = getCellsPosition(getPair(folder, name))
        # get GCL/INL pair exclusion factor
        exclusion_factor = getExclusionFactor(cells_position, inl_cells_position)

        return nb_of_cells, round(ri, 2), exclusion_factor

    return nb_of_cells, round(ri, 2)

#--------------------------------------------------------------------------#
def getShortestDistList(coord_list):
    """ Return a list containing shortest neighbour distance for each cell
    @param:
        tab of cells coordinates """
    shortest_dist_list = []
    # for each cell coordinate in this list
    for i in range(0, len(coord_list)):
        distance_list = []
        cell_coord = coord_list[i]
        # for each other cell in this list
        for j in range(0, len(coord_list)):
            other_cell_coord = coord_list[j]
            # calculate Euclidean distance between these two cells
            tempsDistance = np.sqrt(np.square(cell_coord[0] - other_cell_coord[0]) + np.square(cell_coord[1] - other_cell_coord[1]))
            # if cell is not itself
            if tempsDistance != 0:
                distance_list.append(tempsDistance)
        # add shortest distance
        shortest_dist_list.append(min(distance_list))

    return shortest_dist_list

#--------------------------------------------------------------------------#
def getPair(folder_data, name):
    """ Return SAC_INL file path corresponding to the SAC_GCL pair file 
    @params:
        SAC_GCL folder path
        SAC_GCL file name """
    files_couple = []
    pair_name = name[:len(name)-7]+"INL.txt"
    for file_name in os.listdir(folder_data):
        if file_name == pair_name:
            return folder_data+"/"+file_name

#--------------------------------------------------------------------------#
def getCellsPosition(coord_file):
    """ Return a list containing positions of each cell of a coordinate file
    @params:
        coordinate file """
    cells_position = []
    fichier = open(coord_file, "r")
    file_lines = fichier.readlines()
    for line in file_lines:
        # ...,x,y
        m = re.search( r'[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t[0-9.]+\t([0-9.]+)\t([0-9.]+)', line, re.M|re.I)
        if m :
            cells_position.append([round(float(m.group(1)), 2), round(float(m.group(2)), 2)])

    return cells_position

#--------------------------------------------------------------------------#
def getExclusionFactor(tab_1, tab_2):
    """ Return a double between [0,1]: 1 if mosaics completely excluse one each other, 0 if they are totaly superimposed
    @params:
        [x, y] coordinate list 1
        [x, y] coordinate list 2 """
    redondancy = 0
    r = 10
    for cell_1 in tab_1:
        for cell_2 in tab_2:
            if (np.square(cell_1[0] - cell_2[0]) + np.square(cell_1[1] - cell_2[1]) <= np.square(r)):
                redondancy += 1
                # exit loop to avoid multiple positif for one cell
                break

    return round(1 - redondancy/len(tab_1+tab_2)*2, 2)

#--------------------------------------------------------------------------#
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 0, length = 100, fill = "#"):
    """ Create progress bar
    @params:
        iteration   - Required: current iteration (Int)
        total       - Required: total iterations (Int)
        decimals    - Optional: positive number of decimals in percent complete (Int)
        length      - Optional: character length of bar (Int)
        fill        - Optional: bar fill character (Str)"""
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s [%s] %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print new line on complete
    if iteration == total: 
        print("\n")

#--------------------------------------------------------------------------#
def figures(figures_tab):
    """ Create figures for cell density and RI evolution through development
    @params:
        list containing all info for figure creation. format: dev_day, nb_of_cells, ri, cell_type """

    # get tab for figure creation
    data = createFiguresTabs(figures_tab)
    dev_days = data[0]
    rgc_nb = data[1]; rgc_nb_std = data[2]; rgc_ri=data[3]; rgc_ri_std = data[4]
    sac_gcl_nb = data[5]; sac_gcl_nb_std = data[6]; sac_gcl_ri=data[7]; sac_gcl_ri_std = data[8]
    sac_inl_nb = data[9]; sac_inl_nb_std = data[10]; sac_inl_ri=data[11]; sac_inl_ri_std = data[12]

    # RGC figure
    createFigure(dev_days, rgc_nb, rgc_nb_std, "Postnatal day", "Cell density", "RGC density during development")
    createFigure(dev_days, rgc_ri, rgc_ri_std, "Postnatal day", "RI", "RGC Regularity Index during development")
    createFigureRiEvolution(rgc_nb, rgc_ri)

    # SAC GCL figure
    createFigure(dev_days, sac_gcl_nb, sac_gcl_nb_std, "Postnatal day", "Cell density", "GCL SAC density during development")
    createFigure(dev_days, sac_gcl_ri, sac_gcl_ri_std, "Postnatal day", "RI", "GCL SAC Regularity Index during development")
    # SAC INL figure
    createFigure(dev_days, sac_inl_nb, sac_inl_nb_std, "Postnatal day", "Cell density", "INL SAC density during development")
    createFigure(dev_days, sac_inl_ri, sac_inl_ri_std, "Postnatal day", "RI", "INL SAC Regularity Index during development")
    
#--------------------------------------------------------------------------#
# TODO: improve: do everything in one loop
def createFiguresTabs(figures_tab):
    """ Create a dictionary to store and sort data for figures. Return sorted data for figure creation 
    @params:
        list containing non sorted data obtained from the main function """
    temps_dict = {}

    for file_info in figures_tab:
        day = file_info[0]; nb = file_info[1]; ri = file_info[2]; c_type = file_info[3]
        # if this developmental day is not in the dictionnary
        if not day in temps_dict:
            # add this dev day as a dictionnary key
            temps_dict.setdefault(day, [])
            # add data to this key
            temps_dict[day].append([nb, ri, c_type])
        else:
            # add data to date key
            temps_dict[day].append([nb, ri, c_type])

    return sortData(temps_dict)

#--------------------------------------------------------------------------#
# TODO: improve
def sortData(temps_dict):
    """ Return sorted data for figure creation
    @params:
        dictionnary containing non sorted data for figure creation. day as key, cell nb and ri as values """
    dev_days = []
    rgc_nb = []; rgc_nb_std = []; rgc_ri = []; rgc_ri_std = []
    sac_gcl_nb = []; sac_gcl_nb_std = []; sac_gcl_ri = []; sac_gcl_ri_std = []
    sac_inl_nb = []; sac_inl_nb_std = []; sac_inl_ri = []; sac_inl_ri_std = []

    for each_day in temps_dict.keys():
        dev_days.append(int(each_day))
        data=temps_dict[each_day]

        day_rgc_nb = []; day_rgc_ri = []
        day_gcl_nb = []; day_gcl_ri = []
        day_inl_nb = []; day_inl_ri = []

        for this_data in data:
            if this_data[2] == 'RGC':
                day_rgc_nb.append(this_data[0]); day_rgc_ri.append(this_data[1])
            if this_data[2] == 'SAC_GCL':
                day_gcl_nb.append(this_data[0]); day_gcl_ri.append(this_data[1])
            if this_data[2] == 'SAC_INL':
                day_inl_nb.append(this_data[0]); day_inl_ri.append(this_data[1])

        if len(day_rgc_nb) > 0 and len(day_rgc_ri) > 0:
            rgc_nb.append(round(np.average(day_rgc_nb), 2)); rgc_nb_std.append(round(np.std(day_rgc_nb), 2))
            rgc_ri.append(round(np.average(day_rgc_ri), 2)); rgc_ri_std.append(round(np.std(day_rgc_ri), 2))
        else:
            rgc_nb.append(0); rgc_nb_std.append(0)
            rgc_ri.append(0); rgc_ri_std.append(0)

        if len(day_gcl_nb) > 0 and len(day_gcl_ri) > 0:
            sac_gcl_nb.append(round(np.average(day_gcl_nb), 2)); sac_gcl_nb_std.append(round(np.std(day_gcl_nb), 2))
            sac_gcl_ri.append(round(np.average(day_gcl_ri), 2)); sac_gcl_ri_std.append(round(np.std(day_gcl_ri), 2))
        else:
            sac_gcl_nb.append(0); sac_gcl_nb_std.append(0)
            sac_gcl_ri.append(0); sac_gcl_ri_std.append(0)

        if len(day_inl_nb) > 0 and len(day_inl_ri) > 0:
            sac_inl_nb.append(round(np.average(day_inl_nb), 2)); sac_inl_nb_std.append(round(np.std(day_inl_nb), 2))
            sac_inl_ri.append(round(np.average(day_inl_ri), 2)); sac_inl_ri_std.append(round(np.std(day_inl_ri), 2))
        else:
            sac_inl_nb.append(0); sac_inl_nb_std.append(0)
            sac_inl_ri.append(0); sac_inl_ri_std.append(0)

    # TODO: create generic function for sorting
    # insertion sort: set correct order for direct figure creation
    for i in range(1, len(dev_days)):
        day_temp = dev_days[i]
        rgc_nb_temps = rgc_nb[i]; rgc_nb_std_temps = rgc_nb_std[i]; rgc_ri_temps = rgc_ri[i]; rgc_ri_std_temps = rgc_ri_std[i]
        sac_gcl_nb_temps = sac_gcl_nb[i]; sac_gcl_nb_std_temps = sac_gcl_nb_std[i]; sac_gcl_ri_temps = sac_gcl_ri[i]; sac_gcl_ri_std_temps = sac_gcl_ri_std[i]
        sac_inl_nb_temps = sac_inl_nb[i]; sac_inl_nb_std_temps = sac_inl_nb_std[i]; sac_inl_ri_temps = sac_inl_ri[i]; sac_inl_ri_std_temps = sac_inl_ri_std[i]

        j = i-1
        # search index from 0 to i-1 where to insert value[i]=temp
        while j >= 0 and dev_days[j] > day_temp:
            dev_days[j+1] = dev_days[j]
            rgc_nb[j+1] = rgc_nb[j]; rgc_nb_std[j+1] = rgc_nb_std[j]; rgc_ri[j+1] = rgc_ri[j]; rgc_ri_std[j+1] = rgc_ri_std[j]
            sac_gcl_nb[j+1] = sac_gcl_nb[j]; sac_gcl_nb_std[j+1] = sac_gcl_nb_std[j]; sac_gcl_ri[j+1] = sac_gcl_ri[j]; sac_gcl_ri_std[j+1] = sac_gcl_ri_std[j]
            sac_inl_nb[j+1] = sac_inl_nb[j]; sac_inl_nb_std[j+1] = sac_inl_nb_std[j]; sac_inl_ri[j+1] = sac_inl_ri[j]; sac_inl_ri_std[j+1] = sac_inl_ri_std[j]
            j -= 1
        # insert value[i] to its position
        dev_days[j+1] = day_temp
        rgc_nb[j+1] = rgc_nb_temps; rgc_nb_std[j+1] = rgc_nb_std_temps; rgc_ri[j+1] = rgc_ri_temps; rgc_ri_std[j+1] = rgc_ri_std_temps
        sac_gcl_nb[j+1] = sac_gcl_nb_temps; sac_gcl_nb_std[j+1] = sac_gcl_nb_std_temps; sac_gcl_ri[j+1] = sac_gcl_ri_temps; sac_gcl_ri_std[j+1] = sac_gcl_ri_std_temps
        sac_inl_nb[j+1] = sac_inl_nb_temps; sac_inl_nb_std[j+1] = sac_inl_nb_std_temps; sac_inl_ri[j+1] = sac_inl_ri_temps; sac_inl_ri_std[j+1] = sac_inl_ri_std_temps

    return dev_days, rgc_nb, rgc_nb_std, rgc_ri, rgc_ri_std, sac_gcl_nb, sac_gcl_nb_std, sac_gcl_ri, sac_gcl_ri_std, sac_inl_nb, sac_inl_nb_std, sac_inl_ri, sac_inl_ri_std

#--------------------------------------------------------------------------#
def createFigure(x, y, std, x_label, y_label, title):
    """ Create error bar figure
    @params:
        x data list
        y data list
        standard deviation list
        label for x axis (Str)
        label for y axis (Str)
        figure title (Str) """
    plt.figure()

    plt.errorbar(x, y, std)
    plt.xlim(1.5, 10.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

#--------------------------------------------------------------------------#
def createFigureRiEvolution(rgc_nb, rgc_ri):
    """ Create figure of ri evolution vs death rate """
    ri_diff = []; death_rate = []
    for i in range(0, len(rgc_nb)-4):
        death_rate.append( (12000 - rgc_nb[i])/12000*100 )
        if i > 0:
            ri_diff.append(rgc_ri[i]-rgc_ri[i-1])

    plt.figure()
    plt.plot(death_rate[0:len(death_rate)-1], ri_diff)
    plt.axhline(y=0, color='gray', linestyle='--')
    plt.xlabel("death rate")
    plt.ylabel("RI evolution")

#--------------------------------------------------------------------------#
def figureExclusion(tab):
    """ Prepare data for GCL/INL exclusion score figure and call figure creator
    @params:
        list from the main function containing each developmental day and exclusion score [ [day, score], [day, score] ]"""
    # get dev_days, score and score_std for exclusion figure
    data = sort(getFigureExclusionTab(tab))
    # create figure
    createFigure(data[0], data[1], data[2], "Postnatal day", "Exclusion score", "Exclusion score of GCL and INL mosaics during development")

#--------------------------------------------------------------------------#
def getFigureExclusionTab(tab):
    """Return (non sorted) data in list format for exclusion figure creation
    @params:
        list from figureExclusion function containing each developmental day and exclusion score [ [day, score], [day, score] ] """
    dev_days = []; score = []; score_std = []; temps = []
    for cell in tab:
        # if this developmental day is not in the list dev_days
        if not cell[0] in dev_days:
            # add it
            dev_days.append(cell[0])
            # if temps is not empty, ie if this is not the first time we enter the loop and so not the first developmental day we add to dev_days
            if not len(temps) == 0:
                # compute previous developmental day average and std score and add them to corresponding lists
                score.append(np.average(temps))
                score_std.append(np.std(temps))
            # empty temps list as this is a new developmental day
            temps = []
            # add this exclusion score to temps
            temps.append(cell[1])
        # if this developmental day is already in the list
        else:
            # just add this exclusion score to temps
            temps.append(cell[1])
    # add average and std score of the last developmental day to lists
    score.append(np.average(temps))
    score_std.append(np.std(temps))

    return dev_days, score, score_std

#--------------------------------------------------------------------------#
# TODO: create generic function for sorting
def sort(data):
    """ Return sorted data depending on development days list
    @params:
        list of non sorted data """
    dev_days = data[0]; score = data[1]; score_std = data[2]
    # insertion sort: set correct order for direct figure creation
    for i in range(1, len(dev_days)):
        day_temp = dev_days[i]; score_temps = score[i]; score_std_temps = score_std[i]
        j = i-1
        # search index from 0 to i-1 where to insert value[i]=temp
        while j >= 0 and dev_days[j] > day_temp:
            dev_days[j+1] = dev_days[j]; score[j+1] = score[j]; score_std[j+1] = score_std[j]
            j -= 1
        # insert value[i] to its position
        dev_days[j+1] = day_temp; score[j+1] = score_temps; score_std[j+1] = score_std_temps

    return dev_days, score, score_std

#--------------------------------------------------------------------------#
if len(sys.argv) != 3:
    exit("arg error - need 2 arg: [folder with position files] [absolute path of output file]")
else:
    main(sys.argv[1], sys.argv[2])
    print("done")
