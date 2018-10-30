#!/usr/bin/env python
import sys, os, re, numpy, peakutils, shutil, random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D

global tab_of_tab_z
tab_of_tab_z=[]

#TODO: clean and comment!

#--------------------------------------------------------------------------#
def main(path_dir):
    nb_file=0; bary_x=[]; bary_y=[]; bary_z=[]; list_branching=[]
    list_isometrie=[]; list_distance_branching=[]; list_cell_type=[]; list_distance_terminal=[]
    list_thickness=[]; list_disc_diam=[]; list_nb_dend_root=[]
#    list_peak_distance=[]

    os.chdir(path_dir)
#    if not os.path.isdir("./ON"):
#        os.mkdir("./ON"); os.mkdir("./OFF"); os.mkdir("./ON_OFF"); os.mkdir("./OTHER")
    if os.path.isfile("measures_output.txt"):
        os.remove("measures_output.txt")
    sortie=open("measures_output.txt", "a")
#    sortie.write("cell_ID barry_center_x barry_center_y barry_center_z nb_branching anisometry ave_branching_dist ave_branching_dist_std ave_tip_dist ave_tip_dist_std disc_thickness disc_diam_ave disc_diam_95 nb_dend_root\n")
    sortie.write("cell_ID"+" "+"anisometry"+" "+"disc_diam"+" "+"branching_dist"+" "+"tip_dist"+" "+"cell_type"+" "+"branching_nb"+" "+"soma_x_coord"+" "+"soma_y_coord"+"\n")

    print "processing swc files.."
    number_of_file = len([name for name in os.listdir('.') if os.path.isfile(name)])
    for file in os.listdir("./"):
        if file.endswith(".swc"):
            nb_file+=1
            sys.stdout.write("\r"+str(nb_file)+"/"+str(number_of_file))
            sys.stdout.flush()
            temps_lecture=lecture(file)
            if temps_lecture:
                if len(temps_lecture)==2: # if only a cell body
                    sortie=open("measures_output.txt", "a")
                    sortie.write(file+" "+str(0)+" "+str(0)+" "+str(0)+" "+str(0)+" "+str(4)+" "+str(0)+" "+str(temps_lecture[0])+" "+str(temps_lecture[1])+"\n")
                if len(temps_lecture)==12: # if no branching point in this cell
                    bary_x.append(temps_lecture[0]); bary_y.append(temps_lecture[1]); bary_z.append(temps_lecture[2]); list_branching.append(float(temps_lecture[3]))
                    list_isometrie.append(temps_lecture[4]); list_cell_type.append(temps_lecture[5])
                    list_distance_terminal.append(temps_lecture[6]); list_thickness.append(temps_lecture[7]); list_disc_diam.append(temps_lecture[8])
        #            list_peak_distance.append(temps_lecture[9])
                    list_nb_dend_root.append(temps_lecture[9])
                    #TODO write in measure_output the measures when the cell has no branching point.
                    #sortie=open("measures_output.txt", "a")
                    #sortie.write(file+" "+str(temps_lecture[4])+" "+str(temps_lecture[9])+" "+str(temps_lecture[5])+" "+str(temps_lecture[7])+" "+str(temps_lecture[6])+" "+str(temps_lecture[3])+" "+str(temps_lecture[11])+" "+str(temps_lecture[12])+"\n")
#                else:
#                    sortie.write(file+" "+str(temps_lecture[0])+" "+str(temps_lecture[1])+" "+str(temps_lecture[2])+" "+str(temps_lecture[3])+" "+str(temps_lecture[4])+" "+str(temps_lecture[5])+" "+str(temps_lecture[6])+" "+str(temps_lecture[7])+" "+str(temps_lecture[8])+" "+str(temps_lecture[9])+" "+str(temps_lecture[10])+" "+str(temps_lecture[11])+" "+str(temps_lecture[12])+"\n")
                if len(temps_lecture)==13:
                    bary_x.append(temps_lecture[0]); bary_y.append(temps_lecture[1]); bary_z.append(temps_lecture[2]); list_branching.append(float(temps_lecture[3]))
                    list_isometrie.append(temps_lecture[4]); list_distance_branching.append(temps_lecture[5]); list_cell_type.append(temps_lecture[6])
                    list_distance_terminal.append(temps_lecture[7]); list_thickness.append(temps_lecture[8]); list_disc_diam.append(temps_lecture[9])
        #            list_peak_distance.append(temps_lecture[10])
                    list_nb_dend_root.append(temps_lecture[10])
                    sortie=open("measures_output.txt", "a")
                    sortie.write(file+" "+str(temps_lecture[4])+" "+str(temps_lecture[9])+" "+str(temps_lecture[5])+" "+str(temps_lecture[7])+" "+str(temps_lecture[6])+" "+str(temps_lecture[3])+" "+str(temps_lecture[11])+" "+str(temps_lecture[12])+"\n")
    else:
        print "\ncould not open a file"

    nb_b=round(numpy.average(list_branching), 3); std_nb_b=round(numpy.std(list_branching), 3)
    bar_x=round(numpy.average(bary_x), 3); bar_y=round(numpy.average(bary_y), 3); bar_z=round(numpy.average(bary_z), 3)
    std_bar_x=round(numpy.std(bary_x), 3); std_bar_y=round(numpy.std(bary_y), 3); std_bar_z=round(numpy.std(bary_z), 3)
    iso=round(numpy.average(list_isometrie), 3); std_iso=round(numpy.std(list_isometrie), 3)
    dist_branch=round(numpy.average(list_distance_branching), 3); std_dist_branch=round(numpy.std(list_distance_branching), 3)
    dist_terminal=round(numpy.average(list_distance_terminal), 3); dist_terminal_std=round(numpy.std(list_distance_terminal), 3)
    thickness_ave=round(numpy.average(list_thickness), 3); thickness_std=round(numpy.std(list_thickness), 3)
    disc_diam_ave=round(numpy.average(list_disc_diam), 3); disc_diam_std=round(numpy.std(list_disc_diam), 3)
    ON_OFF_cell=list_cell_type.count(1); ON_cell=list_cell_type.count(2); OFF_cell=list_cell_type.count(3); OTHER_cell=list_cell_type.count(4)
    nb_dend_root_ave=round(numpy.average(list_nb_dend_root), 3); nb_dend_root_std=round(numpy.std(list_nb_dend_root), 3);

    string_output="---- summary ----\n"+str(nb_file)+" file opened\naverage branching number for this directory: "+str(nb_b)+" (real 94.5)\n\twith standard deviation of "+str(std_nb_b)+" (real 54.7)\naverage coordinates of the barycenter: x= "+str(bar_x)+"; y= "+str(bar_y)+"; z= "+str(bar_z)+" (real x= 2.0; y= 24.4; z= 24.5)\n\twith standard deviation of "+str(std_bar_x)+" "+str(std_bar_y)+" and "+str(std_bar_z)+" (real 26.059 28.162 and 11.919)\naverage anisometry measure: "+str(iso)+" (real 0.3)\n\twith standard deviation of "+str(std_iso)+" (real 0.19)\naverage distance of branching: "+str(dist_branch)+" (real 67.4)\n\twith standard deviation of "+str(std_dist_branch)+" (real 22.1)\naverage tip distance: "+str(dist_terminal)+" (real 96.1)\n\twith standard deviation of "+str(dist_terminal_std)+" (real 29.9)\naverage disc thickness: "+str(thickness_ave)+" (real 21.9)\n\twith standard deviation of "+str(thickness_std)+" (real 7.5)\naverage disc span: "+str(disc_diam_ave)+" (real 125.6)\n\twith standard deviation of "+str(disc_diam_std)+" (real 36.3)\naverage dendrite nember attached to the soma: "+str(nb_dend_root_ave)+" (real 4.5)\n\twith standard deviation of "+str(nb_dend_root_std)+" (real 1.2)\nON-OFF_cell: "+str(ON_OFF_cell)+" - ON_cell: "+str(ON_cell)+" - OFF_cell :"+str(OFF_cell)+" - OTHER_cell: "+str(OTHER_cell)

    print string_output
    sortie.write(string_output)
    sortie.close()

#    print "\n average peak distance:", round(numpy.average(list_peak_distance), 3), "with std of", round(numpy.std(list_peak_distance), 3)
#    print sorted(list_peak_distance)
    
    if sys.argv[3]=="plot":
        plot(nb_b, std_nb_b, iso, std_iso, dist_branch, std_dist_branch, dist_terminal, dist_terminal_std, list_isometrie)

    # ---- cloud figures of dif parameters ---- #
#    figure_cloud(list_branching, list_isometrie)
#    figure_cloud(list_distance_branching, bary_z)
#    figure_cloud(list_isometrie, bary_z)
#    figure_cloud(list_distance_terminal, bary_z)
#    figure_cloud(list_distance_terminal, list_isometrie)

#    figure_cloud(list_cell_type, bary_z)
#    figure_cloud(list_cell_type, list_isometrie)
    # ---- IN 3D!!! ---- #
#    figure_cloud_3d(list_isometrie, list_distance_terminal, list_thickness)
#    figure_cloud_3d(list_isometrie, list_disc_diam, list_thickness)
#    pca_analysis(list_isometrie, list_distance_terminal, list_thickness)
#    pca_analysis(list_isometrie, list_thickness, list_disc_diam, list_branching)

    # ---- I don't remember ---- #
#    data=np.array([list_isometrie, list_thickness, list_disc_diam])
#    resultats=find_centers(data, 3)

    # ---- pca attempt ---- #
#    pca_analysis(list_isometrie, list_distance_terminal, list_distance_branching, list_branching)

#--------------------------------------------------------------------------#
def lecture(file_name):
    #TODO: with open("./"+file_name, "r") as fichier:
    fichier=open("./"+file_name, "r")
    if 'fichier' in locals():
        file_lines=fichier.readlines()
        nb_dend_seg=0; nb_ligne=0; branching_index=[]; prev_parent=0; terminasion_index=[]; prev_id=0; line_non_cell=0
        mean_x=0; mean_y=0; mean_z=0; iso_part=[]
        distance_branching_tab=[]; distance_terminal_tab=[]; tab_z=[]; nb_dend_root=0
        x_soma=0; y_soma=0; z_soma=0
    #    prev_x_coord=0; prev_y_coord=0; prev_z_coord=0; dico_branching_point_coord={}; total_dend_length=0
        for  line in file_lines:
            m=re.search( r' ?([0-9]+) ([0-9]+) (.+) (.+) (.+) (.+) ([-0-9]+).?', line, re.M|re.I)
            if m and (line[0]!="#" or line[0]!="\n"):
                nb_ligne+=1
                if int(m.group(7))==-1:# if it's the soma
                    x_soma=float(m.group(3)); y_soma=float(m.group(4)); z_soma=float(m.group(5))
    #                mosaic_tab_x.append(round(x_soma, 1)), mosaic_tab_y.append(round(y_soma, 1)), mosaic_tab_z.append(round(z_soma, 1))
                if int(m.group(7))==1: # if the parent is the soma
                    nb_dend_root+=1
    #                prev_x_coord=0; prev_y_coord=0; prev_z_coord=0
                if int(m.group(2))==3 or int(m.group(2))==5 or int(m.group(2))==6: # if the point is a dendrite or a branching point or a terminal point
                    nb_dend_seg+=1
                    x_coord=float(m.group(3))-x_soma; y_coord=float(m.group(4))-y_soma; z_coord=float(m.group(5))-z_soma
                    tab_z.append(round(z_coord, 1)) # count z
                    if int(m.group(2))==5: # if it is a branching point
                        distance_branching_tab.append(distance_point(x_coord, y_coord, z_coord))
    #                    dico_branching_point_coord.update({int(m.group(1)): [x_coord, y_coord, z_coord]})
                    if int(m.group(2))==6: # if it is a terminal point
                        distance_terminal_tab.append(distance_point(x_coord, y_coord, z_coord))
                    iso_part.append(iso_distributor(x_coord, y_coord))
                    mean_x+=float(m.group(3)); mean_y+=float(m.group(4)); mean_z+=float(m.group(5)) # center of mass
                if (int(m.group(7))<prev_parent) and (int(m.group(7))>1): # if the point is not a son of the previous point and don't belong to the soma
                    branching_index.append(int(m.group(7))) # it means that the parent of that point is a branching point
                    terminasion_index.append(prev_id) # and that the previous point is a terminaison point
                prev_parent=int(m.group(7))
                prev_id=int(m.group(1))
    #            if int(m.group(7) in dico_branching_point_coord.keys():
    #                 coord=dico_branching_point_coord[int(m.group(7)]; prev_x_coord=coord[0]; prev_y_coord=coord[1]; prev_z_coord=coord[2]
    #            total_dend_length+=distance_point(x_coord, y_coord, z_coord)-distance_point(prev_x_coord, prev_y_coord, prev_z_coord)
    #            prev_x_coord=x_coord; prev_y_coord=y_coord; prev_z_coord=z_coord;
            elif line[0]=="#" or line[0]=="\n":
                line_non_cell+=1

        if nb_dend_seg==0:
            print "\tfile only contains a cell body!"
    #        return None
            return int(round(x_soma)), int(round(y_soma))
        isometry_val=isometry(iso_part)
        type_finder_results=type_finder(tab_z); cell_type=type_finder_results[0]; thickness=type_finder_results[1]#; peak_distance=max(type_finder_results[2])
        disc_diam=disc_span(distance_terminal_tab); disc_diam_95=disc_span_95(distance_terminal_tab); average_terminal_distance=round(numpy.average(distance_terminal_tab), 3)
        average_terminal_distance_std=round(numpy.std(distance_terminal_tab), 3)

    #    copy_by_type(cell_type, file_name) # copying cells into diffrent type folders
        if len(distance_branching_tab)==0: # if no branching point in this cell
            print "\tno branching point for this cell!"
            return mean_x/nb_dend_seg, mean_y/nb_dend_seg, numpy.absolute(mean_z/nb_dend_seg), len(branching_index), isometry_val, cell_type, average_terminal_distance, thickness, disc_diam, nb_dend_root, x_soma, y_soma#, peak_distance
            fichier.close()
        else:
            average_branch_distance=round(numpy.average(distance_branching_tab), 3)
            average_branch_distance_std=round(numpy.std(distance_branching_tab), 3)
    #        return round(mean_x/nb_dend_seg, 2), round(mean_y/nb_dend_seg, 2), round(numpy.absolute(mean_z/nb_dend_seg), 2), len(branching_index), isometry_val, average_branch_distance, average_branch_distance_std, average_terminal_distance, average_terminal_distance_std, thickness, disc_diam, disc_diam_95, nb_dend_root, int(round(x_soma)), int(round(y_soma))
            fichier.close()
            return mean_x/nb_dend_seg, mean_y/nb_dend_seg, numpy.absolute(mean_z/nb_dend_seg), len(branching_index), isometry_val, average_branch_distance, cell_type, average_terminal_distance, thickness, disc_diam, nb_dend_root, int(round(x_soma)), int(round(y_soma))#, peak_distance

#--------------------------------------------------------------------------#
#TODO:total dendritic arbor lenght (sum of every dist between 2 point)
def dendrit_lenght():
    return None

#--------------------------------------------------------------------------#
def distance_point(x_coord, y_coord, z_coord):
#TODO: soma is not 0 anymore
    return round(np.sqrt(np.square(x_coord)+np.square(y_coord)+np.square(z_coord)), 3)

#--------------------------------------------------------------------------#
def iso_distributor(x, y):
    if x<0:
        if y<0:
            return 1 if abs(x)>abs(y) else 2
        else:
            return 3 if abs(x)>y else 4
    else:
        if y<0:
            return 5 if x>abs(y) else 6
        else:
            return 7 if x>y else 8

#--------------------------------------------------------------------------#
def isometry(tab):
    theo_val=float(len(tab)/8)
    return round(((abs(tab.count(1)-theo_val)/theo_val)+(abs(tab.count(2)-theo_val)/theo_val)+(abs(tab.count(3)-theo_val)/theo_val)+(abs(tab.count(4)-theo_val)/theo_val)+(abs(tab.count(5)-theo_val)/theo_val)+(abs(tab.count(6)-theo_val)/theo_val)+(abs(tab.count(7)-theo_val)/theo_val)+(abs(tab.count(8)-theo_val)/theo_val))/14, 3)

#--------------------------------------------------------------------------#
def count(tab, inf_boundary, sup_boundary):
    count=0
    for i in tab:
        if abs(i)>=inf_boundary and abs(i)<sup_boundary:
            count+=1
    return count

#--------------------------------------------------------------------------#
def type_finder(tab_z):
    tab_z_count=[]; interval=[]
    for i in frange(0,60, 0.5):
        tab_z_count.append(count(tab_z, i, i+0.5))
        interval.append(i)
    tab_z_count=np.asarray(tab_z_count)
    indexes=peakutils.indexes(tab_z_count, thres=0.2, min_dist=16) # thres=[0., 1.] #real 0.3
#    print "peak(s) distance:", indexes/2
#TODO use the peak distance for type differentiation + std (if std is high, may have more than 1 type) # add std stratification: mean for each cell / peak for intra classe
    thickness=disc_thickness(tab_z_count)
#    thickness=disc_thickness_mean(tab_z)
#    plt.plot(interval, tab_z_count); plt.show() # plot point distance from soma
    tab_of_tab_z.append(tab_z_count)
    if indexes.size==0:
        indexes=numpy.argmax(tab_z_count)
    if indexes.size>1 and np.amin(indexes)<50 and np.amax(indexes)>50: # ON-OFF # 50=25 due to step of 0.5
        return 1, thickness, indexes/2
    elif np.amax(indexes)<50: # ON
        return 2, thickness, indexes/2
    elif np.amin(indexes)>50: # OFF
        return 3, thickness, indexes/2
    else: # other
        return 4, thickness, indexes/2

#--------------------------------------------------------------------------#
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump

#--------------------------------------------------------------------------#
def disc_thickness_mean(tab):
    mean_plus_std=numpy.average(tab)+numpy.std(tab); mean_minus_std=numpy.average(tab)-numpy.std(tab); new_tab=[]
    for dist in tab:
        if dist <= mean_plus_std or dist >= mean_minus_std:
            new_tab.append(dist)
    return max(new_tab)-min(new_tab)

#--------------------------------------------------------------------------#
def disc_thickness(tab):
    tab=tab.tolist()
    tab=filter(lambda a: a!=0, tab)
    dend_nb=sum(tab)
    del tab[0]
    return cut_function_z(tab, dend_nb)

#--------------------------------------------------------------------------#
def cut_function_z(tab, dend_nb):
    if tab[0] > tab[len(tab)-1]:
        del tab[len(tab)-1]
    if tab[0] <= tab[len(tab)-1]:
        del tab[0]
    if (float(sum(tab))/dend_nb)>0.95:#0.95
        return cut_function_z(tab, dend_nb)
    else:
        return float(len(tab))/2 # /2 -> tab step of 0.5

#--------------------------------------------------------------------------#
def disc_span_95(tab):
    tab_length=len(tab); new_tab=[]
    tab.sort()
    for dist in tab:
        if len(new_tab)<0.95*tab_length:
            new_tab.append(dist)
        else:
            return max(new_tab)

#--------------------------------------------------------------------------#
def disc_span(tab):
    mean_plus_std=numpy.average(tab)+numpy.std(tab); new_tab=[]
    for dist in tab:
        if dist <= mean_plus_std:
            new_tab.append(dist)
    return max(new_tab)

#--------------------------------------------------------------------------#
def copy_by_type(cell_type, file_name):
    if cell_type==1:
        shutil.copyfile('./'+file_name, './ON_OFF/'+file_name)
#        print "moved to \"ON_OFF\" folder"
    elif cell_type==2:
        shutil.copyfile('./'+file_name, './ON/'+file_name)
#        print "moved to \"ON\" folder"
    elif cell_type==3:
        shutil.copyfile('./'+file_name, './OFF/'+file_name)
#        print "moved to \"OFF\" folder"
    elif cell_type==4:
        shutil.copyfile('./'+file_name, './OTHER/'+file_name)
#        print "moved to \"OTHER\" folder"

#--------------------------------------------------------------------------#
def plot(nb_b, std_nb_b, iso, std_iso, dist_branch, std_dist_branch, dist_terminal, dist_terminal_std, list_isometrie):
    print "\nconstruction of figures.."

    # ---- classique bar chart ---- #
    figure([0.307], [0.192], [iso], [std_iso], 'score', 'score of anisometry', 1, "y")
    figure([94.504], [54.741], [nb_b], [std_nb_b], 'number', 'number of branching points', 0, "n")
    figure([67.387], [22.17], [dist_branch], [std_dist_branch], 'micrometer', 'branching distance from soma', 0, "n")
    figure([96.077], [29.942], [dist_terminal], [dist_terminal_std], 'micrometer', 'tip distance from soma', 0, "n")

    # ---- anisonetry histrogram and boxplot ---- #
    fig_aniso=plt.figure(5)
    sub1=fig_aniso.add_subplot(2,1,1);
    plt.hist(list_isometrie, orientation="horizontal"); plt.axhline(y=iso, color='r', linestyle='-')
    axes = plt.gca(); axes.set_ylim([0,1]) # set y axis to 0-1
    sub2=fig_aniso.add_subplot(2,1,2); plt.boxplot(list_isometrie, 0, '')
    axes = plt.gca(); axes.set_ylim([0,1]) # set y axis to 0-1

    # ---- average z distance from soma ---- #
    plt.figure()
    tab_mean=[]; interval=[]
    for i in range (0, len(tab_of_tab_z[0])): #120 -> x axis, interval
        temps=0
        for j in range (0, len(tab_of_tab_z)): #y axis, nb of cell
            temps+=tab_of_tab_z[j][i]
        tab_mean.append(round(temps/len(tab_of_tab_z), 1))
    for i in frange(0,60, 0.5):
        interval.append(i)
    plt.plot(interval, tab_mean);

    # ---- show every plots ---- #
    print "construction done"
    plt.show()

#--------------------------------------------------------------------------#
def figure(seung_cells, seung_std, modeled_cells, modeled_std, unit, legende, ylim, legend):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    N=len(seung_cells)

    ind=np.arange(N) # the x locations for the groups
    width=0.35 # the width of the bars
    # the bars
    rects1=ax.bar(ind, seung_cells, width,
            color='black',
            yerr=seung_std,
            error_kw=dict(elinewidth=2,ecolor='red'))

    rects2=ax.bar(ind+width, modeled_cells, width,
            color='red',
            yerr=modeled_std,
            error_kw=dict(elinewidth=2,ecolor='black'))
    # axes and labels
    max_y=max(max(seung_cells)+max(seung_std), max(modeled_cells)+max(modeled_std))
#    ax.set_xlim(-width,len(ind)+(width/2))
    ax.set_xlim(-width, width*2)
    if ylim==1:
        ax.set_ylim(0,1)
    else:
        ax.set_ylim(0,max_y+max_y/10)
    ax.set_ylabel(unit, fontsize=24)
    ax.set_title(legende, fontsize=24)

    xTickMarks=[''] # =['x', 'y', ' z']
    ax.set_xticks(ind+width)
    xtickNames=ax.set_xticklabels(xTickMarks)
#    plt.setp(xtickNames, fontsize=16) #plt.setp(xtickNames, rotation=45, fontsize=16)

    if legend=="y":
        ax.legend((rects1[0], rects2[0]), ('real', 'modeled'), fontsize=20)
    
#--------------------------------------------------------------------------#
def figure_cloud(param1, param2):
#TODO: 3D plot. plot z density peak vs list_isometry vs list_distance_branching
    if len(param1)==len(param2):
        plt.plot(param1, param2, 'o')
#        plot.axis
        plt.show()
    else:
        print "not the same lenght"
#TODO: principal component analysis

#--------------------------------------------------------------------------#
def figure_cloud_3d(param1, param2, param3):
    fig=plt.figure()
    ax=Axes3D(fig)
    pltData=[param1, param2, param3]
    ax.scatter(pltData[0], pltData[1], pltData[2], 'o')
    xAxisLine=((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
    yAxisLine=((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
    zAxisLine=((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.

    ax.set_xlabel("x-axis label") 
    ax.set_ylabel("y-axis label")
    ax.set_zlabel("z-axis label")
    ax.set_title("")
    plt.show()

#--------------------------------------------------------------------------#
def pca_analysis(param1, param2, param3, param4):
    zip_param1=np.array(zip(param1))
    zip_param2=np.array(zip(param2))
    zip_param3=np.array(zip(param3))
    zip_param4=np.array(zip(param4))
    matrix=[]
    for i in range (0, len(zip_param1)):
#        matrix.append([float(zip_param1[i]), float(zip_param2[i]), float(zip_param3[i])])
        matrix.append([float(zip_param1[i]), float(zip_param2[i]), float(zip_param3[i]), float(zip_param4[i])])
    matrix=np.array(matrix)
    results=PCA(matrix)
    print "PCA - variance percentages for each component:", results.fracs
#    print results.Y
    plot3d_pca(results)

#--------------------------------------------------------------------------#
def plot3d_pca(results):
    x=[]; y=[]; z=[]
    for item in results.Y:
        x.append(item[0])
        y.append(item[1])
        z.append(item[2])
    figure_cloud_3d(x,y,z)

#--------------------------------------------------------------------------#
if len(sys.argv)==4 and (sys.argv[2]=="rewrite" or sys.argv[2]=="-") and (sys.argv[3]=="plot" or sys.argv[3]=="-"):
    if not os.path.isdir(sys.argv[1]):
        exit("given directory doesn't exist")
    if sys.argv[2]=="rewrite":
        import correction
        correction.main(sys.argv[1])
        arg=sys.argv[1]+"/rewrite"
        main(arg)
    else:
        main(sys.argv[1])
    print "\ndone"
else:
    exit("arg error - need 3 arg: [absolute path to directory containing .swc files] [rewrite,-] [plot,-]")

