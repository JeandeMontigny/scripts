#!/usr/bin/env python3
import sys, os, re, numpy
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------#
def main(output_file):
# dico={}#TODO: create a dico with the cell sub-type as key
# dico[cell_type]=soma_x()
	tab_x_On=[]; tab_x_Off=[]; tab_x_OnOff=[]; tab_x_other=[]
	tab_y_On=[]; tab_y_Off=[]; tab_y_OnOff=[]; tab_y_other=[]
	fichier=open(output_file, "r")
	file_lines=fichier.readlines()
	for  line in file_lines:
		if line[0]!="c":
            # 1: cell type; 2: x; 3: y
			m=re.search( r'.+ .+ .+ .+ .+ ([1-4]+) .+ (.+) (.+)', line, re.M|re.I)
			if m:
				if int(m.group(1))==1: # ON-OFF
					tab_x_OnOff.append(int(m.group(2))); tab_y_OnOff.append(int(m.group(3)))
				elif int(m.group(1))==2: # ON
					tab_x_On.append(int(m.group(2))); tab_y_On.append(int(m.group(3)))
				elif int(m.group(1))==3: # OFF
					tab_x_Off.append(int(m.group(2))); tab_y_Off.append(int(m.group(3)))
				elif int (m.group(1))==4: # other
					tab_x_other.append(int(m.group(2))); tab_y_other.append(int(m.group(3)))

    # cells distance for every cells
	on_cells=distance(tab_x_On, tab_y_On)
	off_cells=distance(tab_x_Off, tab_y_Off)
	on_off_cells=distance(tab_x_OnOff, tab_y_OnOff)
	other_cells=distance(tab_x_other, tab_y_other)

	if (sys.argv[2]=="all"):
		print(len(on_cells), "ON cells:")
		figure_all(on_cells)
		print(len(off_cells), "OFF cells:")
		figure_all(off_cells)
		print(len(on_off_cells), "ON-OFF cells:")
		figure_all(on_off_cells)
	figure_ave(on_cells, off_cells, on_off_cells, other_cells)

#--------------------------------------------------------------------------#
def figure_ave(on, off, on_off, other):
	fig_ave=plt.figure(1)
	if len(on)!=0:
		sub1=fig_ave.add_subplot(2,2,1); plt.hist(ave_sorted(on)[0]); plt.xlim(0, plt.xlim()[1]); sub1.set_title("ON cells")
	if len(off)!=0:
		sub2=fig_ave.add_subplot(2,2,2); plt.hist(ave_sorted(off)[0]); plt.xlim(0, plt.xlim()[1]); sub2.set_title("OFF cells")
	if len(on_off)!=0:
		sub3=fig_ave.add_subplot(2,2,3); plt.hist(ave_sorted(on_off)[0]); plt.xlim(0, plt.xlim()[1]); sub3.set_title("ON-OFF cells")
	if len(other)!=0:
		sub4=fig_ave.add_subplot(2,2,4); plt.hist(ave_sorted(other)[0]); plt.xlim(0, plt.xlim()[1]); sub4.set_title("other cells")
	plt.show()

#--------------------------------------------------------------------------#
def figure_all(cells_tab):
	tab_all_closest=[]
	fig=plt.figure(1)
#	nb=np.ceil(np.sqrt(len(cells_tab))) # how many row and column to create for every cell to fit in the plot
	for i in range(0, len(cells_tab)):
		min_dist=min(cells_tab[i])
		print("closest cell:", min_dist, "for cell", i+1)
		tab_all_closest.append(min_dist)
		# plot cell distance histogram for every cell
#		fig.add_subplot(nb, nb, i+1) # row, column, nb of subplot
#		plt.hist(cells_tab[i])
#		plt.xlim(0, np.amax(cells_tab))
	plt.hist(tab_all_closest)
	plt.show()

#--------------------------------------------------------------------------#
def ave_sorted(cells_tab):
	cells_tab_ave=np.average(np.sort(cells_tab), axis=0)
	cells_tab_std=np.std(np.sort(cells_tab), axis=0)
	print("average closest cell:", round(cells_tab_ave[0], 1), "with std of", round(cells_tab_std[0], 1), "with", len(cells_tab), "cells")
	print("regularity index (mean/std):", round(cells_tab_ave[0]/cells_tab_std[0], 2), "(higher than 2-3 indicates non random distibution)")
	return cells_tab_ave, cells_tab_std

#--------------------------------------------------------------------------#
def distance(tab_x, tab_y):
	if (len(tab_x)!= len(tab_y)):
		exit("tab x and y should be of the same size")
	tab_cells_dist=[]
	for i in range(0, len(tab_x)):
		tab_cell=[]
		for j in range(0, len(tab_x)):
			temps_distance=round(np.sqrt(np.square(abs(tab_x[i]-tab_x[j]))+np.square(abs(tab_y[i]-tab_y[j]))), 1)
			if temps_distance!=0.0:
				tab_cell.append(temps_distance)
		tab_cells_dist.append(tab_cell)
	return tab_cells_dist

#--------------------------------------------------------------------------#
if len(sys.argv)==3 and (sys.argv[2]=="all" or sys.argv[2]=="-"):
	main(sys.argv[1])
	print("done")
else:
	exit("arg error - need 2 arg: [measure_output file] [all,-]")

