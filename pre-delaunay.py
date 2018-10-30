#!/usr/bin/env python3
import sys, os, re

#--------------------------------------------------------------------------#
def main(path_dir):
	os.chdir(path_dir)
	if os.path.isfile("ON_cells.txt"):
		os.remove("ON_cells.txt")
		os.remove("OFF_cells.txt")
	sortie_on=open("ON_cells.txt", "a")
	sortie_off=open("OFF_cells.txt", "a")

	fichier=open("measures_output.txt", "r")
	file_lines=fichier.readlines()
	for  line in file_lines:
		if line[0]=="r":
            # 1: ON/OFF; 2: x; 3: y
			m=re.search(".+type=([ONF]+).+ .+ .+ .+ .+ .+ (.+) (.+)", line)
			if m:
				if m.group(1)=="ON":
					sortie_on.write(str(m.group(2))+" "+str(m.group(3))+"\n")
				else:
					sortie_off.write(str(m.group(2))+" "+str(m.group(3))+"\n")

#--------------------------------------------------------------------------#
if len(sys.argv)==2:
	main(sys.argv[1])
	print("done")
else:
	exit("arg error - need 2 arg: [directory to measure_output file]")
