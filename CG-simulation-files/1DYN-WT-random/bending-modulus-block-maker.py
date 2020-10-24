#	Code to create blocks of spectra data by shifting window method
#   
#       Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 
#	For any queries regarding usage of the script please feel 
#       free to contact me on 
#       krishnakanth.baratam@gmail.com, krishnakanth@iisc.ac.in

#  Load the required packages
import os
import matplotlib.pylab as plt
import numpy as np

# input .xyz file
traj = "centered-production_po4.xyz"

# number of frames available in the .xyz file
n_frames = 30228

# correction factor on structure factor obtained from cpp code
c_factor = 10.24 # possibly s(q)*n_lipids/1000

# shifting window jump
jump = 5000

#window_size
w_size = 10000

os.system("mkdir block-analysis")

i = 0
counter = 0
while (i+w_size) <= n_frames:
	os.system("./a.out "+traj+" "+"block-analysis/block-"+str(counter) +".dat "+str(i)+" "+str(i+w_size))
	i = i+jump
	counter  = counter + 1

	
	

