#	Code to plot the distribution of curvatures
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 15 February 2021

#  Load the required packages
import os
import numpy as np

# systems
msystems = ["pure-membrane","with-phd-random","with-phd-collar"]

count = 0
for i in msystems:
	
	mdata = np.loadtxt("./mean-curvature-data-"+i+"/frame-"+str(0)+".csv",delimiter=",",skiprows=1)
	mean_curv = mdata[:,0]
	xyz = mdata[:,[1,2,3]]

	
	for j in range(1,1001):
		mdata = np.loadtxt(i+"/mean-curvature-data-"+i+"/frame-"+str(j)+".csv",delimiter=",",skiprows=1)
		mean_curv = np.column_stack([mean_curv,mdata[:,0]])
	mean_mean_curv = np.mean(mean_curv,axis = 1)
	mean_mean_curv = np.column_stack([xyz,mean_mean_curv])
	np.savetxt(i+"-time-averaged-mean-curvature.csv",mean_mean_curv,delimiter=',',header
='x,y,z,mean_curvature')
