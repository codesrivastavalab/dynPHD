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
	
	mdata = np.loadtxt("./gaussian-curvature-data-"+i+"/frame-"+str(0)+".csv",delimiter=",",skiprows=1)
	gauss_curv = mdata[:,0]
	xyz = mdata[:,[1,2,3]]
	
	for j in range(1,1001):
		mdata = np.loadtxt(i+"/gaussian-curvature-data-"+i+"/frame-"+str(j)+".csv",delimiter=",",skiprows=1)
		gauss_curv = np.column_stack([gauss_curv,mdata[:,0]])
	mean_gauss_curv = np.mean(gauss_curv,axis = 1)
	mean_gauss_curv = np.column_stack([xyz,mean_gauss_curv])
	np.savetxt(i+"-time-averaged-gaussian-curvature.csv",mean_gauss_curv,delimiter=',',header
='x,y,z,gaussian_curvature') 
