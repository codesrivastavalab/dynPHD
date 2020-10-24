#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt


#rescolors = ['C0','C0','C1','C2','C3',"C4","C5","C6"]
#req_systems = ["F579A","F579A-2","F579A-S581K","I533A-F579A","I533A","WT","K583A"]
rescolors = ['C4','C0','C1','C2']
req_systems = ["WT","F579A-S581K","F579A-K583A"]
#plt.rcParams.update({'font.size':25})
#plt.rcParams['axes.linewidth'] = 3
#plt.figure(figsize=(12,10))
#plt.rcParams['xtick.major.size'] = 15
#plt.rcParams['xtick.major.width'] = 4
#plt.rcParams['ytick.major.size'] = 15
#plt.rcParams['ytick.major.width'] = 4

#small function to find the start of the datavalues
#-------------------------------------#
def start_finder(x,mpattern):
	loc=0
	for j in range(50):
		if(str.find(x[j],mpattern)!=-1):
			loc=j
	return(loc+1)
#-------------------------------------#
counter = 0
for i in req_systems:
    temp = open(i+"-BS-output.xvg").readlines()
    temp=start_finder(temp,"@")
    us_data = np.loadtxt(i+"-BS-output.xvg",skiprows = temp)
    us_data[:,1] = us_data[:,1]/ 4.184 # kcal correction on profile
    us_data[:,2] = us_data[:,2]/ 4.184 # kcal correction on errors
    #print(np.min(us_data[:,1]))
    us_data[:,1] = us_data[:,1] - np.min(us_data[:,1])
    plt.plot(us_data[:,0],us_data[:,1],color=rescolors[counter],label =i)
    plt.fill_between(us_data[:,0],us_data[:,1]-us_data[:,2],us_data[:,1]+us_data[:,2],color =rescolors[counter],alpha = 0.5)
    counter = counter + 1



plt.xlabel('CV (nm) ')
plt.ylabel('Energy (kcal/mol)')
plt.legend(loc='lower right')
plt.savefig("double-mutants-pmf-plots.tiff",dpi  = 300)






