# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:37:15 2016

@author: matthias
"""

import csv
import numpy as np
import matplotlib.pyplot as plt

path="/Users/matthias/Documents/popdyn/botero-model/Output to analyze/botero compare/R=1000/P=0.1/pop1_std_genes.csv"
std=[]
with open(path) as f:
    reader = csv.reader(f,delimiter=",")
    for (i,row) in enumerate(reader):
        if row:
            if row[0]=="n":
                data = np.genfromtxt(path,skip_header=i+1,delimiter=",") #reads mean genes and n , nperPos from csv file
                break

std_genes = data[:,:-1]  #all genes
std_genes = np.fabs(std_genes) #gives absolute value
std_genes=np.delete(std_genes,1,1) #delete env
std_genes=np.delete(std_genes,3,1) #delete M
std_genes=np.delete(std_genes,7,1) #delete m
std_genes=np.delete(std_genes,7,1) #delete ma

plt.figure()
labels=["n","$I_0$","$I_0'$","$a$","$b$","$b'$","$h$","$s$"]
for i in range(1,len(std_genes[0])):
    plt.plot(std_genes[:,0],std_genes[:,i],label=("{0}".format(labels[i])),linewidth='1')
plt.legend()
plt.xlabel("Generation")
plt.ylabel("Standards deviation")
#plt.yscale("log")
#plt.xscale("log")
plt.savefig("std_plot.pdf")


    
std_genes=np.delete(std_genes,2,1) #delete I0'
std_genes=np.delete(std_genes,4,1) #delete b'
#std_genes=np.delete(std_genes,5,1) #delete s

plt.figure()
labels=["n","$I_0$","$a$","$b$","$h$","$s$"]
for i in range(1,len(std_genes[0])):
    plt.plot(std_genes[:,0],std_genes[:,i],label=("{0}".format(labels[i])),linewidth='1')
plt.legend()
plt.xlabel("Generation")
plt.ylabel("Standard deviation")
#plt.yscale("log")
#plt.xscale("log")
plt.savefig("std_plot2.pdf")


    
