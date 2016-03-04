# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:37:15 2016

@author: matthias
"""

import csv
import numpy as np
import matplotlib.pyplot as plt

path="/Users/matthias/Documents/popdyn/botero-model/fair_stop/output/16-03-04_22-39-55/pop1_std_genes.csv"
std=[]
with open(path) as f:
    reader = csv.reader(f,delimiter=",")
    for (i,row) in enumerate(reader):
        if row:
            if row[0]=="n":
                data = np.genfromtxt(path,skip_header=i+1,delimiter=",") #reads mean genes and n , nperPos from csv file
                break
names=np.genfromtxt(m,skip_header=0,delimiter=",",max_rows=1,dtype=str) #reads mean genes and n , nperPos from csv file
std = np.genfromtxt(s,skip_header=1,delimiter=",") #reads mean genes and n , nperPos from csv file

names=np.delete(names,8)
names=np.delete(names,8)
std_genes = data[:,:]  #all genes
std_genes = np.fabs(std_genes) #gives absolute value
std_genes=np.delete(std_genes,1,1) #delete env
std_genes=np.delete(std_genes,3,1) #delete M
std_genes=np.delete(std_genes,7,1) #delete m
std_genes=np.delete(std_genes,7,1) #delete ma
std_genes=np.delete(std_genes,8,1) #delete nperPos

plt.figure()
labels=["n","$I_0$","$I_0'$","$a$","$b$","$b'$","$h$","$s$","Lineage"]
for i in range(1,len(std_genes[0])):
    plt.plot(std_genes[:,0],std_genes[:,i],label=("{0}".format(labels[i])),linewidth='1')
plt.legend()
plt.xlabel("Generation")
plt.ylabel("Standards deviation")
#plt.yscale("log")
#plt.xscale("log")
plt.savefig("std_plot.pdf")


    
delete=["I_0","I_0p","a","b","bp","h","s"]
for d in delete:
    for i,n in enumerate(names):
        if n==d:
           names=np.delete(names,i)
           mean_genes=np.delete(mean_genes,i-1,1)
           std_genes=np.delete(std_genes,i-1,1)
           break

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


    
