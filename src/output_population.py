#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#   output_population.py
#   Author: Dion Häfner (dionhaefner@web.de)
#   
#   Responsible for output and plotting
#   
#   Licensed under BSD 2-Clause License
#
#########################################################
"""

# Import third-party libraries
import matplotlib as mpl
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import warnings
import math
#rom numba import jit
# Seaborn makes prettier plots, but is not installed in a fresh Anaconda python
try: 
    import seaborn as sns 
    have_seaborn = True
except ImportError:
    have_seaborn = False
mpl.rcParams["axes.formatter.useoffset"] = False
# Import other parts of the project
#from animal import Animal
#from population import Population
#from environment import Environment
from constants import model_constants
constants = model_constants
zeros=len(str(constants["generations"]))
#@jit
def output_population(population,f1,f2,j,k,path,force_plot,t,env,sizes,variable=False):
    """
    Outputs state of the Population. Inputs:
        population: instance of Population to be output, 
        f1: file handle for mean gene file, f2: file handle for standard dev. of genes
        j: current generation counter, k: current population counter,
        path: output path, timeseries: whether complex output should be saved, #???
        t: current time step, env: list of environment
    """

    animals = np.array(population.animals())

    genes = list(map(lambda x: x.gene_dict, animals))
    #print(genes) 

    for a,animal in enumerate(animals):
        genes[a]["W"]=animal.lifetime_payoff() #add mismatch for plotting later
    n= len(genes)
    data = pd.DataFrame(genes)
    mean = pd.DataFrame(data.mean())
    mean[0]["s"]=data["s"].median() #use median instead of mean for gene s
    mean[0]["t"]=data["t"].median()
    mean=mean.transpose()   
    std = pd.DataFrame(data.std()).transpose()
   


    f1.write(str(j)+",") #generation
    f2.write(str(j)+",")

    if n == 0:
        f1.write("0,0,0,0,0,0,0,0,0,0,0")
        f2.write("0,0,0,0,0,0,0,0,0,0,0")
    else:
        mean.to_csv(f1, header=False, index=False, line_terminator='')
        std.to_csv(f2, header=False, index=False, line_terminator='')
    if variable:     
        f1.write(","+str(n)) #number of animals
        f2.write(","+str(n))
    f1.write(","+str(max(np.bincount(population.lineage())))) #size of biggest family (those with same lineage)

    f1.write("\n")
    f2.write("\n")
    if constants["proc"]>1:
        dtype="png"
    else:
        dtype=constants["format"]
        
    filename = path+'timeseries/pop'+str(k+1)+'_genes_'+str(j).zfill(zeros)+'.'+dtype
    if force_plot:
        plot_situation(j,data,n,env,filename,sizes,variable)
    elif constants["plot_every"] > 0:
        if (j % constants["plot_every"]) == 0: #modulo to plot every n times
            plot_situation(j,data,n,env,filename,sizes,variable)
    elif constants["plot_every"] < 0:
        T=math.ceil(constants["environment"][0]/6)  #if plot_every is set smaller than 0, plot 6 times per environment cycle
        if (j % T) == 0: 
            plot_situation(j,data,n,env,filename,sizes,variable)
    return mean, std

#@jit
def plot_situation(j,data,n,env,filename,sizes,variable=False):
    constants = model_constants
    if constants["verbose"]:
        print("\nPlotting ...")
    fsize=12  #fontsize of axis labels and tick labels
    #palette = sns.color_palette("Set2", 4)
    plt.figure(figsize=(10,18)) #size adjustment to have axis labels visible
    rows=4
    index=0
    if variable:
        rows=5
        index=1
        ax = plt.subplot2grid((rows,1),(0,0))       
        ax.plot(sizes,"-",lw=0.7)
        ax.set_xlim(0,j+1)
        ax.set_xlabel("Time",fontsize=fsize)
        ax.set_ylim(0,constants["size"]+200)
        ax.set_ylabel("Size",fontsize=fsize)
        plt.tick_params(axis='both', which='both', labelsize=fsize)
       


    if (n> 0):
        ax = plt.subplot2grid((rows,1),(index,0),rowspan=2)
    
        if have_seaborn:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sns.violinplot(data=data,ax=ax,scale='width',linewidth=0.6)  #gene plot
        else:
            data.boxplot(ax=ax)

        ax.set_ylim(-2,2)
        ax.set_xlabel("Genes",fontsize=fsize)
        ax.set_ylabel("Value distribution",fontsize=fsize)
        plt.tick_params(axis='both', which='both', labelsize=fsize)
        ax1 = plt.subplot2grid((rows,1),(index+2,0),rowspan=2)
        scale = 5*constants["L"]*env.R #5 whole cycles per plot
        if j <= constants["L"]*constants["generations"]:

            j0 = np.arange(min(max(0,j-scale/2), np.abs(constants["L"]*constants["generations"]-scale)) ,min(constants["L"]*constants["generations"],max(j+scale/2,scale)),0.01*env.R*constants["L"]) #star always in the middle except at beginning and end
        else:
            j0 = np.arange(j-scale/2,j+scale/2,0.01*env.R*constants["L"])
        ax1.plot(j0,np.array(list(map(env.evaluate,j0)))[:,0],linewidth=0.7) #plot E(t)
        ax1.scatter(j,env.evaluate(j)[0],s=100,marker='*') #show actual time as *
        ax1.set_ylim(-2,2)
        ax1.set_xlim(j0[0],j0[-1]) 
        ax1.set_xlabel("Time t",fontsize=fsize)
        ax1.set_ylabel("E",fontsize=fsize)#labels und überschriften
        plt.tick_params(axis='both', which='both', labelsize=fsize)
    plt.suptitle("The situation at generation"+str(j),fontsize=17)
    plt.subplots_adjust(top=0.95)
    plt.savefig(filename)
    plt.close()

#@jit
def plot_size(path,fi,k): #plots the number of animals in each environment for each generation
    constants = model_constants
    f = open(fi)
    for (i,row) in enumerate(f):
        if row[0]=="n":
            data = np.genfromtxt(fi,skip_header=i+1,delimiter=",")
            break

    sizes = data[:,-2] #[:,-2] gives all elements in the last column (n per position)                                
    plt.figure()
    plt.plot(sizes[:],alpha=0.7,label="Environment ",linewidth=0.5)  #sizes[:,i] gives the elements of the ith column
    plt.legend()
    plt.ylim(0,constants["size"]+200)
    plt.xlabel("Generation")
    plt.ylabel("Number of individuals")
    if constants["proc"]>1:
        dtype="png"
    else:
        dtype=constants["format"]
    plt.savefig(str(path)+"sizes_"+str(int(k)+1)+"."+dtype,bbox_inches='tight')
    plt.close()
