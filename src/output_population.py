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

# Import other parts of the project
#from animal import Animal
#from population import Population
#from environment import Environment
from constants import model_constants

#@jit
def output_population(population,f1,f2,j,k,path,force_plot,t,env,sizes,times,variable=False):
    """
    Outputs state of the Population. Inputs:
        population: instance of Population to be output, 
        f1: file handle for mean gene file, f2: file handle for standard dev. of genes
        j: current generation counter, k: current population counter,
        path: output path, timeseries: whether complex output should be saved, #???
        t: current time step, env: list of environments
    """
    constants = model_constants
    animals = np.array(population.animals())

    genes = list(map(lambda x: x.gene_dict, animals))
    #print(genes) 

    for a,animal in enumerate(animals):
        genes[a]["M"]=animal.mismatch  #add mismatch for plotting later
    n= len(genes)
    data = pd.DataFrame(genes)
    mean = pd.DataFrame(data.mean()).transpose()   
    std = pd.DataFrame(data.std()).transpose()
   


    f1.write(str(j)+",") #generation, environment
    f2.write(str(j)+",")

    if n == 0:
        f1.write("0,0,0,0,0,0,0,0,0")
        f2.write("0,0,0,0,0,0,0,0,0")
    else:
        mean.to_csv(f1, header=False, index=False, line_terminator='')
        std.to_csv(f2, header=False, index=False, line_terminator='')
    
    f1.write(","+str(n)+","+str(max(np.bincount(population.lineage())))) #last numbers: number of animals and size of biggest family
    f2.write(","+str(n))

    f1.write("\n")
    f2.write("\n")

    filename = path+'timeseries/pop'+str(k+1)+'_genes_'+str(j)+'.pdf'
    if force_plot:
        plot_situation(t,data,n,env,filename,sizes,times,variable)
    elif constants["plot_every"] > 0:
        if (j % constants["plot_every"]) == 0: #modulo to plot every n times
            plot_situation(t,data,n,env,filename,sizes,times,variable)
    elif constants["plot_every"] < 0:
        T=math.ceil(constants["environments"][0]/6)  #if plot_every is set smaller than 0, plot 6 times per environment cycle
        if (j % T) == 0: 
            plot_situation(t,data,n,env,filename,sizes,times,variable)
    return mean, std

#@jit
def plot_situation(t,data,n,env,filename,sizes,times,variable=False):
    constants = model_constants
    if constants["verbose"]:
        print("\nPlotting ...")
    fsize=37  #fontsize of axis labels and tick labels
    palette = sns.color_palette("Set2", 4)
    plt.figure(figsize=(30,40)) #size adjustment to have axis labels visible
    rows=4
    index=0
    if variable:
        rows=5
        index=1
        ax = plt.subplot2grid((rows,1),(0,0))       
        ax.plot(times,sizes,"-")
        ax.set_xlim(times[0],t+1)
        ax.set_xlabel("Time",fontsize=fsize)
        ax.set_ylim(0,constants["environment_sizes"]+200)
        ax.set_ylabel("Size",fontsize=fsize)
        plt.tick_params(axis='both', which='both', labelsize=fsize)
       


    if (n> 0):
        ax = plt.subplot2grid((rows,1),(index,0),rowspan=2)
    
        if have_seaborn:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sns.violinplot(data=data,ax=ax,scale='width')  #gene plot
        else:
            data.boxplot(ax=ax)

        ax.set_ylim(-2,2)
        ax.set_xlabel("Genes",fontsize=fsize)
        ax.set_ylabel("Value distribution",fontsize=fsize)
        plt.tick_params(axis='both', which='both', labelsize=fsize)
        ax1 = plt.subplot2grid((rows,1),(index+2,0),rowspan=2)
        scale = 5*constants["L"]*env.R #5 whole cycles per plot
        if t <= constants["L"]*constants["generations"]:

            t0 = np.arange(min(max(0,t-scale/2), np.abs(constants["L"]*constants["generations"]-scale)) ,min(constants["L"]*constants["generations"],max(t+scale/2,scale))) #star always in the middle except at beginning and end
        else:
            t0 = np.arange(t-scale/2,t+scale/2)
        ax1.plot(t0,np.array(list(map(env.evaluate,t0)))[:,0]) #plot E(t)
        ax1.scatter(t,env.evaluate(t)[0],s=250,marker='*') #show actual time as *
        ax1.set_ylim(-2,2)
        ax1.set_xlim(t0[0],t0[-1]) 
        ax1.set_xlabel("Time t",fontsize=fsize)
        ax1.set_ylabel("E",fontsize=fsize)#labels und überschriften
        plt.tick_params(axis='both', which='both', labelsize=fsize)

    plt.suptitle("The situation at $t = $"+str(t),fontsize=50)
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
    plt.ylim(0,constants["environment_sizes"]+200)
    plt.xlabel("Generation")
    plt.ylabel("Number of individuals")
    plt.savefig(str(path)+"sizes_"+str(int(k)+1)+".pdf",bbox_inches='tight')
    plt.close()
