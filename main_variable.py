#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#    main_variable.py
#    Author: Dion Häfner (dionhaefner@web.de)
#    
#    Main controller, with variable population size
#
#    Usage:
#        python main_variable.py [options]
#
#    Licensed under BSD 2-Clause License
#
#########################################################
"""

import sys
sys.path.insert(0, './src')

#
# Import third-party packages
#

import numpy as np # For efficient array operations
#import matplotlib.pyplot as plt # For plotting
import time # For timing parts of the script, optimizing run time
#import pandas as pd # Easier data handling
import os # To create directories
import datetime # To access the current time
import sys # To access command line arguments
#import warnings # To warn the user
import csv # For file operations

try: # Seaborn makes prettier plots, but is not installed in a fresh Anaconda python
    import seaborn as sns 
    have_seaborn = True
except ImportError:
    have_seaborn = False

#
# Import other parts of the project
#

from animal import Animal
from population import Population
from environment import Environment
from constants import model_constants
from iterate_population import iterate_population



if __name__ == '__main__':
    # Get model constants
    constants = model_constants
    population_size=constants["environment_sizes"]
    
    if constants["path"] !="":
        path=constants["path"]
    else:  
        p=open("./path.txt","r")
        path=p.readline()
        print(path)
    f_mean = path + "pop{0}_mean_genes.csv".format(constants["use_pop"]) #constants["mean_file"]   #specify path to csv file in constants (from run with constant pop size) 
    f_std = path+ "pop{0}_std_genes.csv".format(constants["use_pop"]) #constants["std_file"] 
    

    # create output directory
    now = datetime.datetime.today()


    start = time.clock()

    # read the csv files

    with open(f_mean) as f:
        reader = csv.reader(f,delimiter=",")
        for (i,row) in enumerate(reader):
            if row:
                if row[0]=="n":
                    data = np.genfromtxt(f_mean,skip_header=i+1,delimiter=",") #reads mean genes and n , nperPos from csv file
                    break
                elif (row[0][0]!="R"): #why 2* [0]?
                    environment=list(map(float,row))  #reads environment values 
    
    std = np.genfromtxt(f_std,skip_header=i+1,delimiter=",")
    if len(std[0])==10:
        std=np.delete(std,3,1) #delete mismatch
    if len(data[0])==11:
        data=np.delete(data,(3,10),1) #delete mismatch,lineage
        std=np.delete(std,(3,10),1) #delete mismatch,lineage
    elif len(data[0])==12: #for files with lineage as last entry    
        data=np.delete(data,(1,8,9),1) #delete environment,m ,ma
        std=np.delete(std,(1,8,9),1) #delete environment,m ,ma
    elif len(data[0])==13:
        data=np.delete(data,(1,4,9,10),1) #delete environment,mismatch, m,ma
        std=np.delete(std,(1,4,9,10),1) #delete environment,mismatch, m,ma
    else:
        raise Exception("Check the delete mechanisms for initial data!")
    
    mean_genes = data[-1,1:-1]  #just last row, that is last generation. All genes, cut off generation and size 
    final_t = data[-1,0]*constants["L"]#*env[0][0]/constants["environments"][0][0] #final time: #generations*lifetime*R_file/R_const  #why R?
    std_genes = data[-1,1:-1] #std of genes from last generation
    std_genes = np.fabs(std_genes) #gives absolute value  
    size=data[-1,-1]    


    # create environments and output information about them

    if constants["trans"]: #for transition runs, change environment parameters in constants file
        env = Environment(*constants["environments"]) 
    else: #or use environment from constant run
        env = Environment(*environment)

    path = "./output_variable/{0:%y}-{0:%m}-{0:%d}_{0:%H}-{0:%M}-{0:%S}-{1}/".format(now,constants["desc"])
   
    
    try: 
        os.makedirs(path)
        os.makedirs(path+"timeseries/")
    except OSError:
        if not os.path.isdir(path):
            raise
    f3 = open(path+"__overview.txt",'w')
    f3.write("initial conditions \n")
    f3.write("R,P,A,B,O\n{0},{1},{2},{3},{5}\n".format(env.R,env.P,env.A,env.B,i,env.O))

    f3.write("Mean genes(I0,I0p,a,b,bp,h,m,ma,s):\n{0}\n".format(mean_genes))
    f3.write("Std genes:\n{0}\n".format(std_genes))
    for key in ['generations','L','kd','ka','tau','q','mu','environments','environment_names','environment_sizes','km','limit','populations','plot_every','verbose',\
'random_choice','std_min','lineage_stop','desc','trans','path','use_pop']:
        f3.write("{0}:\t{1}\n".format(key,constants[key]))
    


    end = time.clock()
    if constants["verbose"]:
        print("Set-up time: {0:.2e}s\n".format(end-start))
    start = time.clock()

    # main loop over multiple populations
    survival_rate = 0
    for k in range(constants["populations"]):
        # write starting genes in files

        f1 = open(path+"pop"+str(k+1)+"_mean_genes.csv",'w')
        f1.write("\nn,I0,I0p,mismatch,a,b,bp,h,s,size,lin\n")        
        f2 = open(path+"pop"+str(k+1)+"_std_genes.csv",'w')
        f2.write("\nn,I0,I0p,mismatch,a,b,bp,h,s,size,lin\n")
            
        # create animals with the mean genes that shall be tested for each environment
        animals=[]
        genes = []
        # unfortunately, the genes are written in a different order as it is used here
        gene_order = [5,6,2,0,1,3,4]
        for j in gene_order:
            if (std_genes[j] > 0):
                genes.append(np.random.normal(size=size,loc=mean_genes[j],scale=std_genes[j])) #if std>0 create size*genes using normal distribution around mean for each environment and gene
            else:
                genes.append(mean_genes[j]*np.ones(size))
        animals.append([Animal(np.array(g),lineage=k) for k,g in enumerate(zip(*genes))])# create animals with genes in environment i
        animals = [item for sublist in animals for item in sublist] # flatten animal list 
        # create a population of population_size animals that have the correct mean genes
        population = Population(population_size,animals)
        
        
        f3.write("Population {0}".format(k+1))
        pop_mean, pop_std, final_gen = iterate_population(k,population,env,f1,f2,path,final_t,True) #start at final time of constant pop run
        end = time.clock()


        if pop_mean is None:
            if not constants["verbose"]:
                print("\n\tDied out! Total time: {0:.2f} min\n".format((end-start)/60)) 
            f3.write(" died at generation {0}!\n".format(final_gen))
        else:
            if not constants["verbose"]:
                print("\n\tSurvived! Total time: {0:.2f} min\n".format((end-start)/60)) 
            survival_rate = survival_rate+1
            f3.write(" survived!\n")

        if constants["verbose"]:
            print("\n---------------------------------------")
            print(" Population {0} done! Total time: {1:.2f} min".format(k+1,(end-start)/60))
            print("---------------------------------------\n")

    f3.write("\n\nIn total, {0}/{1} Populations survived.".format(survival_rate,constants["populations"]))
    f3.close()
