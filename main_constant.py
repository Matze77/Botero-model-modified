#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#    main_constant.py
#    Author: Dion Häfner (dionhaefner@web.de)
#    
#    Main controller, with constant population size
#
#    Usage:
#        python main_constant.py [options]
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
import math
import matplotlib.pyplot as plt # For plotting
import time # For timing parts of the script, optimizing run time
import pandas as pd # Easier data handling
import os # To create directories
import datetime # To access the current time
import sys # To access command line arguments
import warnings # To warn the user

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


path=str()
if __name__ == '__main__':
    # Get model constants
    constants = model_constants
    population_size=constants["environment_sizes"]
    

    if have_seaborn: # initialize seaborn
        sns.set('poster')
        sns.set_palette("deep", desat=.6)
        sns.set_context(rc={"figure.figsize": (10, 7.5)})

    # create output directory
    now = datetime.datetime.today()
    path = "./output/{0:%y}-{0:%m}-{0:%d}_{0:%H}-{0:%M}-{0:%S}-{1}/".format(now,constants["desc"])
    try: 
        os.makedirs(path)
        os.makedirs(path+"timeseries/")
    except OSError:
        if not os.path.isdir(path):
            raise
    #save path to use it in main_variable
    p=open("./path.txt","w")
    p.write(path)
    p.close()
    # write simulation parameters
    f = open(path+"parameters.txt","w")
    for key in constants:
        f.write("{0}:\t{1}\n".format(key,constants[key]))
    f.close()    

    # plot environments
       
    print(constants["environments"])
  
    T=round(5*constants["environments"][0]*constants["L"])  #how many time steps to plot
    step=(constants["environments"][0]*constants["L"])/100 #step size
    M=math.floor(T/step+1) #number of data points
    t0 = np.arange(0,T,step)
    env = Environment(*constants["environments"]) #create new environment

    env_val = np.array(list(map(env.evaluate,t0))) #calculate its values
    
      
    plt.figure()
    plt.plot(t0,env_val[0:M,0],label='E', linewidth=0.5)
    plt.plot(t0,env_val[0:M,1],'.',label='C', markersize=3.0)
    plt.legend()
    plt.xlabel('Time t')
    plt.ylabel('E, C')
    plt.ylim(-2,2)
    plt.savefig(path+'environment.pdf',bbox_inches='tight')


    # main loop over multiple populations
    means, stds = [], []
    error_occured = False
    for k in range(constants["populations"]):
        start = time.clock()
       
        # in case a population dies out, it is repeated
        repeat = True
        animal_list=[]
        while repeat:             
                # create animals in each environment according to environment_sizes that already have the correct random genes
            animal_list.extend([Animal(np.array([]),lineage=_) for _ in range(constants["environment_sizes"])]) #if only one value is given
            # create a Population from animal_list
#            for a in animal_list:
#                print(a.genes)
            population = Population(population_size,animal_list)

            end = time.clock()
            if constants["verbose"]:
                print("Set-up time: {0:.2e}s\n".format(end-start))
            start = time.clock()

            # initial output
            f1 = open(path+"pop"+str(k+1)+"_mean_genes.csv",'w')
            f2 = open(path+"pop"+str(k+1)+"_std_genes.csv",'w')            
            f1.write("R,P,A,B,O\n{0},{1},{2},{3},{4}\n".format(env.R,env.P,env.A,env.B,env.O))
            f2.write("R,P,A,B,O\n{0},{1},{2},{3},{4}\n".format(env.R,env.P,env.A,env.B,env.O))

            f1.write("\nn,I0,I0p,mismatch,a,b,bp,h,m,ma,s,nperPos,lin\n")
            f2.write("\nn,I0,I0p,,mismatch,a,b,bp,h,m,ma,s,nperPos,lin\n")
                    
            # iterate on the population and create outputs
            try:
                #%timeit iterate_population(k,population,environments,f1,f2,path)
                pop_mean, pop_std, _ = iterate_population(k,population,env,f1,f2,path) #create plots and return values for last generation
                repeat = False #if pop dies out no error occurs?
            except RuntimeError:
                error_occured = True
                pass


        end = time.clock()
        if constants["verbose"]:
            print("\n---------------------------------------")
            print(" Population {0} done! Total time: {1:.4f} min".format(k+1,(end-start)/60))
            print("---------------------------------------\n")
        else:
            print("\n\tDone! Total time: {0:.4f} min\n".format((end-start)/60))

        plt.close('all')

        means.append(pop_mean)
        stds.append(pop_std)

    # plot average genes of ALL populations run (always last generation)

    plt.figure()
    average = pd.concat(means)
    if have_seaborn:
        sns.violinplot(data=average,scale='width')
    else:
        average.boxplot()
    plt.ylim(-2,2)
    plt.xlabel("Genes")
    plt.ylabel("Average values")
    plt.savefig(path+"total_average_env_.pdf",bbox_inches='tight')
    plt.close()

    if error_occured:
        warnings.warn("At least one population died out and was repeated!")


