#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#   iterate_population.py
#   Author: Dion Häfner (dionhaefner@web.de)
#   
#   Main controller for a given population
#   
#   Licensed under BSD 2-Clause License
#
#########################################################
"""

# Import third party libraries
import numpy as np 
import time
import sys
#from numba import jit

# Import other parts of the project
#from animal import Animal
#from population import Population
#from environment import Environment
from constants import model_constants
from output_population import output_population, plot_size

#@jit
def iterate_population(k,population,environments,f1,f2,path,t=0,variable=False):
    """ 
    MAIN CONTROLLER
    Inputs:
        k: population counter,  population: the Population instance to be iterated,
        environments: Environment instances to be operated on,
        f1: pointer to output file for gene means,  f2: for gene standard deviations,
        path: path to the output files  t: initial time,   
        variable: variable population size
    """

    constants = model_constants
    nE = len(environments)

    for j in np.arange(constants["generations"]):  
        # MAIN TIME STEP LOOP
     #   start = time.clock()
        t1=time.clock()   
        mean,std=output_population(population,f1,f2,j,k,path,False,t,environments,variable) #creates plots and csv files
     #   t2=time.clock()
      #  print("plot Time: {0:.2e}s\n".format(t2-t1))
        t1=time.clock() 
        for _ in range(constants["L"]): #loop for time steps in each animal's life
            E, C = np.empty(nE), np.empty(nE) #initialze E, C

            for (i,env) in enumerate(environments): #calculate E,C at time t for all environments
                E[i], C[i] = env.evaluate(t)

            population.react(E,C) #animals react to environment
            t = t+1
       # t2=time.clock()
       # print("react Time (lifetime): {0:.2e}s\n".format(t2-t1))
        t1=time.clock() 
        if variable:
            population.breed_variable() #old generation is replaced by new one
        else:
            population.breed_constant()

        if population.size() == 0:
            print("Population died out!\n\n")
            return None, None, j
        #t2=time.clock()
        #print("breedTime: {0:.2e}s\n".format(t2-t1))        
        #t1=time.clock()
        population.react(E,C,1)   #all plastic animals react
        #t2=time.clock()
        #print("react Time: {0:.2e}s\n".format(t2-t1))
        end = time.clock()
        if constants["verbose"]:
            print("Computation time: {0:.2e}s".format(end-start))

	    # Print progress bar
        percent = float(j+1) / constants["generations"]
        hashes = '#' * int(round(percent * 20))
        spaces = ' ' * (20 - len(hashes))
        sys.stdout.write("\rProgress population {2} of {3}: [{0}] {1:.1f}%".format(hashes + spaces, percent * 100,k+1,constants["populations"]))
        sys.stdout.flush()     
        stop=False
        std_min=constants["std_min"]
        if len(std_min)!=0 and std_min:
            stop=True
            for i in range(nE):
                for l in ["I0","I0p","a","b","bp","h","s"]:
                    try:
                        if float(std[i][l])>std_min[i]:   
                            stop=False
                    except:
                        if float(std[i][l])>std_min[0]:
                            stop=False                              
        if max(np.bincount(population.lineage()))==len(population._animals):
            print("\n Common ancestry reached, loop stopped after {0} generations!".format(j))
            break
        elif stop: #if all std above are <std_min: break loop
            print("\n Desired std reached, loop stopped after {0} generations!".format(j))
            break



    # Final outputs for each population
    final_mean, final_std = output_population(population,f1,f2,j,k,path,True,t,environments,variable) #force last plot
    f1.close()
    f2.close()
    if  variable:
        plot_size(path,path+"pop"+str(k+1)+"_mean_genes.csv",k) #plots the number of animals in each environment for each generation

    return final_mean, final_std, j #j: last generation
