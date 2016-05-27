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
from constants import model_constants
from output_population import output_population, plot_size

#@jit
def iterate_population(k,population,environment,f1,f2,path,t=0,variable=False):
    """ 
    MAIN CONTROLLER
    Inputs:
        k: population counter,  population: the Population instance to be iterated,
        environment: Environment instance to be operated on,
        f1: pointer to output file for gene means,  f2: for gene standard deviations,
        path: path to the output files  t: initial time,   
        variable: variable population size
    """

    constants = model_constants
    sizes=[population._size]
    if variable:
        if constants["start_hgt"] and constants["hgt"]:
            for animal in population._animals:
                genes=animal.genes                    
                genes[8]=0.75                    
                animal.genes=genes
            
    E, C = environment.evaluate(t)
    population.react(E,C,1)   #for initial random animals: all plastic animals react
    output_population(population,f1,f2,0,k,path,True,t,environment,sizes,variable) #creates plots and csv files
    for j in np.arange(1,constants["generations"]+1):  
        # MAIN TIME STEP LOOP
        start = time.clock()    
        for _ in range(constants["L"]): #loop for time steps in each animal's life        
            t = t+1
            E, C = environment.evaluate(t) #calculate E,C at time t for all environment
            population.react(E,C) #animals react to environment

          
        output_population(population,f1,f2,j,k,path,False,t,environment,sizes,variable) #creates plots and csv files
                
        if constants["save_all"]:
            f3 = open(path+"all_genes/pop{0}_gen{1}.csv".format(k,j),'w')
            f3.write("h,s,a,I0,I0p,b,bp,mu,t,ta,M,W,I,A,T,lin\n")
            for a in population.animals():
                for i,g in enumerate(a.genes):
                    if i!=0:
                        f3.write(",") 
                    f3.write(str(g))
                f3.write(","+str(a.mismatch))
                f3.write(","+str(a.lifetime_payoff()))
                f3.write(","+str(a.insulation))
                f3.write(","+str(a.adjustments))
                f3.write(","+str(a.transfers))
                f3.write(","+str(a.lineage))            
                f3.write("\n") 
            f3.close()
        
        if variable:
            population.breed_variable() #old generation is replaced by new one
        else:
            population.breed_constant()
                
        sizes.append(population._size)
        
        if float(population.size())/constants["size"] <= constants["stop_below"]:
            break
        population.react(E,C,1)   #all plastic animals react
        end = time.clock()
        if constants["verbose"]:
            print("Computation time: {0:.2e}s".format(end-start))

	     #Print progress bar
        percent = float(j) / constants["generations"]
        hashes = '#' * int(round(percent * 20))
        spaces = ' ' * (20 - len(hashes))
        sys.stdout.write("\rProgress population {2} of {3}: [{0}] {1:.1f}%".format(hashes + spaces, percent * 100,k+1,constants["populations"]))
        sys.stdout.flush()                             



    # Final outputs for each population
    final_mean, final_std = output_population(population,f1,f2,j,k,path,True,t,environment,sizes,variable) #force last plot
    f1.close()
    f2.close()
    if  variable:
        plot_size(path,path+"pop"+str(k+1)+"_mean_genes.csv",k) #plots the number of animals in each environment for each generation
        
    return final_mean, final_std, j #j: last generation
