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

    if float(constants["mutation"][1])>0: #if mu is mutable, also this trait can be chosen for hgt
        if float(constants["mutation"][3])>0: #if sc is mutable                          
            g=np.arange(0,9) #genes to choose from
        else:
            g=np.arange(0,8) 
    elif float(constants["mutation"][3])>0: 
        g=np.delete(np.arange(0,9),7)
    else:
        g=np.arange(0,8)
        
    E, C = environment.evaluate(t)
    population.react(E,C,g,1)   #for initial random animals: all plastic animals react
    output_population(population,f1,f2,0,k,path,True,t,environment,sizes,variable) #creates plots and csv files
    for j in np.arange(1,constants["generations"]+1):  
        # MAIN TIME STEP LOOP
        start = time.clock()    
        for _ in range(constants["L"]): #loop for time steps in each animal's life        
            t = t+1
            E, C = environment.evaluate(t) #calculate E,C at time t for all environment
            population.react(E,C,g) #animals react to environment

          
        mean,std=output_population(population,f1,f2,j,k,path,False,t,environment,sizes,variable) #creates plots and csv files
                
        if constants["save_all"]:
            f3 = open(path+"all_genes/pop{0}_gen{1}.csv".format(k,j),'w')
            f3.write("h,s,a,I0,I0p,b,bp,mu,sc,t,M,W,I,A,T,lin\n")
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
            population.breed_variable(j) #old generation is replaced by new one
        else:
            population.breed_constant(j)
                
        sizes.append(population._size)
        
        if (float(population.size())/constants["size"] <= float(constants["stop_below"][0]) and j>float(constants["stop_below"][1])) or population.size()==0:
            break
        population.react(E,C,g,1)   #all plastic animals react
        end = time.clock()
        if constants["verbose"]:
            print("Computation time: {0:.2e}s".format(end-start))

	     #Print progress bar
        percent = float(j) / constants["generations"]
        hashes = '#' * int(round(percent * 20))
        spaces = ' ' * (20 - len(hashes))
        sys.stdout.write("\rProgress population {2} of {3}: [{0}] {1:.1f}%".format(hashes + spaces, percent * 100,k+1,constants["populations"]))
        sys.stdout.flush()               
              
        stop=False
        std_min=constants["std_min"][1]
        if j>constants["std_min"][0] and std_min>0:
            stop=True
            for c,l in enumerate(["I0","I0p","W","a","b","bp","h","mu","s","sc","t"]):
                if c not in constants["std_min"][2:] and len(constants["std_min"])>2:
                    continue
                #print(std[l])
                if float(std[l])>std_min:   
                    stop=False
                                  
        if constants["lineage_stop"] and max(np.bincount(population.lineage()))==len(population._animals):
            print("\n Common ancestry reached, loop stopped after {0} generations!".format(j))
            break                          
        if stop: #if all std above are <std_min: break loop
            print("\n Desired std reached, loop stopped after {0} generations!".format(j))
            break
        

    # Final outputs for each population
    final_mean, final_std = output_population(population,f1,f2,j,k,path,True,t,environment,sizes,variable) #force last plot
    f1.close()
    f2.close()
    if  variable:
        plot_size(path,path+"pop"+str(k+1)+"_mean_genes.csv",k) #plots the number of animals in each environment for each generation
        
    return final_mean, final_std, j #j: last generation
