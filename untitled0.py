# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 11:18:49 2016

@author: matthias
"""

 
import sys
sys.path.insert(0, './src')

#
# Import third-party packages
#

import numpy as np # For efficient array operations
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


constants=model_constants

#print(constants["environment_sizes"][2])
#for j in range(3):
#    print(j)
#a=range(10)
#print(sum(a))    
#print(max(constants["environment_sizes"]))
#
#
#
#payoff_factor=np.array([1,1,1,1,1,1,1,1])
#offspring = np.random.poisson(lam=payoff_factor) 
#print(offspring)

#population_size=sum(constants["environment_sizes"])
animal_list=[]
nE=len(constants["environments"])
Sizes=[3,2,5,4]
pop_size=sum(constants["environment_sizes"])
for j in range(nE):               
    # create animals in each environment according to environment_sizes that already have the correct random genes
    animal_list.extend([Animal(np.array([]),j) for _ in range(constants["environment_sizes"][j])])
                   
# create a Population from animal_list
for animal in animal_list:
    print(animal.position)

pop = Population(pop_size,animal_list)
#pop._positions=pop.positions()
#print(list(pop._positions))

#for animal in population._animals:
#    print(animal.position)

#print(animal_list)
#print("\n")
#runner_min1=0
#runner_max1=0
#runner_min2=0
#runner_max2=0
#chosen=[]
#insert=0 #insertion index for clones
#
#for j in range(nE):
#    pop._animals=animal_list #refresh animal list
#    pop._positions=pop.positions() #refresh positions
#    print(list(pop._positions))
#    clones=[] #set empty every iteration
#    runner_max+=pop._positions[j]
#    print("runner_max:"+str(runner_max))
#    #chosen.extend(np.random.choice(animal_list[runner_min:runner_max], constants["environment_sizes"][j],replace=False))
#    clones.extend(np.random.choice(animal_list[runner_min:runner_max],constants["environment_sizes"][j], replace=False))
#    print(clones)
#    insert+=pop._positions[j] #index for inserting clones at the end of each environment (for first inserted animal)
#    print("insert:"+str(insert))
#    for i in range(len(clones)):
#        animal_list.insert(insert,clones[i])#insert clones consecutively at the end of each environmen
#        insert+=1 #increase insert position by one
#    #clones = [Animal(x.genes,x.position) for x in clone_candidates]
#    
#    for animal in animal_list:
#        print(animal.position)
#    pop._positions=pop.positions()
#    print("\n")
#    runner_min+=pop._positions[j]
#    runner_max+=len(clones)
#    print("clones:"+str(len(clones)))
#    print("runner_min:"+str(runner_min))
#    print("insert:"+str(insert))
#    print(animal_list)
#    
     
#   
#for j in range(nE):
#    pop._animals=animal_list #refresh animal list
#    pop._positions=pop.positions() #refresh positions
#    print(list(pop._positions))   
#    runner_max1 += pop._positions[j]
#    runner_max2 += pop._positions[j]
#    print("runner_min1:"+str(runner_min1))
#    print("runner_min2:"+str(runner_min2))
#    print("runner_max1:"+str(runner_max1))
#    print("runner_max2:"+str(runner_max2))
#    if pop._positions[j] > pop._constants["environment_sizes"][j]:          
#        animal_list[runner_min1:runner_max1]=np.random.choice(animal_list[runner_min1:runner_max1], constants["environment_sizes"][j],replace=False)           
#        N_rem=pop._positions[j]-pop._constants["environment_sizes"][j]      #number of removed animals
#        runner_max1 -= N_rem
#        runner_max2 -= N_rem #adjust upper boundary         
#    elif pop._positions[j] < pop._constants["environment_sizes"][j]:
#        clones=[] #set empty every iteration
#        N_cloned=pop._constants["environment_sizes"][j]-pop._positions[j] #number of needed clones
#       # print(N_cloned)        
#        clones.extend(np.random.choice(animal_list[runner_min2:runner_max2],N_cloned, replace=True))
#       # print(clones)
#        insert+=pop._positions[j] #index for inserting clones at the end of each environment (for first inserted animal)
#        #print("insert:"+str(insert))
#        for i in range(len(clones)):
#            animal_list.insert(insert,clones[i])#insert clones consecutively at the end of each environmen
#            insert+=1 #increase insert position by one     
#        runner_max1 += len(clones)
#        runner_max2 += len(clones)
##        print("clones:"+str(len(clones)))
##        print("runner_min:"+str(runner_min))
##        print("insert:"+str(insert))
##        print(animal_list)
#    pop._animals=animal_list
#    pop._positions=pop.positions()
#    runner_min1 += pop._positions[j]  
#    runner_min2 += pop._positions[j]      
##    print("runner_min1:"+str(runner_min1))
#       
#for animal in animal_list:
#            print(animal.position)    
#print(animal_list)

calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(pop._positions))
lifetime_payoff = calc_payoff(pop._animals)
mean_payoff     = np.mean(lifetime_payoff)
population_size=sum(pop._constants["environment_sizes"])
if (mean_payoff == 0):
    raise RuntimeError("Mean payoff of population decreased to 0. Check your parameters!")
else:
    payoff_factor = lifetime_payoff/mean_payoff
print("animal_list:")
print(np.array(animal_list))
offspring = np.random.poisson(lam=payoff_factor) #number of offspring drawn from Poisson distr for each animal 
print("offspring")
print(offspring)
born_animals = np.repeat(pop._animals,offspring) # Create list with offspring repeats for each animal (cloned animals)
print("born_animals:")
print(born_animals)
mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.position)) #animals are created with mutated genes and at the position of their parents
new_animals = mutate_pop(born_animals) #create and mutate offspring (use mutated genes as parent genes), ordered with respect to environment
print("mutated animals:")
print(new_animals)
pop._animals = new_animals 
pop._positions = pop.positions() #update positions

nE=len(pop._constants["environments"])
s1=[]
for animal in animal_list:
    s1.append(animal.position)
  
if pop._constants["verbose"]:
    print("\n\nAnimals per environment: {0}".format(pop._positions))
    print("Population size: {0}\tMean payoff: {1:.2f}".format(population_size,mean_payoff))
runner_min1=0 #lower boundary for killing animals in new_animals
runner_max1=0 #upper boundary for killing animals in new_animals
runner_min2=0 #lower boundary for cloning animals in new_animals
runner_max2=0 #upper boundary for cloning animals in new_animals
insert=0 #defines where in the list clones shall be inserted
for j in range(nE):
    pop._animals=new_animals #refresh animal list
    pop._positions=pop.positions() #refresh positions   
    runner_max1 += pop._positions[j] #increase upper boundary for modification in new_animals
    runner_max2 += pop._positions[j]
    if pop._positions[j] > pop._constants["environment_sizes"][j]:          #if environment overcrowded randomly select the correct amount of animals out of the new generation
        new_animals=list(new_animals)
        print("new_animals before removing:")
        print(np.array(new_animals))
        new_animals[runner_min1:runner_max1]=np.random.choice(new_animals[runner_min1:runner_max1], pop._constants["environment_sizes"][j],replace=False)         
        new_animals=np.array(new_animals)   
        print("new_animals after removing:")
        print(new_animals)
        for animal in new_animals:
            print(animal.position)
        N_rem=pop._positions[j]-pop._constants["environment_sizes"][j]      #number of removed animals
        runner_max1 -= N_rem #adjust upper boundary 
        runner_max2 -= N_rem        
        insert += pop._constants["environment_sizes"][j]
        print("insert:")
        print(insert)
    elif pop._positions[j] < pop._constants["environment_sizes"][j]: #if too few, create the needed amount of clones of randomly chosen animals in new_animals in the right environment
        clones=[] #list for cloned animals for each environment, set empty every iteration
        N_cloned=pop._constants["environment_sizes"][j]-pop._positions[j] #number of needed clones  
        if pop._positions[j]!=0:
            clones.extend(np.random.choice(new_animals[runner_min2:runner_max2],N_cloned, replace=True))
        else: #if animals in environment dont have any offspring, create new random Animals
            clones.extend([Animal(np.array([]),j) for _ in range(N_cloned)])
        print("Number of animals to be cloned:")
        print(N_cloned)
        print("clones:")
        print(clones)
        print("new_animals before cloning:")
        print(new_animals)
        insert+=pop._positions[j] #index for inserting clones at the end of each environment (for first inserted animal)
        print("insert:")
        print(insert)
        for i in range(len(clones)):
            new_animals=np.insert(new_animals,insert,clones[i])#insert clones consecutively at the end of each environmen
            insert+=1 #increase insert position by one     
        runner_max1 += len(clones) #adjust upper boundary 
        runner_max2 += len(clones)
        print("new_animals after cloning:")
        print(new_animals)
        for animal in new_animals:
            print(animal.position)
    else:
        insert += pop._constants["environment_sizes"][j]       
    pop._animals=new_animals
    pop._positions=pop.positions()
    runner_min1 += pop._positions[j]   #adjust lower boundary 
    runner_min2 += pop._positions[j]  
    



print(new_animals)
s2=[]
for animal in new_animals:
    s2.append(animal.position)
if s1!=s2:
    raise













