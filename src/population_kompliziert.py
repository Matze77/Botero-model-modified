#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#    population.py
#    Author: Dion Häfner (dionhaefner@web.de)
#    
#    Implements Population class
#    
#    Licensed under BSD 2-Clause License
#
#########################################################
"""

import numbers
import numpy as np
import time

from constants import model_constants # Import model constants
from animal import Animal


class Population:
    def __init__(self,size,animals):
        """Takes a population size and a list of Animal as input"""
        if (isinstance(size,int) & (all(isinstance(x,Animal) for x in animals))):
            if (size == len(animals)):
                self._animals     = np.array(animals)
                self._size = size
                self._constants    = model_constants
                self._positions = self.positions()
            else:
                raise ValueError('The size parameter must be equal to the length of the list of animals.')
        else:
            raise TypeError('First argument must be of type int, second of type list of Animal.')

    def animals(self):
        """Returns the ndarray of animals"""
        return self._animals

    def size(self):
        """Returns the current size of the population"""
        return self._size

    def react(self,E,C,evolve_all=False):
        """Calculates the insulation of each Animal in the Population based on cue C and environment E"""                              
        migrated=[] 
        not_migrated=[]
        for animal in self._animals:           
            p1=animal.position #position before reaction
            animal.react(E,C,evolve_all)#animal reacts to environment, possibly migrates
            p2=animal.position #position after reaction
            if p1!=p2:        # that means, the animal is migrated                        
                migrated.append(animal) #save migrated animal in list                              
            else:
                not_migrated.append(animal) 
                
        self._animals=np.array(not_migrated,ndmin=1)
        if len(self._animals)==0:  #to prevent empty array
            self._animals=np.array(migrated[0],ndmin=1)
            del migrated[0]
        self._positions=self.positions()        
        j=0    #counter
        for mig in migrated:    
            insert=0  #index to insert migrated animal in different environment
            for p in range(mig.position+1):
                insert+=self._positions[p]
            self._animals=np.insert(self._animals,insert,mig) #insert animal in new environment
            self._positions=self.positions() #update positions
            j+=1 
    def breed_constant(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with CONSTANT population size"""
        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(self._positions))
        lifetime_payoff = calc_payoff(self._animals)
       # print(list(lifetime_payoff))
        if lifetime_payoff.any()<0:
            raise Exception("!!!")
        mean_payoff     = np.mean(lifetime_payoff)
        population_size=sum(self._constants["environment_sizes"])
        if (mean_payoff == 0):
            raise RuntimeError("Mean payoff of population decreased to 0. Check your parameters!")
        else:
            payoff_factor = lifetime_payoff/mean_payoff
        offspring = np.random.poisson(lam=payoff_factor) #number of offspring drawn from Poisson distr for each animal 
        born_animals = np.repeat(self._animals,offspring) # Create list with offspring repeats for each animal (cloned animals)
      
        mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.position)) #animals are created with mutated genes and at the position of their parents
        new_animals = mutate_pop(born_animals) #create and mutate offspring (use mutated genes as parent genes), ordered with respect to environment
             
        
        self._animals = new_animals 
        self._positions = self.positions() #update positions
        nE=len(self._constants["environments"])
  
        if self._constants["verbose"]:
            print("\n\nAnimals per environment: {0}".format(self._positions))
            print("Population size: {0}\tMean payoff: {1:.2f}".format(population_size,mean_payoff))
        runner_min1=0 #lower boundary for killing animals in new_animals
        runner_max1=0 #upper boundary for killing animals in new_animals
        runner_min2=0 #lower boundary for cloning animals in new_animals
        runner_max2=0 #upper boundary for cloning animals in new_animals
        insert=0 #defines where in the list clones shall be inserted
        for j in range(nE):  
            runner_max1 += self._positions[j] #increase upper boundary for modification in new_animals
            runner_max2 += self._positions[j]
            if self._positions[j] > self._constants["environment_sizes"][j]:          #if environment overcrowded randomly select the correct amount of animals out of the new generation
                new_animals=list(new_animals)
                new_animals[runner_min1:runner_max1]=np.random.choice(new_animals[runner_min1:runner_max1], self._constants["environment_sizes"][j],replace=False)         
                new_animals=np.array(new_animals)
                N_rem=self._positions[j]-self._constants["environment_sizes"][j]      #number of removed animals
                runner_max1 -= N_rem #adjust upper boundary 
                runner_max2 -= N_rem     
                insert += self._constants["environment_sizes"][j]
            elif self._positions[j] < self._constants["environment_sizes"][j]: #if too few, create the needed amount of clones of randomly chosen animals in new_animals in the right environment
                clones=[] #list for cloned animals for each environment, set empty every iteration
                N_cloned=self._constants["environment_sizes"][j]-self._positions[j] #number of needed clones                  
                if self._positions[j]!=0:
                    clones.extend(np.random.choice(new_animals[runner_min2:runner_max2],N_cloned, replace=True))
                else: #if animals in environment dont have any offspring, create new random Animals
                    clones.extend([Animal(np.array([]),j) for _ in range(N_cloned)])
                insert+=self._positions[j] #index for inserting clones at the end of each environment (for first inserted animal)
                for i in range(len(clones)):
                    new_animals=np.insert(new_animals,insert,clones[i])#insert clones consecutively at the end of each environmen
                    insert+=1 #increase insert position by one     
                    
                runner_max1 += len(clones) #adjust upper boundary 
                runner_max2 += len(clones)
            else:
                insert += self._constants["environment_sizes"][j]      
            self._animals=new_animals
            self._positions=self.positions()                                                                
            runner_min1 += self._positions[j]   #adjust lower boundary 
            runner_min2 += self._positions[j]      
    def breed_variable(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with VARIABLE population size"""
        nE = len(self._constants["environments"])
        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(self._positions))
        lifetime_payoff = calc_payoff(self._animals)

        #max_payoff     = 1/self._constants["q"] #(1-1/nE)/self._constants["q"]
        payoff_factor     = self._constants["q"]*lifetime_payoff
        offspring     = np.random.poisson(lam=payoff_factor)
        born_animals     = np.repeat(self._animals,offspring)

        try: # check if all animals are dead yet
            born_animals[0]
        except IndexError:
            self._size = 0
            return

        mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.position))
        new_animals = mutate_pop(born_animals)
        self._animals = new_animals 
        self._positions = self.positions() #update positions
        N = len(new_animals)
        if self._constants["verbose"]:
            print("\n\nAnimals per environment: {0}".format(self._positions))
            print("Population size: {0}".format(N))
        runner_min1=0 #lower boundary for killing animals in new_animals
        runner_max1=0 #upper boundary for killing animals in new_animals
        for j in range(nE):
            runner_max1 += self._positions[j] #increase upper boundary for modification in new_animals
            if self._positions[j] > self._constants["environment_sizes"][j]:          #if environment overcrowded randomly select the correct amount of animals out of the new generation
                new_animals=list(new_animals)
                new_animals[runner_min1:runner_max1]=np.random.choice(new_animals[runner_min1:runner_max1], self._constants["environment_sizes"][j],replace=False)         
                new_animals=np.array(new_animals)
                N_rem=self._positions[j]-self._constants["environment_sizes"][j]      #number of removed animals
                runner_max1 -= N_rem #adjust upper boundary                                
            self._animals=new_animals
            self._positions=self.positions()
            runner_min1 += self._positions[j]   #adjust lower boundary      
        self._size = len(new_animals)

    def positions(self):
        """Returns the number of animals in each environment"""
        fun = np.vectorize(lambda x: x.position)
        pos = fun(self._animals)
        pos=np.array(pos,ndmin=1) #important if just one animal left
        return np.bincount(pos,minlength=len(self._constants["environments"]))

