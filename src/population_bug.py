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

def mismatch(animal_list):
    '''Returns the mismatch of each animal'''
    fun = np.vectorize(lambda x: x.mismatch)
    mis = fun(animal_list)
    return np.array(mis,ndmin=1) 
    
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
            raise TypeError('First argument must be of type int, second of type list of Animals.')

    def animals(self):
        """Returns the ndarray of animals"""
        return self._animals

    def size(self):
        """Returns the current size of the population"""
        return self._size



    def react(self,E,C,evolve_all=False):
        """Calculates the insulation of each Animal in the Population based on cue C and environment E"""                              
        for animal in self._animals:           
            animal.react(E,C,evolve_all)#animal reacts to environment, possibly migrates                       
    def breed_constant(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with CONSTANT population size"""
        nE=len(self._constants["environments"])
        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(self._positions))
        lifetime_payoff = calc_payoff(self._animals)
        mean_payoff     = np.mean(lifetime_payoff)
        population_size=sum(self._constants["environment_sizes"])
        if (mean_payoff == 0):
            raise RuntimeError("Mean payoff of population decreased to 0. Check your parameters!")
        else:
            payoff_factor = lifetime_payoff/mean_payoff
        offspring = np.random.poisson(lam=payoff_factor) #number of offspring drawn from Poisson distr for each animal 
#        pos=np.vectorize(lambda x: x.position)
#        positions=pos(self._animals)
        d=np.sum(offspring)-self._constants["environment_sizes"]
        if d>0:
            for i in range(d):
                pf=payoff_factor
                pf[offspring==0]=np.max(pf)+1
                m=np.argmin(pf)
                offspring[m]-=1               
                                                        
        born_animals = np.repeat(self._animals,offspring) # Create list with offspring repeats for each animal (cloned animals)      
        mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.position,x.lineage)) #animals are created with mutated genes and at the position of their parents
        new_animals = mutate_pop(born_animals) #create and mutate offspring (use mutated genes as parent genes), ordered with respect to environment        
        self._animals = new_animals 
        self._positions = self.positions() #update positions

  
        if self._constants["verbose"]:
            print("\n\nAnimals per environment: {0}".format(self._positions))
            print("Population size: {0}\tMean payoff: {1:.2f}".format(population_size,mean_payoff))
            
        animals={} #dictionary for animals in different environments
        for j in range(nE):
            a=[]
            for animal in self._animals:

                if animal.position==j:
                    a.append(animal)
            animals[j]=a     
        for j in range(nE):        
            d=self._positions[j] - self._constants["environment_sizes"][j]
            if d>0:          #if environment overcrowded randomly select the correct amount of animals out of the new generation

#                for i in range(d):
#                    m=np.argmax(mismatch(animals[j]))
#                    animals[j]=np.delete(animals[j],m)

                animals[j]=np.random.choice(animals[j], self._constants["environment_sizes"][j],replace=False)         
            elif d<0: #if too few, create the needed amount of clones of randomly chosen animals in new_animals in the right environment
                N_cloned=self._constants["environment_sizes"][j]-self._positions[j] #number of needed clones                  
                if self._positions[j]!=0:
                    animals[j].extend(np.random.choice(animals[j],N_cloned, replace=True))
                else: #if animals in environment dont have any offspring, create new random Animals
                    animals[j].extend([Animal(np.array([]),j) for _ in range(N_cloned)])

        a=[]
        for j in range(nE):
                a.extend(animals[j]) #create list of animals o

        self._animals=np.array(a)
        self._positions=self.positions()   

    def breed_variable(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with VARIABLE population size"""
        nE = len(self._constants["environments"])
        if len(self._constants["q"])>1:
            q=self._constants["q"][self.position]
        else:
            q=self._constants["q"][0]

        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(self._positions))
        lifetime_payoff = calc_payoff(self._animals)
        payoff_factor=np.array([])
        for j,animal in enumerate(self._animals):
             payoff_factor=np.append(payoff_factor,q*lifetime_payoff[j])
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
            
        animals={} #dictionary for animals in different environments
        for j in range(nE):
            a=[]
            for animal in self._animals:
                if animal.position==j:
                    a.append(animal)
            animals[j]=a               
        for j in range(nE):            
            if self._positions[j] > self._constants["environment_sizes"][j]:          #if environment overcrowded randomly select the correct amount of animals out of the new generation
                animals[j]=np.random.choice(animals[j], self._constants["environment_sizes"][j],replace=False)         
       
        a=[]
        for j in range(nE):
                a.extend(animals[j]) #create list of animals to update self._animals           
        self._animals=np.array(a)
        self._positions=self.positions()    
        self._size = len(new_animals)

    def positions(self):
        """Returns the number of animals in each environment"""
        fun = np.vectorize(lambda x: x.position)
        pos = fun(self._animals)
        pos=np.array(pos,ndmin=1) #important if just one animal left
        return np.bincount(pos,minlength=len(self._constants["environments"]))

    def lineage(self):
        '''Returns the lineage of each animal'''
        fun = np.vectorize(lambda x: x.lineage)
        lin = fun(self._animals)
        return np.array(lin,ndmin=1) 

