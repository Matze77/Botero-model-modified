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

#import numbers
import numpy as np
import time

#from numba import jit


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
 #   @jit
    def breed_constant(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with CONSTANT population size"""
        nE=len(self._constants["environments"])
        pos=np.vectorize(lambda x: x.position)
        positions=pos(self._animals)
        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(self._positions))
        lifetime_payoff = calc_payoff(self._animals)
        mean_payoff     = [np.mean(lifetime_payoff[positions==j]) for j in range(nE)]
        population_size=sum(self._constants["environment_sizes"])
        if (mean_payoff == 0):
            raise RuntimeError("Mean payoff of population decreased to 0. Check your parameters!")
        else:
            try:
                payoff_factor = np.array([lifetime_payoff[positions==j]/mean_payoff[j] for j in range(nE)]).flatten()#flatten array
              
                if payoff_factor.shape[0]!=population_size:                   
                    raise                   
            except:
                payoff_factor = np.hstack(np.array([lifetime_payoff[positions==j]/mean_payoff[j] for j in range(nE)]))
                # if array has different dimensions (due to migration) in each row
        offspring = np.random.poisson(lam=payoff_factor) #number of offspring drawn from Poisson distr for each animal 
        for j in range(nE):
            d=np.sum(offspring[positions==j])-self._constants["environment_sizes"][j]
            if self._constants["random_choice"]:    
                Ind=np.arange(0,len(offspring))
                I=Ind[(positions==j) & (offspring>0)] # array of relevant indices to choose from
                if len(I)==0: 
                    I=Ind[positions==j]
                if d>0:#if environment overcrowded reduce offspring randomly
                    add=-1
                else:#if too few, increase offspring randomly
                    add=1
                for i in range(abs(d)):  
                    stop=False
                    while not stop:
                        m=np.random.choice(I)
                        if  offspring[m]!=0 or sum(offspring[positions==j])==0: #animals with zero offspring are disregarded unless there are no animals in the environment
                            offspring[m]+=add
                            stop=True
            else:
                pf=np.array(payoff_factor)
                if d>0:#if environment overcrowded let the least fit animals have less offspring
                    mx=np.max(pf)
                    pf[(offspring==0)|(positions!=j)]=mx+1 #to make sure that no animals with offspring 0 or from other environment are selected
                    for i in range(d):                         
                        m=np.argmin(pf)
                        offspring[m]-=1
                        if offspring[m]==0:
                            pf[m]=mx+1
                        
                elif d<0:#if too few, clone the fittest animals
                    pf[positions!=j]=0
                    for i in range(-d):                     
                        m=np.argmax(pf)
                        pf[m]=0 #to ensure that not all additional offspring is from one animal
                        offspring[m]+=1                                               
        born_animals = np.repeat(self._animals,offspring) # Create list with offspring repeats for each animal (cloned animals)      
        mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.position,x.lineage)) #animals are created with mutated genes and at the position of their parents
        new_animals = mutate_pop(born_animals) #create and mutate offspring (use mutated genes as parent genes), ordered with respect to environment        
        self._animals = new_animals 
        self._positions = self.positions() #update positions  
        if self._constants["verbose"]:
            print("\n\nAnimals per environment: {0}".format(self._positions))
            print("Population size: {0}\tMean payoff: {1:.2f}".format(population_size,mean_payoff))

 #   @jit
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
#             if lifetime_payoff[j]<0:
#                 print(j,lifetime_payoff[j])
#                 print(self._animals[j].lifetime_payoff(self._positions))
        offspring     = np.random.poisson(lam=payoff_factor)        
        pos=np.vectorize(lambda x: x.position)
        positions=pos(self._animals)       
        for j in range(nE):
            d=np.sum(offspring[positions==j])-self._constants["environment_sizes"][j]
            if d>0:
                if self._constants["random_choice"]:    
                    Ind=np.arange(0,len(offspring))
                    I=Ind[(positions==j) & (offspring>0)] # array of relevant indices to choose from
                    if len(I)==0: 
                        I=Ind[positions==j]
                    for i in range(d):  
                        stop=False
                        while not stop:
                            m=np.random.choice(I)
                            if  offspring[m]!=0: #animals with zero offspring are disregarded unless there are no animals in the environment
                                offspring[m]-=1
                                stop=True
                else:
                    pf=np.array(payoff_factor)
                    mx=np.max(pf)
                    pf[(offspring==0)|(positions!=j)]=mx+1 #to make sure that no animals with offspring 0 or from other environment are selected
                    for i in range(d):                         
                        m=np.argmin(pf)
                        offspring[m]-=1
                        if offspring[m]==0:
                            pf[m]=mx+1                        
        born_animals   = np.repeat(self._animals,offspring)
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

