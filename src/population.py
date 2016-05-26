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

from numba import jit


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
        if not evolve_all:
            r2=np.random.rand(self._size)        
        else:
            r1=np.random.rand(self._size)
            for i,animal in enumerate(self._animals): 
                animal.choose_set(r1[i]) #After birth: choose primed or unprimed gene set for each animal
            r2=np.empty(self._size) #random numbers for adjustments not needed, as adjustment after birth is obligatory 
        if self._constants["hgt"]:  
            r3=np.random.rand(self._size)
            for i,animal in enumerate(self._animals):    
                if animal.t > 0.5:
                    if (r3[i] <= animal.ta):
                        j=np.random.randint(0,self._size) #animal from which to choose from
                        if self._constants["check"]: #Payoff of donor animal is checked (directed HGT)
                            if self._animals[j].lifetime_payoff()< animal.lifetime_payoff():
                                animal.react(E,C,r2[i],evolve_all)#animal reacts to environment 
                                continue    #stop this iteration if current payoff of donor animal is smaller than for the recipient animal
                                
                        if float(self._constants["mutation"][3])>0: #if mu is mutable, also this trait can be chosen for hgt
                            g=np.random.randint(0,7) #gene to choose 
                        else: 
                            g=np.random.randint(0,6)
                        gene=self._animals[j].genes[g] #value of this gene (from chosen animal)                
                        genes=animal.genes                    
                        genes[g]=gene                    
                        animal.genes=genes
                        animal.transfers+=1
                        
                animal.react(E,C,r2[i],evolve_all)#animal reacts to environment               
        else:
            for i,animal in enumerate(self._animals):  
                animal.react(E,C,r2[i],evolve_all)#animal reacts to environment
                
    def breed_constant(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with CONSTANT population size"""
        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff())
        lifetime_payoff = calc_payoff(self._animals)
        mean_payoff  = np.mean(lifetime_payoff)
        population_size=self._constants["size"]
        if (mean_payoff == 0):
            raise RuntimeError("Mean payoff of population decreased to 0. Check your parameters!")
        else:
            payoff_factor = lifetime_payoff/mean_payoff                                        
        offspring = np.random.poisson(lam=payoff_factor) #number of offspring drawn from Poisson distr for each animal 
        d=np.sum(offspring)-population_size
        
        if self._constants["random_choice"]:    
            Ind=np.arange(0,len(offspring)) #all indices
            if d>0:#if environment overcrowded reduce offspring randomly
                add=-1
                I=Ind[(offspring>0)] #animals with zero offspring are disregarded
            else:#if too few, increase offspring randomly
                add=1
                I=Ind
            for i in range(abs(d)):  
                m=np.random.choice(I)
                offspring[m]+=add
                if d>0:
                    I=Ind[(offspring>0)] #update index list to prevent picking an individual with 0 offspring

        else:
            pf=np.array(payoff_factor)
            if d>0:#if environment overcrowded let the least fit animal(s) have less offspring
                mx=np.max(pf)
                pf[offspring==0]=mx+1 #to make sure that no animals with offspring 0  are selected
                for i in range(d):                         
                    m=np.argmin(pf)
                    offspring[m]-=1
                    if offspring[m]==0:
                        pf[m]=mx+1
                
            elif d<0:#if too few, let the fittest animals have more offspring          
                for i in range(-d):                     
                    m=np.argmax(pf)
                    pf[m]=0 #to ensure that not all additional offspring is from one animal
                    offspring[m]+=1  
        born_animals = np.repeat(self._animals,offspring) # Create list with offspring repeats for each animal (cloned animals)      
        mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.lineage)) #animals are created with mutated genes of their parents
        new_animals = mutate_pop(born_animals) #create and mutate offspring (use mutated genes as parent genes)      
        self._animals = new_animals 
        if self._constants["verbose"]:
            print("Population size: {0}\tMean payoff: {1:.2f}".format(population_size,mean_payoff))
        return d
 #   @jit
    def breed_variable(self):
        """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with VARIABLE population size"""
        calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff())
        lifetime_payoff = calc_payoff(self._animals)
        payoff_factor=np.array([])
        for j,animal in enumerate(self._animals):
             payoff_factor=np.append(payoff_factor,self._constants["q"]*lifetime_payoff[j])
        offspring     = np.random.poisson(lam=payoff_factor)             
        d=np.sum(offspring)-self._constants["size"]    
        if d>0:
            if self._constants["random_choice"]:    
                Ind=np.arange(0,len(offspring)) #all indices
                I=Ind[(offspring>0)] #animals with zero offspring are disregarded
                for i in range(d):  
                    m=np.random.choice(I)
                    offspring[m]-=1
                    I=Ind[(offspring>0)] #update index list to prevent picking an individual with 0 offspring

            else:
                pf=np.array(payoff_factor)
                if d>0:#if environment overcrowded let the least fit animals have less offspring
                    mx=np.max(pf)
                    pf[offspring==0]=mx+1 #to make sure that no animals with offspring 0 
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
        mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.lineage))
        new_animals = mutate_pop(born_animals)
        self._animals = new_animals        
        N = len(new_animals)
        self._size=N
        if self._constants["verbose"]:
            print("Population size: {0}".format(N))
        return d          
    def lineage(self):
        '''Returns the lineage of each animal'''
        fun = np.vectorize(lambda x: x.lineage)
        lin = fun(self._animals)
        return np.array(lin,ndmin=1) 

