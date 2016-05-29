# -*- coding: utf8 -*-
"""
#########################################################
#
#    animal.pyx
#    Author: Dion Häfner (dionhaefner@web.de)
#    
#    Implements animal class and methods
#    (cythonized version of animal.py and genome.py)
#    
#    YOU WILL NEED TO RUN setup.py AFTER MAKING CHANGES HERE
#
#    Licensed under BSD 2-Clause License
#
#########################################################
"""

import numpy as np

from constants import model_constants

# imports from c libraries for speed
cimport numpy as np
from libc.math cimport abs as c_abs
from libc.math cimport log as c_log
from libc.math cimport exp as c_exp
from libc.stdlib cimport rand as c_rand
from libc.stdlib cimport RAND_MAX
from libc.math cimport fmax as c_max
from libc.math cimport fmin as c_min
from cpython cimport bool #for boolean type

# custom type definitions
DTYPE = np.int
BTYPE = np.uint8
ctypedef np.int_t DTYPE_t
ctypedef np.uint8_t BTYPE_t

constants = model_constants
cdef bool hgt=constants["hgt"]
cdef double tau=constants["tau"]
cdef double kd=constants["kd"]
cdef double ka=constants["ka"]
cdef double kt=constants["kt"]
cdef double mut1=float(constants["mutation"][0]) 
cdef double scale_mu=float(constants["mutation"][1]) #standard deviation for mutation of mutation rate
cdef double sc1=float(constants["mutation"][2]) 
cdef double scale_sc=float(constants["mutation"][3]) 
cdef bool random_a_b=constants["random_a_b"]


cdef class Animal:
    """Implements a cython class Animal, that is also available outside of this module """
    # properties of Animal, typed as C variables for speed
    cdef object _constants
    cdef bool primed
    # public keyword makes the variable accessible to python
    cdef public double h,s,a,I0,I0p,b,bp,mu,sc,t
    cdef public int adjustments
    cdef public int transfers
    cdef public double insulation
    cdef public double mismatch
    cdef public int lineage

    def __init__(self,np.ndarray[double,ndim=1] parent_genes=np.array([]),int lineage=1):
        """Constructor"""
        self._constants = model_constants
        if parent_genes.size: # empty argument -> random genes (default)
            self.genes = parent_genes

        else:
            self.genes = random_genes()

            

        
        self.mismatch = 0
        self.adjustments = 0 
        self.transfers=0
        self.insulation = self.genes[3]
        self.lineage=lineage
        self.primed = False #defines whiche gene set to use (normal or alternative)
        

# PUBLIC METHODS
    cpdef choose_set(self,double r):        
        """At birth of animal:Determine if primed or unprimed gene set (I0,b) is used"""
        if (r <= self.h): 
            self.primed=False
        else:
            self.primed=True
            self.insulation=self.I0p

    cpdef react(self,double E, double C, double r,BTYPE_t evolve_all=0):
        """Animal migrates and reacts to environment E and cue C. If evolve_all is set, reaction takes place for all animals, regardless of gene 'a'."""
        cdef float new_insulation        
        if self.s > 0.5:         
            if ((r<= self.a) | evolve_all):#draws random number to determine whether the animal adjusts its insulation (just for plastic ones)
                if self.primed:
                    new_insulation = self.I0p+self.bp*C
                else:
                    new_insulation = self.I0+self.b*C
                self.insulation = new_insulation 
                self.adjustments+=1  #raise adjustments
        self.mismatch += c_abs(self.insulation-E)

    cpdef lifetime_payoff(self):
        """Assembles the lifetime payoff of the animal"""
        cdef double costs
        costs=0
        if (self.s > 0.5):
            costs+=kd + (self.adjustments-1) * ka+self.transfers * kt  #ka is payed for every phenotype adjustment (except the obligatory first one), kd is payed every round 
        p= c_max(c_exp(-tau*self.mismatch)-costs, 0) 
        return p

    cpdef mutate(self):
        """Causes the Animal's genes to mutate"""
        cdef np.ndarray[double,ndim=1] new_genes = self.genes
        cdef int k
        cdef double r, mutation_step  
        cdef np.ndarray[double,ndim=1] r1      
        if scale_mu>0.0:
            r = np.random.rand()  
            if (r<=self.mu):
                mutation_step = np.random.normal(loc=0,scale=scale_mu)  #mutate mutation rate with much lower step size
                new_genes[7] += mutation_step
                self.mu=new_genes[7]
                
        if scale_sc>0.0:
            r = np.random.rand()  
            if (r<=self.mu):
                mutation_step = np.random.normal(loc=0,scale=scale_sc)  #mutate scale of mutation 
                new_genes[8] += mutation_step
                self.sc=new_genes[8]
                
        r1=np.random.rand(3)
        for i,k in enumerate([0,3,4]): #genes modified for all individuals: h, I0, I0'
            if (r1[i]<=self.mu):
                mutation_step = np.random.normal(loc=0,scale=self.sc)
                new_genes[k] += mutation_step
                
        r=np.random.rand()  
        if (r<=self.mu):
                if constants["discrete_s"]:
                    new_genes[1]=(new_genes[1]+1)%2
                else:
                    mutation_step = np.random.normal(loc=0,scale=self.sc)
                    new_genes[1] += mutation_step
                if random_a_b:
                    if self.s<0.5 and new_genes[1]>0.5: #transition to plasticity happened: choose a,b,b' at random
                        r1=np.random.rand(3)
                        new_genes[2]=r1[0]
                        new_genes[5]=r1[1]*4-2
                        new_genes[6]=r1[2]*4-2
                
        
        if new_genes[1] > 0.5: # other genes modified if individual is plastic (s>0.5): a , b, b'
            r1=np.random.rand(3)    
            for i,k in enumerate([2,5,6]):
                if (r1[i]<=self.mu):
                    mutation_step = np.random.normal(loc=0,scale=self.sc)
                    new_genes[k] += mutation_step                             
        else:
            new_genes[2], new_genes[5], new_genes[6] = 0, 0, 0 # for non-plastic individuals, set a,b,b' to 0
            
        if hgt:  
            r=np.random.rand()  
            if (r<=self.mu):
                mutation_step = np.random.normal(loc=0,scale=self.sc)
                new_genes[9] += mutation_step

        return new_genes



    property gene_dict:
        """Allows the genes of the animal to be read from python as a dict by calling animal.gene_dict"""
        def __get__(self):
            return {"h":self.h,"s":self.s,"a":self.a,"I0":self.I0,"I0p":self.I0p,"b":self.b,"bp":self.bp,"mu":self.mu,"sc":self.sc,"t": self.t}

    property genes:
        """Allows the genes of the animal to be read and written from python as a list by calling animal.genes"""
        def __get__(self):
            return np.array([self.h,self.s,self.a,self.I0,self.I0p,self.b,self.bp,self.mu,self.sc,self.t])

        def __set__(self, object genes):
            self.set_genes(genes)

    cdef set_genes(self,np.ndarray[double,ndim=1] genes):
        """Sets genes and prevents mutation of them out of the respective sensible interval """

        self.h = c_max(0,c_min(1,genes[0]))
        if not constants["force_plast"]:
            self.s = c_max(0,c_min(1,genes[1]))
        else:
            self.s = c_max(0.501,c_min(1,genes[1])) #forces plasticity s>0.5
        self.a = c_max(0,c_min(1,genes[2]))    
        self.I0 = genes[3]
        self.I0p = genes[4]
        self.b = genes[5]
        self.bp = genes[6]
        self.mu = c_max(mut1,c_min(1,genes[7]))  #to prevent mutation from falling below initial value
        self.sc=c_max(0,genes[8])   
        self.t=c_max(0,c_min(1,genes[9]))       



# PROTECTED FUNCTIONS



cdef inline np.ndarray[double,ndim=1] random_genes():
    """Returns random values for the 10 genes in the chosen intervals:
    h: 1, s: [0,1], a: [0,1], I0: [-1,1], I0p: [-1,1], b: [-2,2], bp: [-2,2] ,mu: mut1, sc: sc1, t: [0,1], ta: [0,1]"""
    cdef np.ndarray[double,ndim=1] rand_numbers, rand_genes 
    cdef str distr_mut
    cdef double mut2,r,s1,s2,t1
    rand_numbers = np.random.rand(8)
    s1=1 #values for plasticity trait
    s2=0
    if constants["discrete_s"]:
        rand_numbers[1]=np.random.randint(2)
    if constants["force_plast"]: #to force s in [0.5,1]
        s1=0
        s2=1
    t1=0
    if constants["hgt"]:
        t1=1
    rand_genes = [0,s1,1,2,2,4,4,t1]*rand_numbers+[1,s2,0,-1,-1,-2,-2,0]

    if (rand_genes[1]<=0.5):
        rand_genes[2], rand_genes[5], rand_genes[6]  = 0, 0, 0
     
    rand_genes=np.insert(rand_genes,7,mut1)  #insert mutation rate in rand_genes before t and ta
    rand_genes=np.insert(rand_genes,8,sc1)
        
    return rand_genes


    

