#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#    main_variable.py
#    Author: Dion Häfner (dionhaefner@web.de)
#    
#    Main controller, with variable population size
#
#    Usage:
#        python main_variable.py [options]
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
from multiprocessing import Process,Pool
import numpy as np # For efficient array operations
#import matplotlib.pyplot as plt # For plotting
import time # For timing parts of the script, optimizing run time
#import pandas as pd # Easier data handling
import os # To create directories
import datetime # To access the current time
import sys # To access command line arguments
#import warnings # To warn the user
import csv

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



if __name__ == '__main__':
    # Get model constants
    constants = model_constants
    population_size=constants["size"]
    
    if constants["path"] !="":
        path=constants["path"]

    else:  
        p=open("./path.txt","r")
        path=p.readline()

    final_state=path+ "pop{0}_final_state.csv".format(constants["use_pop"])  #specify path to csv file in constants (from run with constant pop size) 
    f_mean = path + "pop{0}_mean_genes.csv".format(constants["use_pop"]) #specify path to csv file in constants (from run with constant pop size) 

    # create output directory
    now = datetime.datetime.today()
    start = time.clock()

    # read the csv files
  
    with open(f_mean) as f:
        reader = csv.reader(f,delimiter=",")
        for (i,row) in enumerate(reader):
            if row and row[0]!="":
                if row[0]=="n":
                    names=row
                    data = np.genfromtxt(f_mean,skip_header=i+1,delimiter=",") #reads mean genes and n , nperPos from csv 
                    break
                elif (row[0][0]!="R"): 
                    e=row[:5]
                    environment=list(map(float,e))  #reads environment values 
          

    if constants["trans"]:
        factor=constants["environment"][0]/environment[0] #to ensure that environment value E is continous (no jump from constant run)
    else:
        factor=1
    gen=data[-1,0]
    final_t = gen*constants["L"]*factor #final time of constant run
    data = np.genfromtxt(final_state,skip_header=1,delimiter=",") #reads genes and n , nperPos from csv file of the whole final population 
    genes1=data[:,:-1] #last column mismatch is removed
    if constants["discrete_s"]:
        genes1[:,1]=np.round(genes1[:,1]) #convert continous values for s to discrete values
    if names[-2]=="ta":
        genes1=np.delete(genes1,-2,axis=1) #delete t
        genes1=np.insert(genes1,8,constants["mutation"][2],axis=1) #insert sc

    mean_genes=np.mean(genes1,axis=0)
    std_genes=np.std(genes1,axis=0)             

    # create environment and output information about them

    if constants["trans"]: #for transition runs, change environment parameters in constants file
        env = Environment(*constants["environment"]) 
    else: #or use environment from constant run
        env = Environment(*environment)
        
    if constants["folder"]!="":       
        path = "./output_variable/{0}/".format(constants["folder"])
    else:
        path = "./output_variable/"
        
    if constants["time_tag"] or constants["desc"]=="":        
        path=path+"{0:%y}-{0:%m}-{0:%d}_{0:%H}-{0:%M}-{0:%S}-{1}/".format(now,constants["desc"])  
    else:
        path=path+constants["desc"]+"/"
    
    try: 
        os.makedirs(path+"timeseries/")
        if constants["save_all"]:
                os.makedirs(path+"all_genes/")
    except OSError:
        if not os.path.isdir(path):
            raise
    f3 = open(path+"__overview.txt",'w')
    f3.write("initial conditions \n")
    f3.write("R,P,A,B,O\n{0},{1},{2},{3},{5}\n".format(env.R,env.P,env.A,env.B,i,env.O))
    f3.write("Mean genes(h,s,a,I0,I0p,b,bp,mu,sc,t):\n{0}\n".format(mean_genes))
    f3.write("Std genes:\n{0}\n".format(std_genes))
    for key in ['generations','L','kd','ka','tau','q','mutation','random_a_b',"discrete_s",'environment','environment_name','size','populations','plot_every','verbose',\
'random_choice','desc','force_plast',"hgt",'check','kt',"proc",'trans','path','use_pop','stop_half',"stop_below","survival_goal"]:
        f3.write("{0}:\t{1}\n".format(key,constants[key]))
    
    end = time.clock()
    if constants["verbose"]:
        print("Set-up time: {0:.2e}s\n".format(end-start))


    # main loop over multiple populations
    survival_rate = 0
    rand1=np.random.randint(0,10**3) #different random number seed for every time the program is started
    
    def main(k):
        start = time.clock()
        np.random.seed(rand1+k) #to prevent the same seeds for simultaneously run populations 
            # write starting genes in files

        f1 = open(path+"pop"+str(k+1)+"_mean_genes.csv",'w')
        f1.write("\nn,I0,I0p,payoff,a,b,bp,h,mu,s,sc,t,size,lin\n")        
        f2 = open(path+"pop"+str(k+1)+"_std_genes.csv",'w')
        f2.write("\nn,I0,I0p,payoff,a,b,bp,h,mu,s,sc,t,size\n")
            
        # create animals with the mean genes that shall be tested for each environment
        animals=[]
        animals.append([Animal(np.array(g),lineage=k) for k,g in enumerate(genes1[:constants["size"]])])
        animals = [item for sublist in animals for item in sublist] # flatten animal list 
        
        # create a population of population_size animals that have the correct mean genes
        population = Population(population_size,animals)
        
        pop_mean, pop_std, final_gen = iterate_population(k,population,env,f1,f2,path,final_t,True) #start at final time of constant pop run
        end = time.clock()


        if population.size()==0:
            print("\t Population {0} died out! Time needed: {1:.2f} min".format(k+1,(end-start)/60)) 
            return False,final_gen
        elif float(population.size())/constants["size"] <= constants["stop_below"]:
            print("\t Population {0} fell below critical limit! Time needed: {1:.2f} min".format(k+1,(end-start)/60)) 
            return False,final_gen
        else:
            print("\t Population {0} survived! Time needed: {1:.2f} min".format(k+1,(end-start)/60))
            return True,final_gen

    
    p=constants["proc"]
    survival_rate=0
    if p>1:
        '''Create list of lists of p elements to be used as arguments (population number) in pool.map '''
        list1=[]
        list2=[] 
        for k in range(constants["populations"]):
            list2.append(k)
            if (k % p)==p-1:
                list1.append(list2)
                list2=[]
        list1.append(list2)
        
        '''Run proc number of populations simultaneously'''    
        a=[]
        counter=0
        for l in list1:
            if len(l)!=0:
                pool=Pool(processes=len(l))
                out=pool.map(main,l)
                a.extend(out)
                for x,y in enumerate(out):
                    if y[0]:
                        survival_rate+=1
                        f3.write("Population {0} survived!\n".format(x+counter+1))
                    else:
                        f3.write("Population {0} died at generation {1}!\n".format(x+counter+1,y[1]))  
                if survival_rate>=constants["populations"]/2 and constants["stop_half"]:
                    break
                pool.terminate()
                counter+=p
                
    else:
        a=[]
        for k in range(constants["populations"]):
            out=main(k)
            a.append(out)
            if out[0]:
                survival_rate+=1
                f3.write("Pop {0} survived!\n".format(k+1))
            elif constants["stop_below"]==0:
                f3.write("Pop {0} died at generation {1}!\n".format(k+1,out[1]))  
            else:
                f3.write("Pop{0} fell below critical limit at generation {1}!\n".format(k+1,out[1]))  
            if float((survival_rate+constants["populations"]-k-1))/constants["populations"]<constants["survival_goal"]:
                break
            if survival_rate>=constants["populations"]/2 and constants["stop_half"] and k<constants["populations"]/2:
                break
                  
    f3.write("\n\nIn total,   {0}/{1} Populations survived.".format(survival_rate,len(a)))
    f3.close()





















































