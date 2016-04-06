# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:42:59 2016

@author: matthias
"""
import os
import time

file="main_variable.py"
path= "/Users/matthias/Documents/popdyn/botero-model/single_runs/"
with open(path+"tau_variable.txt","r") as f:
    line=f.readline()
    while line:
        print(line)
        #start=time.time()
        c="python "+path+file+" "+line[0:-1]
        os.system(c)
       # end=time.time()        
        #print("Time needed: {0:4f} min\n".format((end-start)/60))
        line=f.readline()
            

 




