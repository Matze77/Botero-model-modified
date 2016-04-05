# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:42:59 2016

@author: matthias
"""
import os
import time

path= "/Users/matthias/Documents/popdyn/botero-model/fair_stop/"
lines=8
params=[]
with open(path+"setting.txt","r") as f:
    line=f.readline()
    while line!="\n":
        print(line)
        #start=time.time()
        c="python "+path+"main_constant.py "+line[0:-1]
        os.system(c)
       # end=time.time()        
        #print("Time needed: {0:4f} min\n".format((end-start)/60))
        line=f.readline()
            


 




