# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:42:59 2016

@author: matthias
"""
import os
from multiprocessing import Pool

p=1
file="main_constant.py"
path= "/Users/matthias/Documents/popdyn/botero-model/single_runs/"

def run(line):
        c="python "+path+file+" "+line
        print(c)
        os.system(c)

lines=[]

with open(path+"/runs/constant_run2.txt","r") as f:
    line=f.readline()
    lines.append(line)
    while line and line!="":
        line=f.readline()
        lines.append(line)       
    
'''Create list of lists of p elements to be used as arguments in pool.map '''
list1=[]
list2=[] 
for i,l in enumerate(lines):
    list2.append(l[0:-1])
    if (i % p)==p-1:
        list1.append(list2)
        list2=[]
list1.append(list2)

     
  

    
'''Run p number of processes simultaneously'''

a=[]
for l in list1:
    if len(l)!=0:
        pool=Pool(processes=len(l))
        a.extend(pool.map(run,l))    
        pool.terminate()

