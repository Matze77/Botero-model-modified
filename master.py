# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:42:59 2016

@author: matthias
"""
import os
from multiprocessing import Pool
import numpy as np

p=4
file="main_constant.py"
path= "/Users/matthias/Documents/popdyn/botero-model/single_runs/"

def run(lines1):
    for line in lines1:
        if len(line)!=0:
          #  print(line)
            c="python "+path+file+" "+str(line)
            print(c)
            os.system(c)

lines=[]

with open(path+"runs/constant_cbh+plast.txt","r") as f:
    #file=f.readline()
    line=f.readline()
    lines.append(line)
    while line and line!="":
        line=f.readline()
        lines.append(line)       
    
'''Create list of lists of p elements to be used as arguments in pool.map '''
list1=[]
list2=[] 
for i,l in enumerate(lines):
    list2.append(l[:-1])
    if (i % p)==p-1:
        list1.append(list2)
        list2=[]
for i in range(p-len(list2)):
    list2.append('')
list1.append(list2)
list1=np.array(list1).T
if len(list1)!=p:
    raise
'''Run p number of processes simultaneously'''


pool=Pool(processes=p)
pool.map(run,list1) 
pool.terminate()