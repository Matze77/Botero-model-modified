# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:24:33 2016

@author: matze
"""
import numpy as np
import os
from multiprocessing import Pool
import shutil
p=2
file="main_variable.py"
path= "/Users/matthias/Documents/popdyn/botero-model"
#f1=open("./runs/base_lowest_q.csv",'w') 
#f1.write("R,P,q\n")
#f2.write("Transition,q,normal,mut(std=0.001),mut(std=0.01),mut(std=0.1),hgt_random,hgt_check\n")
   # 

def run(params):
    for param in params:
        if param==[]:
            continue
        R=param[0]
        P=param[1]
        
        desc="R{0:.2f}_P{1:.2f}".format(R,P)
        path1=path+"/Output_to_analyze/botero_compare/new/"+desc+"/"
        h=np.genfromtxt(path1+"final_genes.csv",delimiter=",",skip_header=1)[0,6]
       # print(h)
        if h<0.8 and h>0.2:
            ngen=3000
        else:
            ngen=500
        q=1.0
        stop=False
        while not stop:                    
            p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --folder base_extinct --path {3} --q {4}\
            --plot_every 0 --stop_half 0 --populations 1 --generations {5} --stop_below 0.5  --mutation normal 0.001 0.0 0.0 --survival_goal 1".format(R,P,desc,path1,q,ngen) 
            c="python "+path+"/single_runs/"+file+p
            print(c)
            os.system(c)
            path2=path+"/single_runs/output_variable/base_extinct/"+desc+"_"+str(q)
            with open(path2+"/__overview.txt") as f3:
                b=f3.read()               
            s=b[-27:-22]
            s1=float(s.rsplit("/")[0])
            s2=float(s.rsplit("/")[1])
            if s1/s2==1:
                stop=True
                
            else:  
                shutil.move(path2,path+"/single_runs/output_variable/base_extinct/extinct/"+desc+"_"+str(q))
                q+=0.1
                q=round(q,1) 

params=[]
for r in np.arange(0,4.5,0.5):
    R=10**r
    for P in np.arange(0,1.1,0.1):
        params.append([R,P])


'''Create list of lists of p elements to be used as arguments in pool.map '''
list1=[]
list2=[] 
for i,l in enumerate(params):
    list2.append(l)
    if (i % p)==p-1:
        list1.append(list2)
        list2=[]
for i in range(p-len(list2)):
    list2.append([])
list1.append(list2)
list1=np.array(list1).T
if len(list1)!=p:
    raise
'''Run p number of processes simultaneously'''


pool=Pool(processes=p)
pool.map(run,list1) 
pool.terminate()


#f1.close()



