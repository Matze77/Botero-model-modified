# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:24:33 2016

@author: matze
"""
import numpy as np
import os
file="main_constant.py"
path= "/Users/matthias/Documents/popdyn/botero-model/fair_stop/"


R=3.0
dir=1
stop=False
step0=0.32
step=step0
counter=0
while not stop:
    dirs=os.listdir(path+"output/cbh_dbh")
    desc="R{0:.2f}".format(R)
    if desc in dirs:
        desc=desc+"_2"

#    
    p=" --environment {0} 0 1 0 0 --desc {1} --folder cbh_dbh --h_random 1 --plot_every 80 --populations 1 --generations 8000 ".format(round(R,2),desc) 
     
                         
    c="python "+path+file+p
       # print(p)
    os.system(c)
    genes=np.genfromtxt(path+"output/cbh_dbh/"+desc+"/final_genes.csv",delimiter=",",skip_header=1,skip_footer=1) 
    print(genes)
    
    if (genes[6]>0.2 and genes[6]<0.8 and dir==1 and np.abs(genes[0]-genes[1])>0.3)or ((genes[6]<0.2 or genes[6]>0.8 or np.abs(genes[0]-genes[1])<0.3)and dir==-1):
        step/=2
        dir*=-1
        counter=0
        if step<0.05:
            stop=True
    elif counter>=1 and step!=step0:
        step/=2
        if step<0.01:
            stop=True
    R+=step*dir
    counter+=1
        

     
          



