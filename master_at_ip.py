# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:24:33 2016

@author: matze
"""
import numpy as np
import os
file="main_constant.py"
path= "/home/matze/Dokumente/popdyn/Botero-model/"


R=1000
P=0.8
step0=1024

f=open("ip_at.csv","w")
f.write("Size,R_trans\n")
size=5000
tau=0.25
#for size in np.arange(5000,1000,-1000):
for tau in np.arange(0.25,0.05,-0.05):
    dir=1
    stop=False
    step=step0
    counter=0
    while not stop:
        #dirs=os.listdir(path+"output/cbh_dbh")
        desc="R{0:.2f}".format(R)
        #if desc in dirs:
         #   desc=desc+"_2"
    
    #    
        p=" --environment {0} {1} 1 0 0 --desc {2} --folder ip_at --size {3} --plot_every 200 --populations 4 --tau {4} --generations 1000".format(round(R,2),P,desc,size,tau) 
         
                             
        c="python "+path+file+p
           # print(p)
        os.system(c)
        genes=np.genfromtxt(path+"output/ip_at/"+desc+"/final_genes.csv",delimiter=",",skip_header=1,skip_footer=1) 
        #print(genes)
        
        if (genes[8]<=0.5 and dir==1)or (genes[8]>0.5 and dir==-1):
            step/=2
            dir*=-1
            counter=0
            if step<16:
                break
        elif counter>=1 and step!=step0:
            step/=2
            if step<16:
                break
        R+=step*dir
        counter+=1
    if dir==1:
        R+=step
    f.write("{0},{1}\n".format(size,round(R,2)))
    if genes[8]<=0.5:
        R-=2*step
            
f.close()
     
          



