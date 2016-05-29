# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:05:02 2016

@author: matthias
"""


import csv
import numpy as np
import scipy as sp
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import os
import shutil
from multiprocessing import Pool

path="/Users/matthias/Documents/popdyn/botero-model/single_runs/"
folder="trans_cbh_plast"
r1=np.arange(0.41,1.31,0.01)
r2=np.arange(1.35,3.55,0.05)
r=np.array(list(r1)+list(r2))
pops=7
p=4

inc=0.01
f1=open("trans_to_plast_new.csv","w")
f1.write("R,Ptrans\n")

def run(r1):
    tr=[]
    for i,R in enumerate(r1):

        if R==0:
            continue
        dir=1 #direction of change of P: 1 increasing, -1 decreasing
        P=0.24
        if R>=2.5:
            P=0.21
        stop=False
        counter=0
        while not stop:
            
            k=0 #counter for plasticity for the 5 populations
            desc="R{0:.2f}_P{1:.2f}".format(R,P)
                
            p=" --environment {0} {1} 1 0 0 --desc {2} --folder {3} \
                --plot_every 0 --populations {4} --generations 800".format(R,P,desc,folder,pops) 
            c="python "+path+"main_constant.py"+p
            print(c)
            os.system(c)
            for pop in range(1,pops+1):
                pathf=path+"output/{3}/R{0:.2f}_P{1:.2f}/pop{2}_final_state.csv".format(R,P,pop,folder)
                data = np.genfromtxt(pathf,skip_header=1,delimiter=",") 
                s=data[:,1]
                if len(s[s>0.5])>len(s)/2:
                    k+=1
           
            if k>pops/2: 
                if counter==0:
                    dir=-1  #transition below: decrease 
                elif dir==1:
                    tr.append([round(R,2),round(P,2)])
                    break
            elif dir==-1:
                tr.append([round(R,2),round(P+inc,2)])
                break
            P+=inc*dir
            counter+=1
            
    return tr                     

f1=open("trans_to_plast_new.csv","w")
f1.write("R,Ptrans\n")
           
list1=[]
list2=[] 
for i,l in enumerate(r):
    list2.append(l)
    if (i % p)==p-1:
        list1.append(list2)
        list2=[]
for i in range(p-len(list2)):
    list2.append(0)
list1.append(list2)
list1=np.array(list1).T
if len(list1)!=p:
    raise
'''Run p number of processes simultaneously'''


pool=Pool(processes=p)
trans=pool.map(run,list1) 
pool.terminate()
trans=[a for sublist in trans for a in sublist]         
for i in trans:
    #print(str(i)[1:-1])
    f1.write(str(i)[1:-1]+"\n")
f1.close()

plt.figure()
plt.plot(r,trans,"-*",markersize=3,lw=0.5)
plt.ylabel("$P_{\mathrm{trans}}$")
plt.ylim(0.2,0.38)
plt.xlabel("R")  
plt.grid()
plt.savefig("trans_to_plast1.pdf")   

#plt.figure()
#plt.plot(r[91:],trans[91:],"-*",markersize=3,lw=0.5)
#plt.ylabel("$P_{\mathrm{trans}}$")
#plt.ylim(0.17,0.26)
#plt.xlabel("R")  
#plt.grid()
#plt.savefig("trans_to_plast2.pdf")   
#
#ra=[]


            
