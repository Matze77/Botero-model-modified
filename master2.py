# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:24:33 2016

@author: matze
"""
import numpy as np
import os
file="main_variable.py"
path= "/Users/matthias/Documents/popdyn/botero-model/"

params1=[[[0,0.10,"CBH",1,2.1],[0,0.90,"RP",1,1.2]]]   #cbh-rp
params2=[[[2,0.10,"DBH",1,3.80],[2.0,0.90,"IP",1,1.1]]] #dbh-ip        
params3= [[[2.5,0.70,"IP",2,1.3],[3.5,0.70,"AT",1,1.1]]]  #ip-at
params4=[[[2,0.90,"IP",1,1.11],[0,0.90,"RP",1,2.0]]] #rp-ip
params5 =[[[0,0.10,"CBH",1,2.3],[2,0.10,"DBH",1,7.0]]] #cbh-dbh
params6=[[[2,0.10,"DBH",1,7.3],[3.5,0.10,"AT",1,1.5]]] #dbh-at



for s,pair in enumerate(params4):
    R=[]
    P=[]
    strat=[]
    desc=[]
    path1=[]
    use_pop=[]
    Q=[]
    for i,param in enumerate(pair):
        R.append(10**param[0])
        P.append(param[1])  
        strat.append(param[2])
        desc.append("R{0:.2f}_P{1:.2f}".format(R[i],P[i]))
        path1.append(path+"Output_to_analyze/botero_compare/normal/"+desc[i]+"/")
        use_pop.append(param[3])
        Q.append(param[4])
        
    for i,param in enumerate(pair): 
        if i==1:
            continue
#        q=1.0
#        stop=False
#        while not stop:                    
#            p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --folder base_special --path {3} --q {4}\
#            --plot_every 0 --stop_half 0 --populations 5 --generations 5000 --use_pop {5} --survival_goal 1".format(R[i],P[i],desc[i],path1[i],q,param[3]) 
#            c="python "+path+"single_runs/"+file+p
#            os.system(c)
#            with open(path+"/single_runs/output_variable/base_special/"+desc[i]+"_"+str(q)+"/__overview.txt") as f3:
#                b=f3.read()               
#            s=b[-27:-22]
#            s1=float(s.rsplit("/")[0])
#            s2=float(s.rsplit("/")[1])
#            if s1/s2==1:
#                stop=True
#            else:  
#                q+=0.1
#                q=round(q,1)
#    
        transition=strat[(i+1)%2]+"-"+strat[i]
        if strat[i] in ["RP","IP"] and strat[(i+1)%2] not in ["RP","IP"]:
            modes=[["rand_disc",1,1],["normal",0,0],["rand",1,0]]
        elif strat[(i+1)%2] in ["RP","IP"]:
            modes=[["normal",0,0],["rand_disc",1,1]]
        else: 
            modes=[["normal",0,0]]
        for k in modes:            
             
            for j in [["normal",0.,0.],["hgt_c",1,0],["hgt_r",0,0]]:#[["mut_0.001",0.001,0.0],["sc_0.05",0.0,0.05],["mut_0.001_sc_0.05",0.001,0.05]]
                q=Q[i] 
                stop=False
                while not stop:  
                
                    if j[0][0]!="h":
                        p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --trans 1 --folder {9}/trans_{5} --path {3} --q {4}\
                        --plot_every 100 --populations 20 --generations 5000 --stop_half 1 --survival_goal 0.5 --stop_below 0.5 4000 --use_pop {6} --mutation 0.001 {7} 0.05 {8} --random_a_b {10} --discrete_s {11}".format(R[i],P[i],transition,path1[(i+1)%2],q,j[0],use_pop[(i+1)%2],j[1],j[2],k[0],k[1],k[2]) 
                    else:  
                        p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --trans 1 --folder {9}/trans_{5} --path {3} --q {4}\
                        --plot_every 100 --populations 20 --generations 5000 --stop_half 1 --survival_goal 0.5 --stop_below 0.5 4000  --hgt 1 --check {6} --use_pop {8}  --kt {7} --random_a_b {10} --discrete_s {11} ".format(R[i],P[i],transition,path1[(i+1)%2],q,j[0],j[1],j[2],use_pop[(i+1)%2],k[0],k[1],k[2]) 
                    
            
#                p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --trans 1 --folder trans_{5} --path {3} --q {4}\
#                --plot_every 50 --populations 10 --stop_below 0.5 800 --generations 1000 --random_a_b {6} --discrete_s {7} --survival_goal 0.5  --stop_half 1 --use_pop {8}".format(R[i],P[i],transition,path1[(i+1)%2],q,k[0],k[1],k[2],use_pop[(i+1)%2]) 
                    c="python "+path+"single_runs/"+file+p
                   # print(p)
                    os.system(c)
                    with open(path+"single_runs/output_variable/"+k[0]+"/trans_"+j[0]+"/"+transition+"_"+str(q)+"/__overview.txt") as f3:
                            b=f3.read()               
                    s=b[-27:-22]
                    s1=float(s.rsplit("/")[0])
                    s2=float(s.rsplit("/")[1])
                  #  S.append(s1)
        #            S1=np.array(S[1:4])
        #            S2=np.array(S[4:])
        #            if any(S1>9) and any((S2>9)):
                    if s1/s2>=0.5:
                        stop=True
                    else:
                        q+=0.01
                        q=round(q,2)

     
          



