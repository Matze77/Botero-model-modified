# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:24:33 2016

@author: matze
"""
import numpy as np
import os
file="main_variable.py"
path= "/Users/matthias/Documents/popdyn/botero-model/"

#params=[[[0,0.20,"CBH",1],[0,0.30,"RP",1]],    #cbh-rp
#params =[[[0,0.10,"CBH",1],[1,0.10,"DBH",1]]] #cbh-dbh
#params=[[[1,0.10,"DBH",1],[1.0,0.20,"IP",3]]]  #dbh-ip        
#params=[[[2.5,0.20,"DBH",1],[3,0.20,"AT",1]]] #dbh-at
#params=[[[0.5,0.50,"RP",2],[1.5,0.50,"IP",1]]] #rp-ip
#params= [[[3,0.90,"IP",2],[3.5,0.90,"AT",1]]]  #ip-at

#        
params_flat=[item for sublist in params for item in sublist]

#f1=open("./runs/base_special.csv",'w') 
#f2=open("./runs/trans.csv",'w')
#f1.write("R,P")
#f2.write("Transition,q,normal,mut(std=0.001),mut(std=0.01),mut(std=0.1),hgt_random,hgt_check\n")
#for Q in np.arange(2,4.1,0.1):
  #  f1.write(","+str(Q))
   # 
#f1.write("\n")
for j,pair in enumerate(params):
    R=[]
    P=[]
    strat=[]
    desc=[]
    path1=[]
    use_pop=[]
    for i,param in enumerate(pair):
        R.append(10**param[0])
        P.append(param[1])  
        strat.append(param[2])
        desc.append("R{0:.2f}_P{1:.2f}".format(R[i],P[i]))
        path1.append(path+"/output_to_analyze/botero_compare/new/"+desc[i]+"/")
        use_pop.append(param[3])
    for i,param in enumerate(pair): 
        #(param)
       # print(q)
        q=1.0
    #    f1.write("{0},{1}".format(round(R[i],2),P[i]))
        stop=False
        while not stop:                    
            p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --folder base_special --path {3} --q {4}\
            --plot_every 0 --stop_half 0 --populations 5 --generations 5000 --use_pop {5} --survival_goal 1".format(R[i],P[i],desc[i],path1[i],q,param[3]) 
            c="python "+path+"single_runs/"+file+p
           # print(p)
            os.system(c)
            with open(path+"/single_runs/output_variable/base_special/"+desc[i]+"_"+str(q)+"/__overview.txt") as f3:
                b=f3.read()               
            s=b[-27:-22]
            s1=float(s.rsplit("/")[0])
            s2=float(s.rsplit("/")[1])
         #   f1.write(","+s)
            if s1/s2==1:
                stop=True
            else:  
                q+=0.1
                q=round(q,1)
    
  #      f1.write("\n")
        transition=strat[(i+1)%2]+"-"+strat[i]
 #       f2.write(transition)
    
        stop=False
        while not stop:   
            S=[]
  #          f2.write(","+str(q))
            for j in [["normal",0.0],["mut_0.001",0.001],["mut_0.01",0.01],["mut_0.1",0.1],["hgt_r",0,0,0],["hgt_c",1,0,0],["hgt_c2",1,0.01,0.01]]:
                if j[0][0]=="m" or j[0][0]=="n":
                    p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --trans 1 --folder trans_{5} --path {3} --q {4}\
                    --plot_every 50 --populations 20 --generations 8000 --use_pop {7} --mutation normal 0.001 0.0 {6}".format(R[i],P[i],transition,path1[(i+1)%2],q,j[0],j[1],use_pop[(i+1)%2]) 
                else:  
                    p=" --environment {0} {1} 1 0 0 --desc {2}_{4} --trans 1 --folder trans_{5} --path {3} --q {4}\
                    --plot_every 50 --populations 20 --generations 8000 --hgt 1 --check {6} --use_pop {9} --kh {7} --kt {8}".format(R[i],P[i],transition,path1[(i+1)%2],q,j[0],j[1],j[2],j[3],use_pop[(i+1)%2]) 
                c="python "+path+"single_runs/"+file+p
               # print(p)
                os.system(c)
                with open(path+"single_runs/output_variable/trans_"+j[0]+"/"+transition+"_"+str(q)+"/__overview.txt") as f3:
                        b=f3.read()               
                s=b[-27:-22]
                s1=float(s.rsplit("/")[0])
                s2=float(s.rsplit("/")[1])
                S.append(s1)
        #       f2.write(","+s)
         #   f2.write("\n")
            print(S)
            S1=np.array(S[1:4])
            S2=np.array(S[4:])
            if any(S1>9) and any((S2>9)):
                stop=True
            else:
                q+=0.1
                q=round(q,1)

     
#f1.close()
#f2.close()
            
          



