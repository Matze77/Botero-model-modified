# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 15:47:06 2016

@author: matthias
"""

animal_list=[]
animal_list.extend([Animal(np.array([]),position=0,lineage=_) for _ in range(1000)])
population = Population(1000,animal_list)

@jit
def syst_pop_control(self,d,pf,offspring,positions,j):
    if d>0:#if environment overcrowded let the least fit animals have less offspring
        mx=np.max(pf)
        pf[(offspring==0)|(positions!=j)]=mx+1 #to make sure that no animals with offspring 0 or from other environment are selected
        for i in range(d):                         
            m=np.argmin(pf)
            offspring[m]-=1
            if offspring[m]==0:
                pf[m]=mx+1
            
    elif d<0:#if too few, clone the fittest animals
        pf[positions!=j]=0
        for i in range(-d):                     
            m=np.argmax(pf)
            pf[m]=0 #to ensure that not all additional offspring is from one animal
            offspring[m]+=1                  
    return offspring
    
def breed_constant(self):
    """Iterates the entire Population to a new generation, calculating the number of offspring of each Animal with CONSTANT population size"""
    nE=len(self._constants["environments"])
    pos=np.vectorize(lambda x: x.position)
    positions=pos(self._animals)
    calc_payoff     = np.vectorize(lambda x: x.lifetime_payoff(self._positions))
    lifetime_payoff = calc_payoff(self._animals)
    mean_payoff     = [np.mean(lifetime_payoff[positions==j]) for j in range(nE)]
    population_size=sum(self._constants["environment_sizes"])
    if (mean_payoff == 0):
        raise RuntimeError("Mean payoff of population decreased to 0. Check your parameters!")
    else:
        try:
            payoff_factor = np.array([lifetime_payoff[positions==j]/mean_payoff[j] for j in range(nE)]).flatten()#flatten array
          
            if payoff_factor.shape[0]!=population_size:                   
                raise                   
        except:
            payoff_factor = np.hstack(np.array([lifetime_payoff[positions==j]/mean_payoff[j] for j in range(nE)]))
            # if array has different dimensions (due to migration) in each row
    offspring = np.random.poisson(lam=payoff_factor) #number of offspring drawn from Poisson distr for each animal 
    for j in range(nE):
        d=np.sum(offspring[positions==j])-self._constants["environment_sizes"][j]
        if self._constants["random_choice"]:    
            Ind=np.arange(0,len(offspring))
            I=Ind[(positions==j) & (offspring>0)] # array of relevant indices to choose from
            if len(I)==0: 
                I=Ind[positions==j]
            if d>0:#if environment overcrowded reduce offspring randomly
                add=-1
            else:#if too few, increase offspring randomly
                add=1
            for i in range(abs(d)):  
                stop=False
                while not stop:
                    m=np.random.choice(I)
                    if  offspring[m]!=0 or sum(offspring[positions==j])==0: #animals with zero offspring are disregarded unless there are no animals in the environment
                        offspring[m]+=add
                        stop=True
        else:
            pf=np.array(payoff_factor)
            offspring=self.syst_pop_control(d,pf,offspring,positions,j) 
#            if d>0:#if environment overcrowded let the least fit animals have less offspring
#                mx=np.max(pf)
#                pf[(offspring==0)|(positions!=j)]=mx+1 #to make sure that no animals with offspring 0 or from other environment are selected
#                for i in range(d):                         
#                    m=np.argmin(pf)
#                    offspring[m]-=1
#                    if offspring[m]==0:
#                        pf[m]=mx+1
#                    
#            elif d<0:#if too few, clone the fittest animals
#                pf[positions!=j]=0
#                for i in range(-d):                     
#                    m=np.argmax(pf)
#                    pf[m]=0 #to ensure that not all additional offspring is from one animal
#                    offspring[m]+=1                                               
    born_animals = np.repeat(self._animals,offspring) # Create list with offspring repeats for each animal (cloned animals)      
    mutate_pop = np.vectorize(lambda x: Animal(x.mutate(),x.position,x.lineage)) #animals are created with mutated genes and at the position of their parents
    new_animals = mutate_pop(born_animals) #create and mutate offspring (use mutated genes as parent genes), ordered with respect to environment        
    self._animals = new_animals 
    self._positions = self.positions() #update positions  
%load_ext line_profiler 
%lprun -f breed_constant breed_constant(population)
