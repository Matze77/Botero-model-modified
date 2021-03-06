    #!/usr/bin/env python
# -*- coding: utf8 -*-
"""
#########################################################
#
#    constants.py
#    Author: Dion Häfner (dionhaefner@web.de)
#    
#    Defining the constants of the Botero evolution model
#    and parsing them from command line.
#    
#    Usage:
#        from constants import model_constants
#    model_constants is then a dict containing all the 
#    parameters.
#
#    Licensed under BSD 2-Clause License
#
#########################################################
"""
 
# --------------------------
# SETS AVAILABLE MODEL PARAMETERS, THEIR TYPE, DEFAULT VALUE AND DESCRIPTION
# Change default values here

_PARAMETERS = [
        ("generations",int,5000,"number of generations per run"), 
        ("populations",int,1,"number of identical populations per run"), 
        ("L",int,5,"life time of each animal in time steps"), # default 5
        ("kd",float,0.02,"constant cost of plasticity"), #0.02
        ("ka",float,0.01,"cost of each adaptation"), #0.01
        ("tau",float,0.25,"coefficient of lifetime payoff exponential"), #0.25
        ("environment",float,[31.62,0.20,1,0,0], "parameters of each environment in the form R P A B O"),
        ("environment_name",str,"","displayed name of each environment"),
        ("size",int,5000,"Specifies number of animals in each environment"),  #5000     
        ("folder",str,"No_mut","Create additional folder to put output in"),
        ("format",str,"pdf","Format of the figures in timeseries (png or pdf)"),
        ("random_choice",bool,0,"If animals for cloning/killing should be chosen at random or dependent on fitness"),
        ("verbose",bool,0,"triggers verbose output to command line"),   
        ("desc",str,"","Description of the run appended to the path"),
        ("time_tag",bool,0,"Set current time (+description) as folder name"),
        ("save_all",bool,0,"Saves all animals' genes for each generation"),
        ("force_plast",bool,0,"Forces animals to use plastic strategy"),
        ("proc",int,1,"Number of processes (populations) that are executed at the same time"),
        ("std_min",float,[0,0.01,3,4],"Stop loop when, after the given number of generations, the desired standard deviation for the desired genes is reached; order: I0,I0p,W,a,b,bp,h,mu,s,sc,t"),
         ("lineage_stop",bool,0,"Stop if all animas are related to each other (common ancestor)"),


        ("mutation",float,[0.0,0.0,0.05,0.0],"initial mutation rate and scale of mutation steps, and their stds for mutation"), #0.001 0 0.05 0
        ("stop_mutation",int,0,"Stops mutation after the given number of generations, 0 disables the feature"),
        ("evolve_a_b",bool,0,"Lets a and b evolve also for non-plastic animals"),
        ("random_a_b",bool,0,"a and b are selected randomly after mutation to plasticity"),
        ("discrete_s",bool,0,"trait s (plasticity) is discrete (0 or 1)"),
        ("h_random",bool,0,"if h is distributed randomly at the beginning"),
        ("no_dbh",bool,1, "Disable DBH for constant runs by switching off mutation for h"),
        ("I0off",bool,1,"Set I0,I0' to 0 at the beginning"),
        ("plot_every",int,50,"detailed output is plotted every N generations (0 = never)"),
        ("hgt",bool,0,"If HGT is turned on"),
        ("check",bool,1,"If animals check the fitness of the donor before doing HGT"),
        ("kt",float,0,"cost of each horizontal gene transfer"),
        ("path",str,"","set path for genes to use, if empty: path.txt is used (for variable runs) or random animals are created (constant runs)"), 
        ("use_pop",int,1,"which of the populations to use for genes"),

#for variable runs: 
        ("q",float,2.9,"controls expected number of offspring in variable scenario"), #2.2
        ("trans",bool,1,"if true, use these (changed) constants, if false, use the ones from the file"),
        ("stop_half",bool,0,"Stop after half of the populations survived to save time"),
        ("survival_goal",float,0,"Goal for survival rate; Stop if too many populations died out already"),
        ("stop_below",float,[0,0],"For base extinction runs: Stop if population size falls below the given fraction of the original size after the given number of generations"),

]
# --------------------------

import argparse

class ModelConstants(dict):
    """Implement class for containing constants"""
    def __init__(self):
        super(ModelConstants,self).__init__()
        for param in _PARAMETERS:
            if param[1]==bool:
                self[param[0]] = bool(param[2])
            else:
                self[param[0]] = param[2]
 
    def __setattr__(self,name,value):
        raise Exception("Constants are read-only!")
            
    def change_constant(self,key,val):
    # ModelConstant instances' properties should only be changed through this method
        if key in self:                        
            self[key] = val
        else:
            
            raise KeyError("Key {0} is not a valid model constant identifier!".format(key))
    


# Parse parameters from command line
model_constants = ModelConstants()
parser = argparse.ArgumentParser()

for key in _PARAMETERS:
    if key[0] in ["environment"]: # Setting R,P,A,B,O 
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs=5,help=key[3])
    elif key[0] in ["mutation"]: # Setting mutation parameters
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs=4,help=key[3])
    elif key[0] in ["stop_below"]: 
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs=2,help=key[3])
    elif key[0] in ["std_min"]: 
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs="*",help=key[3])
    elif key[1]==bool: # Flags (true or false, no argument)
        parser.add_argument("--"+key[0],type=int,help=key[3])
    else: # Ordinary, single arguments (all optional)
        
        parser.add_argument("--"+key[0],type=key[1],help=key[3])

# Not included in _PARAMETERS, needs to be parsed outside of the loop
for key in ["R","P","A","B","O"]:
    parser.add_argument("--"+key,type=float,nargs="*",help="Overrides parameter {0} for each environment".format(key))

# Store all read arguments in a dict
args = parser.parse_args().__dict__
# Update model_constants object with read parameters
for key in _PARAMETERS:   
    if args[key[0]] or args[key[0]]==0:
        if key[0] in['environment','mutation','stop_below','std_min']:
            val=args[key[0]][0]
        elif key[1]==bool:
            val=bool(args[key[0]])
        else:      
            val=args[key[0]]
        model_constants.change_constant(key[0],val)
for i,key in enumerate(["R","P","A","B","O"]):
    environment = model_constants["environment"]
    if args[key]:
        environment[i] = args[key][0]
        model_constants.change_constant("environment",environment)

# Print some information
print("\nRunning model with the following parameters:")
for key in _PARAMETERS:
    print("\t{0}: {1}".format(key[0],model_constants[key[0]]))
print("\n")
