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
        ("generations",int,100,"number of generations per run"), #default
        ("L",int,5,"life time of each animal in time steps"), # 5
        ("kd",float,0.02,"constant cost of plasticity"), #0.02
        ("ka",float,0.01,"cost of each adaptation"), #0.01
        ("tau",float,0.25,"coefficient of lifetime payoff exponential"), #0.25
        ("q",float,2.2,"controls expected number of offspring in variable scenario"), #2.2
        ("mutation",str,["normal","0.01","0.0","0.0"],"take initial mutation rate from normal (with mean value and std) or uniform (min and max) distribution; \
        the last value is the standard deviation for the mutation of this gene"), #0.001
        ("environment",float,[100,0.1,1,0,0], "parameters of each environment "+ "in the form R P A B O"),
        ("environment_name",str,"","displayed name of each environment"),
        ("size",int,5000,"Specifies number of animals in each environment"),                
        ("populations",int,4,"number of identical populations per run"), 
        ("plot_every",int,100,"detailed output is plotted every N generations (0 = never)"),
        ("verbose",bool,0,"triggers verbose output to command line"),   
        ("random_choice",bool,1,"If animals for cloning/killing should be chosen at random or dependent on fitness"),
        ("desc",str,"RP","Description of the run appended to the path"),
        ("time_tag",bool,0,"Set current time (+description) as folder name"),
        ("force_plast",bool,0,"Forces animals to use plastic strategy"),
        ("save_all",bool,0,"Saves all animals' genes for each generation"),
        ("proc",int,4,"Number of processes (populations) that are executed at the same time"),
        ("format",str,"pdf","Format of the figures in timeseries (png or pdf)"),
        ("folder",str,"","Create additional folder to put output in"),
        ("hgt",bool,0,"If HGT is turned on"),
        ("check",bool,0,"If animals check the fitness of the donor before doing HGT"),
        ("kh",float,0.02,"constant cost of hgt"),
        ("kt",float,0.01,"cost of each transfer"),
#for variable runs: 
        ("trans",bool,1,"if true, use these (changed) constants, if false, use the ones from the file"),
        ("path",str,"/Users/matthias/Documents/popdyn/botero-model/Output_to_analyze/botero_compare/R0.32_P0.50/","set path for genes to use, if empty: path.txt is used"),
        ("use_pop",int,1,"which of the populations to use for mean_genes"),
        ("stop_half",bool,1,"Stop after half of the populations survived to save time"),

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
    if key[0] in ["environment"]: # Setting R,P,A,B,O for each environment
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs=5,help=key[3])
    elif key[0] in ["mutation"]: # Setting R,P,A,B,O for each environment
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs=4,help=key[3])
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
        if key[0] in['environment',"mutation"]:
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
