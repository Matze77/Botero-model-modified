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
        ("generations",int,20000,"number of generations per run"), #default
            ("L",int,5,"life time of each animal in time steps"), # 5
        ("kd",float,0.02,"constant cost of plasticity"), #0.02
        ("ka",float,0.01,"cost of each adaptation"), #0.01
        ("tau",float,0.25,"coefficient of lifetime payoff exponential"), #0.25
        ("q",float,2.2,"controls expected number of offspring in variable scenario"), #2.2
        ("mu",float,0.001,"mutation rate of the genes"), #0.001
        ("environments",float,[1,1,0,0,0], "parameters of each environment "+ "in the form R P A B O"),
        ("environment_names",str,"","displayed name of each environment"),
        ("environment_sizes",int,5000,"Specifies number of animals in each environment"),                
        ("km",float,0.2,"cost of migration"), #0.2
        ("limit",str,["m","ma","h","a","s"],"names of genes that should be limited to [0,1]"),
        ("populations",int,1,"number of identical populations per run"), 
        ("plot_every",int,100,"detailed output is plotted every N generations (0 = never)"),
        ("verbose",bool,False,"triggers verbose output to command line"),   
        ("random_choice",bool,False,"If animals for cloning/killing should be chosen at random or dependent on fitness"),
        ("std_min",float,[],"Stop loop when desired standard deviation for the genes I0,a,b,h (for each environment) is reached"),
        ("lineage_stop",bool,False,"Stop if all animas are related to each other (common ancestor)"),
        ("desc",str,"","Description of the run appended to the path"),
#for variable runs: 
        ("trans",bool,False,"if true, use these (changed) constants, if false, use the ones from the file"),
        ("path",str,"","set path for genes to use, if empty: path.txt is used"),
        ("use_pop",int,1,"which of the populations to use for mean_genes"),

]
# --------------------------

import argparse

class ModelConstants(dict):
    """Implement class for containing constants"""
    def __init__(self):
        super(ModelConstants,self).__init__()
        for param in _PARAMETERS:
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
    if key[0] in ["environments"]: # Setting R,P,A,B,O for each environment
        Nenv = len(key[2])
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs=5,help=key[3])
    elif key[0] in ["environment_names","limit","ka","kd","mu","q","tau","environment_sizes","std_min"]: # May have arbitrary many arguments
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs="*",help=key[3])
    elif key[0] in ["verbose","scaling","migration","trans","random_choice","lineage_stop"]: # Flags (true or false, no argument)
        parser.add_argument("--"+key[0],action="store_true",help=key[3])
    else: # Ordinary, single arguments (all optional)
        parser.add_argument("--"+key[0],type=key[1],help=key[3])

# Not included in _PARAMETERS, needs to be parsed outside of the loop
for key in ["R","P","A","B","O"]:
    parser.add_argument("--"+key,type=float,nargs="*",help="Overrides parameter {0} for each environment".format(key))

# Store all read arguments in a dict
args = parser.parse_args().__dict__
# Update model_constants object with read parameters
for key in _PARAMETERS:
    if args[key[0]]:
        print(key[0])
        if key[0]=='environments':
            val=args[key[0]][0][0]
        else:
            try:  #for keys with single numerical entry          
                val=args[key[0]][0]
                if type(val)==str:
                    raise
            except: #for keys with several entries
                val=args[key[0]]
        model_constants.change_constant(key[0],val)
for i,key in enumerate(["R","P","A","B","O"]):
    environments = model_constants["environments"]
    if args[key]:
        environments[i] = args[key]
        model_constants.change_constant("environments",environments)

# Print some information
print("\nRunning model with the following parameters:")
for key in _PARAMETERS:
    print("\t{0}: {1}".format(key[0],model_constants[key[0]]))
print("\n")
