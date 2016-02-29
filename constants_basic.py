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
        #("not_opulation_size",int,5000,"number of animals per population"), #default 5000
        ("generations",int,500,"number of generations per run"), #default
            ("L",int,5,"life time of each animal in time steps"), # 5
        ("kd",float,0.02,"constant cost of plasticity"), #0.02
        ("ka",float,0.01,"cost of each adaptation"), #0.01
        ("tau",float,0.25,"coefficient of lifetime payoff exponential"), #0.25
        ("q",float,2.2,"controls expected number of offspring in variable scenario"), #2.2
        ("mu",float,0.001,"mutation rate of the genes"), #0.001
        ("environments",float,[[1E2,0.1,0.,1.0,0.],
                    [1E2,0.1,1.0,0.,0.],
                    [1E2,0.9,1.0,0.,0.],
                    [100,0.9,0.,1.0,0.]], "parameters of each environment "+ "in the form R P A B O"),
        ("environment_names",str,["Random, unpredictable",
                        "Deterministic, unpredictable",
                        "Deterministic, predictable",
                        "Random, predictable"],"displayed name of each environment"),
        ("environment_sizes",int,[1250,1250,1250,1250],"Specifies number of animals in each environment"),                
        ("km",float,0.2,"cost of migration"), #0.2
        ("limit",str,"m","names of genes that should be limited to [0,1]"),
        ("populations",int,1,"number of identical populations per run"), 
        ("plot_every",int,5,"detailed output is plotted every N generations (0 = never)"),
        ("verbose",bool,False,"triggers verbose output to command line"),   
        ("scaling",bool,False,"Decreases gene efficiency for extreme values, by introducing scaling function in the adaption process"),
        ("migration",bool,False,"Allow migration between environments. In constant mode the population is controlled as a whole. If false, environments are completely independent"),

]
# --------------------------

import argparse
import sys

if sys.argv[0] == "/Users/matthias/Documents/popdyn/botero-simultation-dion/main_variable1.py":
    _VARIABLE = True
else:
    _VARIABLE = False
#print(_VARIABLE)
# Required arguments when using main_variable.py
_VARIABLE_PARAMETERS = [
            ("mean_file",str,"./output/16-01-27_12-39-45/pop1_mean_genes.csv","specifies the path to the file containing means of the input genes"),
            ("std_file",str,"./output/16-01-27_12-39-45/pop1_std_genes.csv","specifies the path to the file containing the std. dev. of the input genes")
        ]

class ModelConstants(dict):
    """Implement class for containing constants"""
    def __init__(self):
        super(ModelConstants,self).__init__()
        for param in _PARAMETERS:
            self[param[0]] = param[2]
        if _VARIABLE:
            for param in _VARIABLE_PARAMETERS:
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
    elif key[0] in ["environment_names","limit"]: # May have arbitrary many arguments
        parser.add_argument("--"+key[0],type=key[1],action="append",nargs="*",help=key[3])
    elif key[0] in ["verbose"]: # Flags (true or false, no argument)
        parser.add_argument("--"+key[0],action="store_true",help=key[3])
    else: # Ordinary, single arguments (all optional)
        parser.add_argument("--"+key[0],type=key[1],help=key[3])

# Not includes in _PARAMETERS, needs to be parsed outside of the loop
for key in ["R","P","A","B","O"]:
    parser.add_argument("--"+key,type=int,nargs=Nenv,help="Overrides parameter {0} for each environment".format(key))

# Parse required arguments when variable breeding is used
if _VARIABLE:
    for key in _VARIABLE_PARAMETERS:
        parser.add_argument(key[0], type=key[1], help=key[3])

# Store all read arguments in a dict
args = parser.parse_args().__dict__

# Update model_constants object with read parameters
for key in _PARAMETERS:
    if args[key[0]]:
        model_constants.change_constant(key[0],args[key[0]])
if _VARIABLE:
    for key in _VARIABLE_PARAMETERS:
        if args[key[0]]:
            model_constants.change_constant(key[0],args[key[0]])
for i,key in enumerate(["R","P","A","B","O"]):
    environments = model_constants["environments"]
    if args[key]:
        for env in environments:
            env[i] = args[key][i]

# Print some information
print("\nRunning model with the following parameters:")
for key in _PARAMETERS:
    print("\t{0}: {1}".format(key[0],model_constants[key[0]]))
print("\n")
if _VARIABLE:
    print("Using input files for start genes:")
    for key in _VARIABLE_PARAMETERS:
        print("\t{0}: {1}".format(key[0],model_constants[key[0]]))
    print("\n")
