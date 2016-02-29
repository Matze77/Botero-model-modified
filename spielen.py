# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 09:41:01 2016

@author: matthias
"""
import numpy as np
import sys
sys.path.insert(0, './src')

#
# Import third-party packages
#

import numpy as np # For efficient array operations
import matplotlib.pyplot as plt # For plotting
import time # For timing parts of the script, optimizing run time
import pandas as pd # Easier data handling
import os # To create directories
import datetime # To access the current time
import sys # To access command line arguments
import warnings # To warn the user

try: # Seaborn makes prettier plots, but is not installed in a fresh Anaconda python
    import seaborn as sns 
    have_seaborn = True
except ImportError:
    have_seaborn = False

#
# Import other parts of the project
#

from animal import Animal
from population import Population
from environment import Environment
from constants import model_constants
from iterate_population import iterate_population
l=(1,0.6,0,0.9,0.6,0,0,0,0)
s=np.array(l)
k=Animal(s,1)
f=Animal(s,2)
p=k.mutate()
print(p)


print(k.position)
print(f.position)
