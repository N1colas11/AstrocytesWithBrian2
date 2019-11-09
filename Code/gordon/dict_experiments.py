#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 18:28:40 2019

@author: tnt19a
"""

"""
Figure 5: Astrocytes connected in networks

De Pitta model with multiple compartments, but without diffusion. 
All compartments are the same, but there is a soma. 
"""

# Experiment with a parameter dictionary

import matplotlib      # DELETE
# matplotlib.use('agg')  # DELETE
from brian2 import *
from brian2.units import *
from brian2.units.allunits import ymole, ymol, ymolar, umolar
#from collections import defaultdict

default_dict = { 
    'O_P' : 4.4*umolar/second,      # Maximal Ca^2+ uptake rate by SERCAs (Op) (0.9)
    'K_P' : 0.05 * umolar,          # Ca2+ affinity of SERCAs (Kp)
    'Omega_C' : 6./second,           # Maximal rate of Ca^2+ release by IP_3Rs
    'Omega_L' : 0.11/second         # Maximal rate of Ca^2+ leak from the ER (0.1)
}

print("default_d: ", default_dict)

    
# Set all default values

# I want a series of default parameters, which can be overwritten by values 
# in dictionary passed in an argument to the calilng function.

def setParams(**default_dict):
    print(default_dict)

    # One could access the dictionary directly but the notation is then heavier
    O_P = default_dict['O_P']
    K_P = default_dict['K_P']
    Omega_C = default_dict['Omega_C']
    Omega_L = default_dict['Omega_L']

setParams(O_P=1, K_P=2*second)

quit()
