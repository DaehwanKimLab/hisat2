#!/usr/bin/env python

import sys, os
import math
import random
from copy import deepcopy


"""
"""
def clone_IMGTHLA_database():
    os.system("git clone https://github.com/jrob119/IMGTHLA.git")
    
    # Check out one particular revision just to have the same data across multiple computers        
    revision = "d3b559b34b96ff9e7f0d97476222d8e4cdee63ad" # Revision on November 16, 2016
    # revision = "45c377516bdb7f1b926" # Revision on July 14, 2016
    os.system("cd IMGTHLA; git checkout %s; cd .." % revision)

