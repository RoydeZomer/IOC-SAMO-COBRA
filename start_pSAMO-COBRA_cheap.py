# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 15:32:48 2021

@author: r.dewinter
"""
import numpy as np
import multiprocessing
from multiprocessing import freeze_support
import matplotlib.pyplot as plt

from hypervolume import hypervolume
from cheap_SAMO_COBRA_Init import cheap_SAMO_COBRA_Init
from cheap_SAMO_COBRA_PhaseII import cheap_SAMO_COBRA_PhaseII

from testFunctions.BNH import BNH


# example test function with two variables, two constraints, two objectives


if __name__ == '__main__':  
    freeze_support() 
    seed = 1
    batch = 5
    nCores = multiprocessing.cpu_count()
    
    problem = BNH()
    problem.cheapConstr = [True]*len(problem.cheapConstr) 
    cobra = cheap_SAMO_COBRA_Init(problem, batch=batch, nCores=nCores, computeStartingPoints=16, cobraSeed=seed, iterPlot=False) 
    cobra = cheap_SAMO_COBRA_PhaseII(cobra)