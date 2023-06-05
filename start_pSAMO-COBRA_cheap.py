# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 15:32:48 2021

@author: r.dewinter
"""
import multiprocessing
from multiprocessing import freeze_support
import numpy as np

from cheap_SAMO_COBRA_Init import cheap_SAMO_COBRA_Init
from cheap_SAMO_COBRA_PhaseII import cheap_SAMO_COBRA_PhaseII
from hypervolume import hypervolume

from testFunctions.BNH import BNH

# example test function with two variables, two constraints, two objectives


if __name__ == '__main__':  
    freeze_support() 
    batch = 5
    nCores = 7
    for seed in range(1):
        problem = BNH()
        problem.cheapConstr = [True]*len(problem.cheapConstr)
        problem.cheapObj = [False]*problem.nObj
        cobra = cheap_SAMO_COBRA_Init(problem, batch=batch, nCores=nCores, computeStartingPoints=16, cobraSeed=seed, iterPlot=True) 
        cobra = cheap_SAMO_COBRA_PhaseII(cobra)
        