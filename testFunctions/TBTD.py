# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 15:48:47 2018

@author: r.dewinter
"""
import numpy as np

class TBTD:
    def __init__(self):
        self.lower = np.array([1,0.0005,0.0005])
        self.upper = np.array([3,0.05,0.05])
        self.nConstraints = 3
        self.nObj = 2
        self.ref = np.array([0.1,50000])
        self.nadir = np.array([0.1, 92774.46])
        self.cheapConstr = [True, True, True]
        self.cheapObj = [False, False]

    def evaluate(self, x):
        y = x[0]
        x1 = x[1]
        x2 = x[2]
        
        fvolume = (x1*((16+y**2)**0.5)) + (x2*((1+y**2)**0.5))
        fstress = (20*((16+y**2)**0.5))/(y*x1)
        fstressBC = (80*((1+y**2)**0.5))/(y*x2)
        
        g1 = fvolume - 0.1
        g2 = fstress - 100000
        g3 = fstressBC - 100000
        
        return [np.array([fvolume, fstress]), np.array([g1,g2,g3])]
    
    def cheap_evaluate(self, x):
        y = x[0]
        x1 = x[1]
        x2 = x[2]
        fstressBC = (80*((1+y**2)**0.5))/(y*x2)
        
        if self.cheapObj[0]:
            fvolume = (x1*((16+y**2)**0.5)) + (x2*((1+y**2)**0.5))
        else:
            fvolume = np.nan

        if self.cheapObj[1]:
            fstress = (20*((16+y**2)**0.5))/(y*x1)
        else:
            fstress = np.nan

        if self.cheapConstr[0]:
            g1 = fvolume - 0.1
        else:
            g1 = np.nan
        
        if self.cheapConstr[1]:
            g2 = fstress - 100000
        else:
            g2 = np.nan
        
        if self.cheapConstr[2]:
            g3 = fstressBC - 100000
        else:
            g3 = np.nan
        
        return [np.array([fvolume, fstress]), np.array([g1,g2,g3])]

#rngMin = np.array([1,1e-16,1e-16])
#rngMax = np.array([3,1,1])
#nVar = 3
#parameters = np.empty((1000000,3))
#objectives = np.empty((1000000,2))
#constraints = np.empty((1000000,3))
#objectives[:] = 0
#constraints[:] = 0
#for i in range(1000000):
#    x = np.random.rand(nVar)*(rngMax-rngMin)+rngMin
#    parameters[i] = x
#    obj, cons = TBTD(x)
#    objectives[i] = obj
#    constraints[i] = cons
#
#a = np.sum(constraints<0, axis=1)==3
#sum(a)
    

#iteration time 146.92799997329712
#15948.508999824524
#15992.766000270844
    
#iteration time 162.1710000038147
#13681.631999969482
#13729.1859998703