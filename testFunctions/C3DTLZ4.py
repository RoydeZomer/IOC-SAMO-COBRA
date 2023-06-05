# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:50:49 2018

@author: r.dewinter
"""

import numpy as np

class C3DTLZ4:
    
    def __init__(self):
        self.lower = np.array([0,0,0,0,0,0])
        self.upper = np.array([1,1,1,1,1,1])
        self.nConstraints = 2
        self.nObj = 2 
        self.ref = np.array([3,3])
        self.nadir = np.array([2.00, 2.00])
        self.cheapConstr = [False,False]
        self.cheapObj = [False,False]
    
    def evaluate(self, x):
        x = np.array(x)
        gx = np.sum((x[1:]-0.5)**2)
        
        
        f1 = (1+(gx))*np.cos(x[0]*np.pi/2)
        f2 = (1+(gx))*np.sin(x[0]*np.pi/2)
        
        c1 = (f1**2)/4 + f2**2 - 1
        c2 = (f2**2)/4 + f1**2 - 1
        #-1* constr because of sacobra's constraint handling
        return [ np.array([f1, f2]), -1*np.array([c1,c2]) ]

    def cheap_evaluate(self, x):
        x = np.array(x)
        gx = np.sum((x[1:]-0.5)**2)
        f11 = (1+(gx))*np.cos(x[0]*np.pi/2)
        f22 = (1+(gx))*np.sin(x[0]*np.pi/2)
        
        if self.cheapObj[0]:
            f1 = f11
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = f22
        else:
            f2 = np.nan

        if self.cheapConstr[0]:
            c1 = (f11**2)/4 + f22**2 - 1
        else:
            c1 = np.nan
            
        if self.cheapConstr[1]:
            c2 = (f22**2)/4 + f11**2 - 1
        else:
            c2 = np.nan

        #-1* constr because of sacobra's constraint handling
        return [ np.array([f1, f2]), -1*np.array([c1,c2]) ]