# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 09:19:42 2018

@author: r.dewinter
"""
import numpy as np
import time
import random

#def write_results(constraints, objectives):
#    try:
#        con = np.genfromtxt('conOSY.csv',delimiter=',')
#        obj = np.genfromtxt('objOSY.csv',delimiter=',')
#        con = np.asmatrix(con)
#        obj = np.asmatrix(obj)
#    except:
#        con = np.empty((0,6))
#        obj = np.empty((0,2))
#    con = np.append(con,constraints,axis=0)
#    obj = np.append(obj,objectives,axis=0)
#    np.savetxt('conOSY.csv',con,delimiter=',')
#    np.savetxt('objOSY.csv',obj,delimiter=',')

class OSY:
    def __init__(self):
        self.lower = np.array([0,0,1,0,1,0])
        self.upper = np.array([10,10,5,6,5,10])
        self.nConstraints = 6
        self.nObj = 2
        self.ref = np.array([0,386])
        self.nadir = np.array([-41.81, 76.00])
        self.cheapConstr = [False, False, False, False, False, False]
        self.cheapObj = [False, False]

    def evaluate(self, x):
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]
        x4 = x[3]
        x5 = x[4]
        x6 = x[5]
        
        f1 = -25*(x1-2)**2 - (x2-2)**2 - (x3-1)**2 - (x4-4)**2 - (x5-1)**2
        f2 = x1**2 + x2**2 + x3**2 + x4**2 + x5**2 + x6**2
        
        g1 = x1 + x2 - 2
        g2 = 6 - x1 - x2
        g3 = 2 - x2 + x1
        g4 = 2 - x1 + 3*x2
        g5 = 4 - (x3-3)**2 - x4
        g6 = (x5-3)**2 + x6 -4
        
        objectives = np.array([f1, f2])
        constraints = np.array([g1,g2,g3,g4,g5,g6])
        constraints = -1*constraints #transform for sacobra
        return objectives, constraints
    
    def cheap_evaluate(self, x):
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]
        x4 = x[3]
        x5 = x[4]
        x6 = x[5]

        if self.cheapObj[0]:
            f1 = -25*(x1-2)**2 - (x2-2)**2 - (x3-1)**2 - (x4-4)**2 - (x5-1)**2
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = x1**2 + x2**2 + x3**2 + x4**2 + x5**2 + x6**2
        else:
            f2 = np.nan

        if self.cheapConstr[0]:
            g1 = x1 + x2 - 2
        else:
            g1 = np.nan

        if self.cheapConstr[1]:
            g2 = 6 - x1 - x2
        else:
            g2 = np.nan
            
        if self.cheapConstr[2]:
            g3 = 2 - x2 + x1
        else:
            g3 = np.nan

        if self.cheapConstr[3]:
            g4 = 2 - x1 + 3*x2
        else:
            g4 = np.nan

        if self.cheapConstr[4]:
            g5 = 4 - (x3-3)**2 - x4
        else:
            g5 = np.nan

        if self.cheapConstr[5]:
            g6 = (x5-3)**2 + x6 -4
        else:
            g6 = np.nan
        
        objectives = np.array([f1, f2])
        constraints = np.array([g1,g2,g3,g4,g5,g6])
        constraints = -1*constraints #transform for sacobra
        return objectives, constraints


#iteration time 147.76699995994568
#37929.53300023079
#38003.78700017929
    
#iteration time 163.98299980163574
#14955.501999855042
#14997.977999925613
    
#iteration time 123.14800000190735
#13357.06500005722
#13410.796999931335