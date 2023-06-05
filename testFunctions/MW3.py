# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:58:02 2023

@author: r.dewinter
"""
import numpy as np

class MW3:
    def __init__(self):
        self.lower = np.array([0]*6)
        self.upper = np.array([1]*6)
        self.nConstraints = 2
        self.nObj = 2
        self.ref = np.array([1,7])
        self.nadir = np.array([])
        self.cheapConstr = [False,False]
        self.cheapObj = [False,False]
    
    def evaluate(self, x):
        g3 = 1 
        for i in range(2,len(x)):
            g3 += 2*(x[i] + (x[i-1]-0.5)**2 - 1)**2
        f1 = x[0]
        f2 = g3 * (1-f1/g3)
        l = np.sqrt(2)*f2 - np.sqrt(2)*f1
        c1 = 1.05 - f1 - f2 + 0.45*np.sin(0.75*np.pi*l)**6
        c1 = -1*c1
        c2 = 0.85 - f1 - f2 +0.3*np.sin(0.75*np.pi*l)**2
        
        return np.array([f1,f2]), np.array([c1,c2])
 
    def cheap_evaluate(self,x):
        g3 = 1 
        for i in range(2,len(x)):
            g3 += 2*(x[i] + (x[i-1]-0.5)**2 - 1)**2
        f1 = x[0]
        f2 = g3 * (1-f1/g3)
        l = np.sqrt(2)*f2 - np.sqrt(2)*f1
                        
        f_cheap = []
        if self.cheapObj[0]:
            f_cheap.append(f1)
        else:
            f_cheap.append(np.nan)
            
        if self.cheapObj[1]:
            f_cheap.append(f2)
        else:
            f_cheap.append(np.nan)      
        
        if self.cheapConstr[0]:
            c1 = 1.05 - f1 - f2 + 0.45*np.sin(0.75*np.pi*l)**6
            c1 = -1*c1
        else:
            c1 = np.nan
        if self.cheapConstr[1]:
            c2 = 0.85 - f1 - f2 +0.3*np.sin(0.75*np.pi*l)**2
        else:
            c2 = np.nan
        
        return [ np.array(f_cheap), np.array([c1, c2]) ]
    
problem = MW3()
X = np.random.rand(1000000,6)
obj, con = [], []
for x in X:
    res = problem.evaluate(x)
    obj.append(res[0])
    con.append(res[1])

cons = np.array(con)


