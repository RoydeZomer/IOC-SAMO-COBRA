# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:58:02 2023

@author: r.dewinter
"""
import numpy as np

class MW11:
    def __init__(self):
        self.lower = np.array([0]*6) 
        self.upper = np.array([2**0.5]*6) 
        self.nConstraints = 4
        self.nObj = 2
        self.ref = np.array([30]*self.nObj)
        self.nadir = np.array([])
        self.cheapConstr = [False,False,False,False]
        self.cheapObj = [False,False]
    
    def evaluate(self, x):
        g3 = 1 
        for i in range(2,len(x)):
            g3 += 2*(x[i] + (x[i-1]-0.5)**2 - 1)**2
        f1 = g3*x[0]
        
        p = (f1/g3)**2
        if p > 2: # due to numerical instabilities p can become 2.000000000000001
            p = 2
        f2 = g3 * (2-p)**(1/2)
        
        c1 = (3-f1**2 - f2)*(3-2*f1**2-f2)
        c1 = c1*-1
        
        c2 = (3-0.625*f1**2 - f2) * (3-7*f1**2 - f2)
        
        c3 = (1.62 - 0.18*f1**2 - f2) * (1.125 - 0.125*f1**2 - f2)
        c3 = c3*-1

        c4 = (2.07-0.23*f1**2 -f2) * (0.63-0.07*f1**2 - f2)
        
        return np.array([f1,f2]), np.array([c1,c2,c3,c4])
 
    def cheap_evaluate(self,x):
        g3 = 1 
        for i in range(2,len(x)):
            g3 += 2*(x[i] + (x[i-1]-0.5)**2 - 1)**2
        f1 = g3*x[0]
        
        p = (f1/g3)**2
        if p > 2: # due to numerical instabilities p can become 2.000000000000001
            p = 2
        f2 = g3 * (2-p)**(1/2)
                     
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
            c1 = (3-f1**2 - f2)*(3-2*f1**2-f2)
            c1 = c1*-1
        else:
            c1 = np.nan
            
        if self.cheapConstr[1]:
            c2 = (3-0.625*f1**2 - f2) * (3-7*f1**2 - f2)
        else:
            c2 = np.nan
        
        if self.cheapConstr[2]:
            c3 = (1.62 - 0.18*f1**2 - f2) * (1.125 - 0.125*f1**2 - f2)
            c3 = c3*-1
        else:
            c3 = np.nan
        
        if self.cheapConstr[3]:
            c4 = (2.07-0.23*f1**2 -f2) * (0.63-0.07*f1**2 - f2)
        else:
            c4 = np.nan
        
        return [ np.array(f_cheap), np.array([c1,c2,c3,c4]) ]


# problem = MW11()
# X = np.random.rand(1000000,6)
# obj, con = [], []
# for x in X:
#     res = problem.evaluate(x)
#     obj.append(res[0])
#     con.append(res[1])    
# con = np.array(con)
# sum(np.sum(con<=0,axis=1)==4)/1000000*100
# np.max(obj, axis=0)

