# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:49:00 2021

@author: r.dewinter
"""
import numpy as np

class MW2:
    def __init__(self):
        self.lower = np.array([0]*6) 
        self.upper = np.array([1]*6)
        self.nConstraints = 1
        self.nObj = 2
        self.ref = np.array([1,7])
        self.nadir = np.array([1,1])
        self.cheapConstr = [False]
        self.cheapObj = [False, False]
    
    def evaluate(self,x):
        n = len(x)
        g2 = 1
        for i in range(2,n):
            zi = 1 - np.exp(-10*((x[i]-(i-1)/n))**2)
            g2 += 1.5+(0.1/n)*zi**2 - 1.5*np.cos(2*np.pi*zi)
        
        f1 = x[0]
        f2 = g2*(1-f1/g2)
        
        l = np.sqrt(2)*f2 - np.sqrt(2)*f1
        c1 = 1 - f1 - f2 + 0.5*np.sin(3*np.pi*l)**8
        return [np.array([f1,f2]), -1*np.array([c1])]
    
    def cheap_evaluate(self,x):
        n = len(x)
        g2 = 1
        for i in range(2,n):
            zi = 1 - np.exp(-10*((x[i]-(i-1)/n))**2)
            g2 += 1.5+(0.1/n)*zi**2 - 1.5*np.cos(2*np.pi*zi)
        
        f1 = x[0]
        f2 = g2*(1-f1/g2)
        
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
            c1 = -1*(1 - f1 - f2 + 0.5*np.sin(3*np.pi*l)**8)
        else:
            c1 = np.nan
                    
        return [ np.array(f_cheap), np.array([c1]) ]


# problem = MW2()
# X = np.random.rand(1000000,6)
# obj, con = [], []
# for x in X:
#     res = problem.evaluate(x)
#     obj.append(res[0])
#     con.append(res[1])    
# con = np.array(con)
# sum(con<=0)/1000000*100









