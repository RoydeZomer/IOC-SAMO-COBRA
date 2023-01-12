# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:56:36 2020

@author: r.dewinter
"""


import numpy as np    

class BICOP2:
    def __init__(self):
        self.lower = np.zeros(10)
        self.upper = np.ones(10)
        self.nConstraints = 2        
        self.nObj = 2
        self.ref = np.array([70,70]) 
        self.nadir = np.array([1.1, 1.08])
        self.cheapConstr = [True, True]
        self.cheapObj = [False, False]

    def evaluate(self, x):
        a = 0.1
        b = 0.9
        
        gx = np.sum(x[1:] - np.sin(0.5*np.pi*x[0]))**2
        
        f1 = x[0] + gx
        f2 = 1 - x[0]**2 + gx
        
        g1 = gx - a
        g2 = b - gx
        
        c1 = g1
        c2 = g2
        #-1* constr because of sacobra's constraint handling
        return [ np.array([f1, f2]), -1*np.array([c1,c2]) ]    
    
    def cheap_evaluate(self, x):
        a = 0.1
        b = 0.9
        
        gx = np.sum(x[1:] - np.sin(0.5*np.pi*x[0]))**2

        if self.cheapObj[0]:
            f1 = x[0] + gx
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = 1 - x[0]**2 + gx
        else:
            f2 = np.nan

        if self.cheapConstr[0]:
            g1 = gx - a
        else:
            g1 = np.nan
        
        if self.cheapConstr[1]:
            g2 = b - gx
        else:
            g2 = np.nan

        #-1* constr because of sacobra's constraint handling
        return [ np.array([f1, f2]), -1*np.array([g1,g2]) ]


# amount = 1000000
# x = np.random.rand(amount*10)
# x = np.reshape(x, (amount, 10))
# objs = np.zeros((amount,2))
# cons = np.zeros((amount,2))
# for i in range(len(x)):
#     objs[i], cons[i] = BICOP2(x[i])
    
# import matplotlib.pyplot as plt
# plt.plot(objs[:,0], objs[:,1], 'ro')
# plt.plot(objs[np.sum(cons<=0,axis=1)==2][:,0], objs[np.sum(cons<=0,axis=1)==2][:,1], 'bo')