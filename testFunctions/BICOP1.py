# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:13:45 2020

@author: r.dewinter
"""
import numpy as np

class BICOP1:
    def __init__(self):
        self.lower = np.zeros(10)
        self.upper = np.ones(10)
        self.nConstraints = 1        
        self.nObj = 2
        self.ref = np.array([9,9])     
        self.nadir = np.array([1, 1])
        self.cheapConstr = [True]
        self.cheapObj = [False, False]

    def evaluate(self, x):
        x = np.array(x)
        g1 = 1 + 9*np.sum(x[1:]/9)
    
        f1 = x[0]*g1
        f2 = g1 - np.sqrt(f1/g1)
        
        c1 = g1
        #-1* constr because of sacobra's constraint handling
        return [ np.array([f1, f2]), -1*np.array([c1]) ]
    
    def cheap_evaluate(self, x):
        x = np.array(x)
        g1 = 1 + 9*np.sum(x[1:]/9)

        if self.cheapObj[0]:
            f1 = x[0]*g1
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = g1 - np.sqrt(f1/g1)
        else:
            f2 = np.nan

        if self.cheapConstr[0]:
            c1 = g1
        else:
            c1 = np.nan

        return [ np.array([f1, f2]), -1*np.array([c1]) ]

# amount = 1000000
# x = np.random.rand(amount*10)
# x = np.reshape(x, (amount, 10))
# objs = np.zeros((amount,2))
# cons = np.zeros((amount,1))
# for i in range(len(x)):
#     objs[i], cons[i] = BICOP1(x[i])
    
# import matplotlib.pyplot as plt
# plt.plot(objs[:,0], objs[:,1], 'ro')
