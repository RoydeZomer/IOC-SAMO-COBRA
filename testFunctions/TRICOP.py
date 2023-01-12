# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:15:45 2020

@author: r.dewinter
"""


import numpy as np

class TRICOP:
    def __init__(self):
        self.lower = np.array([-4,-4])
        self.upper = np.array([4,4])
        self.nConstraints = 3 
        self.nObj = 3
        self.ref = np.array([34,-4,90])
        self.nadir = np.array([7.67, -11.77, 25.91])
        self.cheapConstr = [False, False, False]
        self.cheapObj = [False, False, False]

    def evaluate(self, x):
        f1 = 0.5*(x[0]-2)**2 + 1/13*(x[1]+1)**2 + 3
        f2 = 1/175*(x[0]+x[1]-3)**2 + 1/17*(2*x[1]-x[0])**2 - 13
        f3 = 1/8*(3*x[0]-2*x[1]+4)**2 + 1/27*(x[0]-x[1]+1)**2 + 15
        c1 = 4- 4*x[0] - x[1]
        c2 = x[0] + 1
        c3 = -x[0] + x[1] + 2
        #-1* constr because of sacobra's constraint handling
        return [ np.array([f1,f2,f3]), -1*np.array([c1,c2,c3]) ]
    
    def cheap_evaluate(self, x):
        if self.cheapObj[0]:
            f1 = 0.5*(x[0]-2)**2 + 1/13*(x[1]+1)**2 + 3
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = 1/175*(x[0]+x[1]-3)**2 + 1/17*(2*x[1]-x[0])**2 - 13
        else:
            f2 = np.nan
            
        if self.cheapObj[2]:
            f3 = 1/8*(3*x[0]-2*x[1]+4)**2 + 1/27*(x[0]-x[1]+1)**2 + 15
        else:
            f3 = np.nan
            
        if self.cheapConstr[0]:
            c1 = 4- 4*x[0] - x[1]
        else:
            c1 = np.nan

        if self.cheapConstr[1]:
            c2 = x[0] + 1
        else:
            c2 = np.nan

        if self.cheapConstr[2]:
            c3 = -x[0] + x[1] + 2
        else:
            c3 = np.nan

        return [ np.array([f1,f2,f3]), -1*np.array([c1,c2,c3]) ]


# amount = 1000000
# x = np.random.rand(amount*2)
# x = np.reshape(x, (amount, 2))*8 - 4
# objs = np.zeros((amount,3))
# cons = np.zeros((amount,3))
# for i in range(len(x)):
#     objs[i], cons[i] = TRICOP(x[i])
    
    
# sum(np.sum(cons<=0,axis=1)==3)/1000000

# import matplotlib.pyplot as plt
# plt.plot(objs[:,1], objs[:,2], 'ro')
# plt.plot(cobra['paretoFrontier'][:,0],cobra['paretoFrontier'][:,1],'go')
