# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:15:41 2021

@author: r.dewinter
"""
import numpy as np 

def MW3(x):
    g3 = 1 
    for i in range(2,len(x)):
        g3 += 2*(x[i] + (x[i-1]-0.5)**2 - 1)**2
    f1 = x[0]
    f2 = g3 * (1-f1/g3)
    l = np.sqrt(2)*f2 - np.sqrt(2)*f1
    c1 = 1.05 - f1 - f2 + 0.45*np.sin(0.75*np.pi*l)**6
    c2 = 0.85 - f1 - f2 +0.3*np.sin(0.75*np.pi*l)**2
    
    return np.array([f1,f2]), np.array([-1*c1,c2])

# import matplotlib.pyplot as plt
# rngMin = np.zeros(4)
# rngMax = np.ones(4)
# problemCall = MW3
# nconstraints = 2
# ref = np.array([1,3])

# par = len(rngMin)
# runs = 1000000
# objectives = np.empty((runs,len(ref)))
# constraints = np.empty((runs,nconstraints))
# for i in range(runs):
#     x = rngMax*np.random.rand(par)+rngMin
#     objectives[i],constraints[i] = problemCall(x)

# print(np.sum(np.sum(constraints<=0,axis=1)==nconstraints)/runs)

# objectives_feasible = objectives[np.sum(constraints<=0,axis=1)==nconstraints]
# plt.plot(objectives_feasible[:,0], objectives_feasible[:,1], 'ro')

