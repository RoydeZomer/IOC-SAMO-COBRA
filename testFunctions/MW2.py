# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:49:00 2021

@author: r.dewinter
"""
import numpy as np

def MW2(x):
    n = len(x)
    g2 = 1
    for i in range(2,n):
        zi = 1 - np.exp(-10*((x[i]-(i-1)/n))**2)
        g2 += 1.5+(0.1/n)*zi**2 - 1.5*np.cos(2*np.pi*zi)
    
    f1 = x[0]
    f2 = g2*(1-f1/g2)
    
    l = np.sqrt(2)*f2 - np.sqrt(2)*f1
    c1 = 1 - f1 - f2 + 0.5*np.sin(3*np.pi*l)**8
    return np.array([f1,f2]), -1*np.array([c1])

# import matplotlib.pyplot as plt
# rngMin = np.zeros(4)
# rngMax = np.ones(4)
# problemCall = MW2
# nconstraints = 1
# ref = np.array([1,7])

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