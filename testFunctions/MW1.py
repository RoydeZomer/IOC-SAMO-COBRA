# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 11:12:09 2021

@author: r.dewinter
"""
import numpy as np

# m is number of objectives
# p total number of constraints, j is jth constraint
#xII = subvector (x_m till Xn)

def MW1(x):
    n = len(x)
    g1 = 1 #+ np.sum((1-np.exp(-10*(x[2:]-0.5-(i-1)/(2*n))**2)))
    for i in range(2, n):
        g1 += np.sum(1-np.exp(-10*(x[i]-0.5-(i-1)/(2*n))**2))
        
        
    f1 = x[0]
    f2 = g1*(1-0.85*f1/g1)
    
    l = np.sqrt(2)*f2 - np.sqrt(2)*f1
    g1 = 1-f1-f2+0.5*np.sin(2*np.pi*l)**8
    
    return np.array([f1, f2]), -1*np.array([g1])

# import matplotlib.pyplot as plt
# rngMin = np.zeros(4)
# rngMax = np.ones(4)
# problemCall = MW1
# nconstraints = 1
# ref = np.array([1,2])

# par = len(rngMin)
# runs = 1000000
# objectives = np.empty((runs,len(ref)))
# constraints = np.empty((runs,nconstraints))
# for i in range(runs):
#     x = rngMax*np.random.rand(par)+rngMin
#     objectives[i],constraints[i] = problemCall(x)

# print(np.sum(np.sum(constraints<=0,axis=1)==nconstraints)/runs)

# objectives_feasible = objectives[np.sum(constraints<=0,axis=1)==1]
# plt.plot(objectives_feasible[:,0], objectives_feasible[:,1], 'ro')