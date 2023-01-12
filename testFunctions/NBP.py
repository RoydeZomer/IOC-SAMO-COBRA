# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 13:12:23 2020

@author: r.dewinter
"""
import numpy as np

class NBP:
    def __init__(self):
        self.lower = np.array([20, 10])
        self.upper = np.array([250, 50])
        self.nConstraints = 5
        self.nObj = 2
        self.ref = np.array([11150, 12500])
        self.nadir = np.array([12500, 114.09])
        self.cheapConstr = [True, True, True, True, True]
        self.cheapObj = [False, False]

    def evaluate(self, x):
        h = x[0]
        b = x[1]
        
        l = 1500
        F = 5000
        E = 216620 #GPa
        v = 0.27
        G = 86650
        
        A = b*h #cross sectional area of the beam
        sigma = 6*F*l/(b*h**2) # bending stress
        
        delta = (F*(l**3)) / (3*E*((b*h**3)/12)) # maximum tip deflection
        tau = 3*F/(2*b*h) #maximum allowable shear stress
        hb = h/b #height to breadth ratio
        
        F_crit = -1*(4/(l**2)) * ( G*(((b*(h**3))/12) + (((b**3)*h)/12)) * E*(((b**3)*h)/12) /(1-(v**2)) )**0.5 #failure force of buckling
    
        
        
        f1 = A
        f2 = sigma
        
        
        c1 =  delta - 5 #delta <= 5
        c2 =  sigma - 240 #sigma_B <= sigma_y
        c3 = tau - 120 #tau <= (sigma_y/2)
        c4 = hb - 10 #hb<= 10 #height to breadth ratio
        c5 = F_crit + 2*F  #F_crit >= f*F
        
        return [ np.array([f1,f2]), np.array([c1,c2,c3,c4,c5]) ]
    def cheap_evaluate(self, x):
        h = x[0]
        b = x[1]
        
        l = 1500
        F = 5000
        E = 216620 #GPa
        v = 0.27
        G = 86650
        
        A = b*h #cross sectional area of the beam
        sigma = 6*F*l/(b*h**2) # bending stress
        
        delta = (F*(l**3)) / (3*E*((b*h**3)/12)) # maximum tip deflection
        tau = 3*F/(2*b*h) #maximum allowable shear stress
        hb = h/b #height to breadth ratio
        
        F_crit = -1*(4/(l**2)) * ( G*(((b*(h**3))/12) + (((b**3)*h)/12)) * E*(((b**3)*h)/12) /(1-(v**2)) )**0.5 #failure force of buckling
        
        if self.cheapObj[0]:
            f1 = A
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = sigma
        else:
            f2 = np.nan

        if self.cheapConstr[0]:
            c1 =  delta - 5 #delta <= 5
        else:
            c1 = np.nan
            
        if self.cheapConstr[1]:
            c2 =  sigma - 240 #sigma_B <= sigma_y
        else:
            c2 = np.nan

        if self.cheapConstr[2]:
            c3 = tau - 120 #tau <= (sigma_y/2)
        else:
            c3 = np.nan

        if self.cheapConstr[3]:
            c4 = hb - 10 #hb<= 10 #height to breadth ratio
        else:
            c4 = np.nan

        if self.cheapConstr[4]:
            c5 = F_crit + 2*F  #F_crit >= f*F
        else:
            c5 = np.nan


        return [ np.array([f1,f2]), np.array([c1,c2,c3,c4,c5]) ]

# import matplotlib.pyplot as plt
# rngMin = np.array([20, 10])
# rngMax = np.array([250, 50])
# nVar = 2
# ref = np.array([100,100])
# parameters = np.empty((100000,2))
# objectives = np.empty((100000,2))
# constraints = np.empty((100000,5))
# objectives[:] = 0
# constraints[:] = 0
# parameters[:]= 0
# for i in range(100000):
#     x = np.random.rand(nVar)*(rngMax-rngMin)+rngMin
#     parameters[i] = x
#     obj, cons = NBP(x)
#     objectives[i] = obj
#     constraints[i] = cons

# a = np.sum(constraints<=0, axis=1)==5
# print(np.sum(constraints<=0, axis=0))
# #sum(a)

# # plt.plot(objectives[:,0],objectives[:,1],'ro')
# plt.plot(objectives[a][:,0],objectives[a][:,1],'go')
# plt.show()
# plt.close()
# plt.plot(parameters[~a][:,1],parameters[~a][:,0],'ro')
# plt.plot(parameters[a][:,1],parameters[a][:,0],'go')

