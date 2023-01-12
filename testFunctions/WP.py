# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 15:05:43 2018

@author: r.dewinter
"""
import numpy as np

class WP:
    def __init__(self):
        self.lower = np.array([0.01,    0.01,  0.01])
        self.upper = np.array([0.45,    0.1,  0.1])
        self.nConstraints = 7   
        self.nObj = 5
        self.ref = np.array([83000, 1350, 2.85, 15989825, 25000])
        self.nadir = np.array([75684, 1347, 2.85, 7625734, 24896])
        self.cheapConstr = [True, True, True, True, True, True, True]
        self.cheapObj = [False, False, False, False, False]

    def evaluate(self,x):
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]
        f1 = 106780.37 * (x2+x3) + 61704.67
        f2 = 3000 * x1
        f3 = 30570 * 0.02289 * x2 / (0.06*2289)**0.65
        f4 = 250 * 2289 * np.exp(-39.75*x2+9.9*x3 + 2.74)
        f5 = 25*((1.39/(x1*x2)) + 4940.0*x3 - 80.0)
        g1 = 0.00139/(x1*x2) + 4.94*x3 - 0.08 - 1
        g2 = 0.000306/(x1*x2)+1.082*x3 - 0.0986 - 1
        g3 = 12.307/(x1*x2) + 49408.24*x3 + 4051.02 - 50000
        g4 = 2.098/(x1*x2) + 8046.33*x3 - 696.71 - 16000
        g5 = 2.138/(x1*x2) + 7883.39*x3 - 705.04 - 10000
        g6 = 0.417*(x1*x2) + 1721.26*x3 - 136.54 - 2000
        g7 = 0.164/(x1*x2) + 631.13*x3 - 54.48 - 550
            
        return np.array([f1,f2,f3,f4,f5]), np.array([g1,g2,g3,g4,g5,g6,g7])
        
    def cheap_evaluate(self, x):
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]
        
        if self.cheapObj[0]:
            f1 = 106780.37 * (x2+x3) + 61704.67
        else:
            f1 = np.nan

        if self.cheapObj[1]:
            f2 = 3000 * x1
        else:
            f2 = np.nan
        
        if self.cheapObj[2]:
            f3 = 30570 * 0.02289 * x2 / (0.06*2289)**0.65 
        else:
            f3 = np.nan
        
        if self.cheapObj[3]:
            f4 = 250 * 2289 * np.exp(-39.75*x2+9.9*x3 + 2.74) 
        else:
            f4 = np.nan
        
        if self.cheapObj[4]:
            f5 = 25*((1.39/(x1*x2)) + 4940.0*x3 - 80.0)
        else:
            f5 = np.nan
            
        if self.cheapConstr[0]:
            g1 = 0.00139/(x1*x2) + 4.94*x3 - 0.08 - 1
        else:
            g1 = np.nan
        
        if self.cheapConstr[1]:
            g2 = 0.000306/(x1*x2)+1.082*x3 - 0.0986 - 1
        else:
            g2 = np.nan
        
        if self.cheapConstr[2]:
            g3 = 12.307/(x1*x2) + 49408.24*x3 + 4051.02 - 50000
        else:
            g3 = np.nan
        
        if self.cheapConstr[3]:
            g4 = 2.098/(x1*x2) + 8046.33*x3 - 696.71 - 16000
        else:
            g4 = np.nan
            
        if self.cheapConstr[4]:
            g5 = 2.138/(x1*x2) + 7883.39*x3 - 705.04 - 10000
        else:
            g5 = np.nan
            
        if self.cheapConstr[5]:
            g6 = 0.417*(x1*x2) + 1721.26*x3 - 136.54 - 2000
        else:
            g6 = np.nan

        if self.cheapConstr[6]:
            g7 = 0.164/(x1*x2) + 631.13*x3 - 54.48 - 550
        else:
            g7 = np.nan

        return np.array([f1,f2,f3,f4,f5]), np.array([g1,g2,g3,g4,g5,g6,g7])