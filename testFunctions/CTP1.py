# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:18:54 2018

@author: r.dewinter
"""

import numpy as np

class CTP1:
    def __init__(self):
        self.lower = np.array([0,0])
        self.upper = np.array([1,1])
        self.nConstraints = 2
        self.nObj = 2
        self.ref = np.array([1,2])
        self.nadir = np.array([0.99, 1.00])
        self.cheapConstr = [True, True]
        self.cheapObj = [False, False]

    def evaluate(self, x):
        x1 = x[0]
        x2 = x[1]
        
        f1 = x1
        f2 = (1+x2)*np.exp((-1*x1)/(1+x2))
        
        g1 = f2 / (0.858 * np.exp( -0.541 * f1 )) - 1
        g2 = f2 / (0.728 * np.exp( -0.285 * f1 )) - 1
        
        objectives = np.array([f1, f2])
        constraints = np.array([g1,g2])
        constraints = -1*constraints #transform for sacobra
        return np.array([objectives, constraints])    
    
    def cheap_evaluate(self, x):
        x1 = x[0]
        x2 = x[1]
        
        if self.cheapObj[0]:
            f1 = x1
        else:
            f2 = np.nan

        if self.cheapObj[1]:
            f2 = (1+x2)*np.exp((-1*x1)/(1+x2))
        else:
            f1 = np.nan

        if self.cheapConstr[0]:
            g1 = f2 / (0.858 * np.exp( -0.541 * f1 )) - 1
        else:
            g1 = np.nan

        if self.cheapConstr[1]:
            g2 = f2 / (0.728 * np.exp( -0.285 * f1 )) - 1
        else:
            g2 = np.nan

        objectives = np.array([f1, f2])
        constraints = np.array([g1,g2])
        return np.array([objectives, -1*constraints])    

#iteration time 80.75800013542175
#14392.611999988556
#14419.184000015259
    
#iteration time 26.826000213623047
#2550.3870000839233
#2563.927999973297
    
#iteration time 20.795000076293945
#2715.882000207901
#2727.6260001659393