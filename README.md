# IOC-SAMO-COBRA
Efficient Multi-Objective Optimization for Problems with (In)expensive Objectives and Constraints

Expensive objectives and constraints are a key characteristic of real-world multi-objective optimization problems. In practice they often occur jointly with inexpensive objectives and constraints. This repository presents the Self-Adapting Multi-Objective Constraint Optimization algorithm that uses Radial Basis function Approximations for the expensive objectives and constraints and the inexpensive objectives and constraints directly. In this repository you can find the algorithm code, and the supplementary material with all the raw results. 

Here an example on how to use the algorithm on the BNH test problem:
```python
import numpy as np
from multiprocessing import freeze_support

from pSAMO_COBRA_Init import pSAMO_COBRA_Init
from pSAMO_COBRA_PhaseII import pSAMO_COBRA_PhaseII

# example test function class two variables, two constraints, two objectives
class BNH:
    def __init__(self):
        self.lower = np.array([0,0]) # lower bound of the decision variables
        self.upper = np.array([5,3]) # upper bound of the decision variables
        self.nConstraints = 2 # number of constraints
        self.nObj = 2 # number of objectives
        self.ref = np.array([140,50]) # Maximum objective score per objective you are interested in. 
        self.nadir = np.array([136,49.24]) # optional 
        self.cheapConstr = [True,True] # constraint 1 and constraint two are considered inexpensive and are directly used
        self.cheapObj = [False,False] # both objectives are considered as expensive and will be predicted with RBFs. 
        
    # expensive evaluation method, note that SAMO-COBRA assumes both objectives are to be minimized.
    def evaluate(self,x): 
        f1 = 4*x[0]**2+4*x[1]**2
        f2 = (x[0]-5)**2 + (x[1]-5)**2
        
        c1 = -1*((x[0]-5)**2 + x[1]-25)
        c2 = (x[0]-8)**2 + (x[1]-3)**2 - 7.7
        
        return [ np.array([f1, f2]), -1*np.array([c1,c2]) ]
    
    # in case some of the objectives or considered are cheap, we use this evaluation method during the search for candidate solutions.
    def cheap_evaluate(self,x):
        if self.cheapObj[0]:
            f1 = 4*x[0]**2+4*x[1]**2
        else:
            f1 = np.nan
        if self.cheapObj[1]:
            f2 = (x[0]-5)**2 + (x[1]-5)**2
        else:
            f2 = np.nan
        
        if self.cheapConstr[0]:
            c1 = -1*((x[0]-5)**2 + x[1]-25)
        else:
            c1 = np.nan
        if self.cheapConstr[1]:
            c2 = (x[0]-8)**2 + (x[1]-3)**2 - 7.7
        else:
            c2 = np.nan
        
        return [ np.array([f1, f2]), -1*np.array([c1,c2]) ]


if __name__ == '__main__':  
    freeze_support() # required for multiprocessing on windows machines
    # The algorithm starts searching new candidate solutions from random locations in the search space.
    seed = 1 Influence the random lacations by setting the seed
    batch = 5 # how many candidate solutions do you want to evaluate per iteration? should be larger or equal to 1. 
    nCores = multiprocessing.cpu_count() # define how many cores you want to use for the optimization process. 
    
    problem = BNH() # test problem class
    problem.cheapConstr = [True]*len(problem.cheapConstr) # example where the constraints are considered as cheap
    # initialize the algorithm and run the DoE
    cobra = cheap_SAMO_COBRA_Init(problem, batch=batch, nCores=nCores, computeStartingPoints=16, cobraSeed=seed, iterPlot=True) 
    # Start optimizing with IOC-SAMO-COBRA
    cobra = cheap_SAMO_COBRA_PhaseII(cobra)

```
Results on all test functions in table format can be found in the following two figures:
![alt text](https://github.com/RoydeZomer/IOC-SAMO-COBRA/blob/main/suplementary%20material/allresultsHV.PNG?raw=true)
![alt text](https://github.com/RoydeZomer/IOC-SAMO-COBRA/blob/main/suplementary%20material/allresultsIGD.PNG?raw=true)
