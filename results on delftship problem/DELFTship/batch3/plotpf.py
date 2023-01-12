# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 12:33:30 2022

@author: r.dewinter
"""

import json 
import matplotlib.pyplot as plt 
import numpy as np

cobra = json.load(open('cobra2.json'))
obj = np.array(cobra['Fres'])
pf = np.array(cobra['paretofrontier'])
Gres = np.array(cobra['Gres'])
feasible = np.sum(Gres<=0, axis=1)==4

plt.title('Obtained Solutions for i-SAMO-COBRA Batch=3')
plt.plot(obj[:,0]*-1,obj[:,1]*-1,'r.', label='Infeasible')
plt.plot(obj[feasible][:,0]*-1,obj[feasible][:,1]*-1,'bo', label='Feasible')
plt.plot(pf[:,0]*-1,pf[:,1]*-1,'go', label='Pareto Efficient')
plt.xlabel('Attained Index [-]')
plt.ylabel('Cargo Hold [$m^{3}$]')
plt.legend()
plt.savefig('paretofrontier.pdf')
plt.show()

hv_progress_cobra = cobra['hypervolumeProgress']
plt.plot(hv_progress_cobra, '-', color='#D7191C', label='SAMO-COBRA')
plt.title('Hypervolume Convergence Plot')
plt.xlabel('Evaluations')
plt.ylabel('Hypervolume')
plt.legend()
plt.savefig('convergencePlot.pdf')
plt.show()