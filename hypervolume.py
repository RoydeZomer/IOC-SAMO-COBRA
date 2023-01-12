# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 12:17:20 2018

@author: r.dewinter
"""

from pygmo import hypervolume as hvv
from numpy import all

def hypervolume(pointset, ref):
    """Compute the absolute hypervolume of a *pointset* according to the
    reference point *ref*.
    """
    
    pointset = pointset[all(pointset<=ref,axis=1)]
    
    if len(pointset)==0:
        return 0
    hv = hvv(pointset)
    contribution = hv.compute(ref)
    return contribution