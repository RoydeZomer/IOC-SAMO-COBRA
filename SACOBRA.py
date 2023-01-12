# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 10:25:12 2017

@author: r.dewinter
"""
import numpy as np

def verboseprint(verbose, important, message):
    if verbose != 0:
        if verbose==2 or (verbose==1 and important):
            print(message)

def standardize_obj(x):
    mean = np.mean(x)
    std = np.std(x)
    x = (x-mean)/std
    return x, mean, std

def reverseStandardize_obj(x, mean, std):
    x2 = x*std + mean
    return x2

def rescale_constr(x):
    maxG = np.max(x)
    minG = np.min(x)
    divider = maxG - minG
    x = x/divider
    return x, divider

def reverse_rescale_constr(x, divider):
    x = x*divider
    return x

def scaleRescale(xStart, originalL, originalU, newlower, newupper):
    return (newupper-newlower)*((xStart-originalL)/(originalU-originalL))+newlower

def rescaleWrapper(fn,lower,upper,newlower,newupper):
    oldfn = fn
    def newfn(x):
        x = scaleRescale(x,newlower,newupper,lower,upper)
        y = oldfn(x)
        return(y)
    return(newfn)
    
def plog(y, pShift=0.0):
    res = np.full_like(y,0.0,dtype='float')
    i = 0
    for yi in y:
        if yi-pShift>=0:
            res[i] = np.log(1+yi-pShift)
        else:
            res[i] = -np.log(1-(yi-pShift))
        i += 1
    return res

def plogReverse(y,pShift=0):
    if y>0:
        ret = np.exp(y)-1+pShift
    else:
        ret = pShift+1-(1/np.exp(y))
    return ret
  
def inverseRescale(x,cobra):
    z = scaleRescale(x, cobra['lower'], cobra['upper'], cobra['originalL'], cobra['originalU'])
    return(z)            