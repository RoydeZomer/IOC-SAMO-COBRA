# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 16:46:10 2017

@author: r.dewinter
"""
import itertools
import numpy as np
import multiprocessing as mp
import time
from functools import partial

from SACOBRA import scaleRescale
from SACOBRA import rescaleWrapper
from SACOBRA import standardize_obj
from SACOBRA import rescale_constr
from SACOBRA import plog
from lhs import lhs 
from halton import halton
from transformLHS import transformLHS
from paretofrontFeasible import paretofrontFeasible
from hypervolume import hypervolume


def newfn(x, fn=None, newlower=None, newupper=None, lower=None, upper=None):
    x = scaleRescale(x,newlower,newupper,lower,upper)
    y = fn(x)
    return(y)

def cheap_SAMO_COBRA_Init(problem, nCores=7, computeStartingPoints=16,
              feval=None,
              infillCriteria="PHV",
              batch=1,
              oneShot=False,
              iterPlot=False,
              initDesign="HALTON",
              initDesPoints=None,
              seqTol=1e-6,
              epsilonInit=None,
              epsilonMax=None,
              epsilonLearningRate=0.1,
              surrogateUpdateLearningRate=0.1,
              newlower=-1,
              newupper=1,
              cobraSeed=1):
    
    start = time.time()
    if hasattr(problem, 'nConstraints'):
        nConstraints = problem.nConstraints
    else:
        raise ValueError('Problem class must at least have the number of constraints defined in nConstraints ref object')
        
    if hasattr(problem, 'ref'):
        ref = problem.ref
    else:
        raise ValueError('Problem class must at least have a reference point defined in ref object')
    
    if hasattr(problem, 'lower'):
        originalL = problem.lower
    else:
        raise ValueError('Problem class must at least have lower bound for decision variables in lower object')
    
    if hasattr(problem, 'upper'):
        originalU = problem.upper
    else:
        raise ValueError('Problem class must at least have upper boudn for decision variables in upper object')
              
    dimension = len(originalL)
    nObj = len(ref)
    print('start pSAMO_COBRA with seed',cobraSeed)
    if initDesPoints is None:
        initDesPoints = max(batch,dimension+1)
    if feval is None:
        feval = dimension*40
    if oneShot:
        initDesPoints = int(np.ceil(feval/2))
        batch = initDesPoints

    originalL = np.array(originalL)
    originalU = np.array(originalU)    
    
    seqFeval = (dimension*batch*batch)*50
    if computeStartingPoints is None:
        computeStartingPoints = (dimension)*batch
    
    phase = 'init'
    
    dimension = dimension #number of parameters
    lower = np.array([newlower]*dimension)
    upper = np.array([newupper]*dimension)
    
    l = newlower-newupper
    if epsilonInit is None:
        epsilonInit = [0.02*l]*nConstraints
    if epsilonMax is None:
        epsilonMax = [0.04*l]*nConstraints
        
    if initDesPoints>=feval:
        raise ValueError('feval should be larger then initial sample size')
        
    I = np.empty((1,1))
    I[:] = np.NaN
    Gres = np.empty((1,1))
    Gres[:] = np.NaN
    Fres = []
    
    np.random.seed(cobraSeed)
    if initDesign == 'RANDOM':
        I = np.random.uniform(low=lower,high=upper,size=(initDesPoints,dimension))
           
    elif initDesign =='LHS':
        I = lhs(dimension, samples=initDesPoints, criterion="center", iterations=5)
        I = transformLHS(I, lower, upper)

    elif initDesign == 'HALTON':
        I = halton(dimension, initDesPoints)
        I = scaleRescale(I, 0, 1, lower, upper)
        
    elif initDesign == 'BOUNDARIES':
        inputdata = [[newlower,newupper]] * dimension
        result = np.array(list(itertools.product(*inputdata)))
        resultLength = len(result)
        indicator = [True]*(initDesPoints) + [False]*(resultLength - initDesPoints)
        np.random.shuffle(indicator)
        I = result[indicator]
        
    else:
        raise ValueError('not yet implemented or invalid init design')
    
    if hasattr(problem, 'evaluate'):
        fnn = partial(newfn, fn=problem.evaluate, newlower=lower, newupper=upper, lower=originalL, upper=originalU)#rescaleWrapper(fn,originalL,originalU,newlower,newupper)
    else:
        raise ValueError('Problem class must at least have evaluate function')
        
    if hasattr(problem, 'cheap_evaluate'):
        fnnCheap = partial(newfn, fn=problem.cheap_evaluate, newlower=lower, newupper=upper, lower=originalL, upper=originalU)#rescaleWrapper(fn,originalL,originalU,newlower,newupper)
        resF, resG = fnnCheap(lower)
        cheapObj = ~np.isnan(resF)
        cheapCon = ~np.isnan(resG)
    else:
        fnnCheap = None
        cheapObj = np.array([False]*nObj)
        cheapCon = np.array([False]*nConstraints)
    
    if nCores>1:
        parallel = nCores
        pool = mp.Pool(processes=parallel)
    else:
        pool = None
    Fres, Gres, res_time = randomResultsFactory(I,fnn,nConstraints,nObj,pool,batch)

    hypervolumeProgress = np.empty((len(Fres),1))
    hypervolumeProgress[:] = np.NAN
    for i in range(len(Fres)):
        paretoOptimal = np.array([False]*(len(Fres)))
        paretoOptimal[:i+1] = paretofrontFeasible(Fres[:i+1,:],Gres[:i+1,:])
        paretoFront = Fres[paretoOptimal]
        hypervolumeProgress[i] = [hypervolume(paretoFront, ref)]
    
    numViol = np.sum(Gres>0,axis=1)
    
    maxViol = np.max([np.zeros(len(Gres)),np.max(Gres, axis=1)],axis=0)
    
    A = I #contains all evaluated points
    
    pff = paretofrontFeasible(Fres,Gres)
    pf = Fres[pff]
    hv = hypervolumeProgress[-1].item()
    
    FresStandardized = np.full_like(Fres, 0)
    FresStandardizedMean = np.zeros(nObj)
    FresStandardizedStd = np.zeros(nObj)
    FresPlogStandardized = np.full_like(Fres, 0)
    FresPlogStandardizedMean = np.zeros(nObj)
    FresPlogStandardizedStd = np.zeros(nObj)
    for obji in range(nObj):
        res, mean, std = standardize_obj(Fres[:,obji])        
        FresStandardized[:,obji] = res
        FresStandardizedMean[obji] = mean 
        FresStandardizedStd[obji] = std
        
        plogFres = plog(Fres[:,obji])
        res, mean, std = standardize_obj(plogFres)
        FresPlogStandardized[:,obji] = res
        FresPlogStandardizedMean[obji] = mean 
        FresPlogStandardizedStd[obji] = std
    
    GresRescaled = np.full_like(Gres, 0)
    GresRescaledDivider = np.zeros(nConstraints)
    GresPlogRescaled = np.full_like(Gres, 0)
    GresPlogRescaledDivider = np.zeros(nConstraints)
    for coni in range(nConstraints):
        GresRescaled[:,coni], GresRescaledDivider[coni] = rescale_constr(Gres[:,coni])
        plogGres, GresPlogRescaledDivider[coni] = rescale_constr(Gres[:,coni])
        GresPlogRescaled[:,coni] = plog(plogGres)


    end = time.time()
    cobra = dict()
    cobra['ref'] = ref
    cobra['nObj'] = nObj
    cobra['currentHV'] = hv
    cobra['hypervolumeProgress'] = hypervolumeProgress
    cobra['paretoFrontier'] = pf
    cobra['paretoFrontierStandardized'] = FresStandardized[pff]
    cobra['paretoFrontierFeasible'] = pff
    cobra['fn'] = fnn
    cobra['fnCheap'] = fnnCheap
    cobra['cheapObj'] = cheapObj
    cobra['cheapCon'] = cheapCon
    cobra['batch'] = batch
    cobra['oneShot'] = oneShot
    cobra['pool'] = pool
    cobra['dimension'] = dimension
    cobra['nConstraints'] = nConstraints
    cobra['nVar'] = len(lower)
    cobra['lower'] = lower
    cobra['upper'] = upper
    cobra['initDesPoints'] = initDesPoints
    cobra['feval'] = feval
    cobra['A'] = A
    cobra['Fres'] = Fres
    cobra['FresStandardized'] = FresStandardized
    cobra['FresStandardizedMean'] = FresStandardizedMean
    cobra['FresStandardizedStd'] = FresStandardizedStd
    cobra['FresPlogStandardized'] = FresPlogStandardized
    cobra['FresPlogStandardizedMean'] = FresPlogStandardizedMean 
    cobra['FresPlogStandardizedStd'] = FresPlogStandardizedStd
    cobra['Gres'] = Gres
    cobra['GresRescaled'] = GresRescaled
    cobra['GresRescaledDivider'] = GresRescaledDivider
    cobra['GresPlogRescaled'] = GresPlogRescaled
    cobra['GresPlogRescaledDivider'] = GresPlogRescaledDivider
    cobra['numViol'] = numViol
    cobra['maxViol'] = maxViol
    cobra['epsilonInit'] = epsilonInit
    cobra['epsilonMax'] = epsilonMax
    cobra['epsilonLearningRate'] = epsilonLearningRate
    cobra['ptail'] = True
    cobra['seqFeval'] = seqFeval
    cobra['computeStartPointsStrategy'] = 'multirandom'
    cobra['computeStartingPoints'] = computeStartingPoints
    cobra['surrogateUpdateLearningRate'] = surrogateUpdateLearningRate
    cobra['seqTol'] = seqTol
    cobra['cobraSeed'] = cobraSeed
    cobra['RBFmodel'] = ['CUBIC','GAUSSIAN','MULTIQUADRIC','INVQUADRIC','INVMULTIQUADRIC','THINPLATESPLINE']
    cobra['bestPredictor'] = [{'objKernel':[cobra['RBFmodel'][0]]*nObj, 'objLogStr': ['Standardized']*nObj, 'conKernel':[cobra['RBFmodel'][0]]*nConstraints, 'conLogStr':['Rescaled']*nConstraints}]
    cobra['phase'] = [phase]*initDesPoints
    cobra['plot'] = iterPlot
    cobra['optimizationTime'] = np.array([(end-start)])
    cobra['evaluationTime'] = res_time
    
    if cobra['oneShot']:
        cobra['EPS'] = [0.0]*cobra['nConstraints']
    else:
        cobra['EPS'] = epsilonInit

    
    if infillCriteria == "PHV" or infillCriteria == "SMS":
        cobra['infillCriteria'] = infillCriteria # "PHV" or "SMS"
    else:
        raise ValueError("This infill criteria is not implemented")
        
    surrogateErrors = {}
    for kernel in cobra['RBFmodel']:
        for obji in range(cobra['nObj']):
            surrogateErrors['OBJ'+str(obji)+kernel] = [0]*cobra['initDesPoints']
            surrogateErrors['OBJ'+str(obji)+'PLOG'+kernel] = [0]*cobra['initDesPoints']
        for coni in range(cobra['nConstraints']):
            surrogateErrors['CON'+str(coni)+kernel] = [0]*cobra['initDesPoints']
            surrogateErrors['CON'+str(coni)+'PLOG'+kernel] = [0]*cobra['initDesPoints']
    cobra['SurrogateErrors'] = surrogateErrors
    
    return(cobra)

def randomResultsFactory(I, fn, nConstraints, nObj, pool, batch):
    objs = np.empty((len(I),nObj))
    constr = np.empty((len(I),nConstraints))
    objs[:] = np.NaN
    constr[:] = np.NaN
    res_time = []
    if pool is None:
        res = []
        for row in I:
            stt = time.time()
            res.append(fn(row))
            res_time.append(time.time()-stt)
    else:
        res = []
        for i in range(int(len(I)/batch)):    
            stt = time.time()
            res += pool.map(fn, I[i*batch:(i+1)*batch])
            res_time.append(time.time()-stt)
    
    i = 0
    for result in res:
        objectiveScore, constaintScore = result
        objs[i,:] = objectiveScore
        constr[i,:] = constaintScore
        i += 1
        
    return objs, constr, res_time
        


            
