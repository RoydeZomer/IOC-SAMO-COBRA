# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:48:42 2017

@author: r.dewinter
"""

from SACOBRA import plog
from SACOBRA import plogReverse
from SACOBRA import standardize_obj
from SACOBRA import rescale_constr

from RbfInter import trainRBF
from RbfInter import interpRBF
from RbfInter import distLine

from lhs import lhs 
from hypervolume import hypervolume
from paretofrontFeasible import paretofrontFeasible
from visualiseParetoFront import visualiseParetoFront

from scipy import optimize
import numpy as np
import time
import warnings
import copy
import os
from functools import partial
import json
import pygmo as pg

def make_json_serializable(cobra):
    for key in cobra:
        if type(cobra[key]) is np.ndarray or type(cobra[key]) is np.matrix:
            cobra[key] = cobra[key].tolist()
    cobra['fn'] = None
    cobra['fnCheap'] = None
    cobra['pool'] = None
    return cobra

def getConstraintPrediction(x, surrogateModels, bestPredictor, nConstraints, GresPlogRescaledDivider, GresRescaledDivider, fnCheap, cheapCon, EPS=None):
    constraintPredictions = np.zeros(nConstraints)
    if any(cheapCon):
        _ , cheap_constr = fnCheap(x)
    for coni in range(nConstraints):
        if cheapCon[coni]:
            if EPS is None:
                constraintPredictions[coni] = cheap_constr[coni]
            else:
                constraintPredictions[coni] = cheap_constr[coni]+EPS[coni]**2
        else:
            conKernel = bestPredictor['conKernel'][coni]
            conLogStr = bestPredictor['conLogStr'][coni]
            surrogate = surrogateModels[conKernel]['Constraints'][conLogStr][coni]
            
            if conLogStr=='PLOGRescaled':
                constraintPrediction = interpRBF(np.array(x), surrogate)
                constraintPrediction = plogReverse(constraintPrediction)
                constraintPrediction = GresPlogRescaledDivider[coni] * constraintPrediction
            else:
                constraintPrediction = interpRBF(np.array(x), surrogate)
                constraintPrediction = GresRescaledDivider[coni] * constraintPrediction                
            
            if EPS is None:
                constraintPredictions[coni] = constraintPrediction
            else:
                constraintPredictions[coni] = constraintPrediction+EPS[coni]**2

    return constraintPredictions

def getEquallyImportantConstraintPrediction(x, surrogateModels, bestPredictor, nConstraints, GresRescaledDivider, fnCheap, cheapCon, EPS=None):
    constraintPredictions = []
    if any(cheapCon):
        _ , cheap_constr = fnCheap(x)
    for coni in range(nConstraints):
        if cheapCon[coni]:
            if EPS is None:
                constraintPredictions.append(-1*(cheap_constr[coni]/GresRescaledDivider[coni]))
            else:
                constraintPredictions.append(-1*((cheap_constr[coni]/GresRescaledDivider[coni]) + EPS[coni]**2))
        else:
            conKernel = bestPredictor['conKernel'][coni]
            conLogStr = bestPredictor['conLogStr'][coni]
            surrogate = surrogateModels[conKernel]['Constraints'][conLogStr][coni]
            
            constraintPrediction = interpRBF(np.array(x), surrogate)
            if conLogStr=='PLOGRescaled':
                constraintPrediction = plogReverse(constraintPrediction)
            
            if EPS is None:
                constraintPredictions.append(-1*(constraintPrediction))
            else:
                constraintPredictions.append(-1*(constraintPrediction+EPS[coni]**2))

    return constraintPredictions

def gCOBRA(x, A, lower, upper, surrogateModels, bestPredictor, nConstraints, GresRescaledDivider, fnCheap, cheapCon, EPS):
    h = 0
    distance = distLine(x, A)
    if any(distance<=(len(A[0])/1e4)):
        if min(distance) != 0:
            h = min(distance)**-2
        else:
            h = np.finfo(np.float64).max
    if not all(np.isfinite(distance)):
        h = np.finfo(np.float64).max
        
    constraintPrediction = getEquallyImportantConstraintPrediction(x, surrogateModels, bestPredictor, nConstraints, GresRescaledDivider,  fnCheap, cheapCon, EPS)
    
    if np.any(np.isnan(constraintPrediction)):
        warnings.warn('gCOBRA: constraintPrediction value is NaN, returning Inf',DeprecationWarning)
        return([np.finfo(np.float64).min]*(len(lower)*2+len(constraintPrediction)+1))
    
    boundaries = []
    for i in range(len(lower)):
        boundaries.append(x[i] - lower[i])
    for i in range(len(upper)):
        boundaries.append(upper[i]- x[i])
    
    h = [-1*h]
    h = h + constraintPrediction
    h = h + boundaries
    return(h)

def batch_gCOBRA(x, batch, A, lower, upper, surrogateModels, bestPredictor, nConstraints, GresRescaledDivider, fnCheap, cheapCon, EPS=None):
    h = []
    size = int(len(x)/batch)
    for i in range(batch):
        xi = x[i*size:(i+1)*size]
        h += gCOBRA(xi, A, lower, upper, surrogateModels, bestPredictor, nConstraints, GresRescaledDivider, fnCheap, cheapCon, EPS)
    return(h)    

def get_potentialSolution(x, surrogateModels, bestPredictor, nObj, infillCriteria, FresPlogStandardizedStd, FresStandardizedStd, FresPlogStandardizedMean, FresStandardizedMean, fnCheap, cheapObj):
    potentialSolution = np.zeros(nObj)
    if any(cheapObj):
        cheap_obj, _ = fnCheap(x)
        
    for obji in range(nObj):
        if cheapObj[obji]:
            potentialSolution[obji] = cheap_obj[obji]
        else:
            objKernel = bestPredictor['objKernel'][obji]
            objLogStr = bestPredictor['objLogStr'][obji]
            surrogate = surrogateModels[objKernel]['Objectives'][objLogStr][obji]
            
            potsol = 0
            uncertainty = 0
            if objLogStr=='PLOGStandardized':
                if infillCriteria == 'PHV':
                    potsol = interpRBF(x, surrogate, uncertainty=False)                
                    potsol = potsol*FresPlogStandardizedStd[obji] + FresPlogStandardizedMean[obji]
                    potsol = plogReverse(potsol)
                elif infillCriteria == 'SMS':
                    potsol, uncertainty = interpRBF(x, surrogate, uncertainty=True)
                    uncertainty = uncertainty * FresStandardizedStd[obji]
                    potsol = potsol*FresPlogStandardizedStd[obji] + FresPlogStandardizedMean[obji]
                    potsol = plogReverse(potsol)
                else:
                    raise ValueError("This infill criteria is not implemented")
            else:
                if infillCriteria == 'PHV':
                    potsol = interpRBF(x, surrogate, uncertainty=False)
                    potsol = potsol*FresStandardizedStd[obji] + FresStandardizedMean[obji]
                elif infillCriteria == 'SMS':
                    potsol, uncertainty = interpRBF(x, surrogate, uncertainty=True)
                    potsol = potsol*FresStandardizedStd[obji] + FresStandardizedMean[obji]
                    uncertainty = uncertainty * FresStandardizedStd[obji]
                else:
                    raise ValueError("This infill criteria is not implemented")
            
            potentialSolution[obji] = potsol - np.abs(uncertainty)
    return potentialSolution
    
def compute_infill_criteria_score(x, surrogateModels, bestPredictor, nObj, infillCriteria, FresPlogStandardizedStd, FresStandardizedStd, FresPlogStandardizedMean, FresStandardizedMean, currentHV, paretoFrontier, ref, fnCheap, cheapObj):
    if np.any(np.isnan(x)):
        return np.finfo(np.float64).max
    if not all(np.isfinite(x)):
        return np.finfo(np.float64).max

    potentialSolution = get_potentialSolution(x, surrogateModels, bestPredictor, nObj, infillCriteria, FresPlogStandardizedStd, FresStandardizedStd, FresPlogStandardizedMean, FresStandardizedMean, fnCheap, cheapObj)
    if not all(np.isfinite(potentialSolution)):
        return np.finfo(np.float64).max
    
    penalty = 0
    ##### add epsilon?
    logicBool = np.all(paretoFrontier<= potentialSolution, axis=1)
    for j in range(paretoFrontier.shape[0]):
        if logicBool[j]:
            p = - 1 + np.prod(1 + (potentialSolution-paretoFrontier[j,:]))
            penalty = max(penalty, p)
    if penalty == 0: #non-dominated solutions
        potentialFrontier = np.append(paretoFrontier, [potentialSolution], axis=0)
        myhv = hypervolume(potentialFrontier, ref)
        f = currentHV - myhv
    else:
        f = currentHV + penalty
    return f

def batch_infill_criteria_score(x, batch, surrogateModels, bestPredictor, nObj, infillCriteria, FresPlogStandardizedStd, FresStandardizedStd, FresPlogStandardizedMean, FresStandardizedMean, paretoFrontier, ref, fnCheap, cheapObj):
    if np.any(np.isnan(x)):
        return np.finfo(np.float64).max
    if not all(np.isfinite(x)):
        return np.finfo(np.float64).max
    
    
    size = int(len(x)/batch)
    solutions = np.zeros((batch,nObj))
    for i in range(batch):
        xi = x[i*size:(i+1)*size]
    
        potentialSolution = get_potentialSolution(xi, surrogateModels, bestPredictor, nObj, infillCriteria, FresPlogStandardizedStd, FresStandardizedStd, FresPlogStandardizedMean, FresStandardizedMean, fnCheap, cheapObj)
        
        if not all(np.isfinite(potentialSolution)):
            return np.finfo(np.float64).max
        
        solutions[i] = potentialSolution

    penalties = []        
    
    ##### add epsilon?
    tempPF = np.vstack((paretoFrontier,solutions))
    ndf, _, _, _ = pg.fast_non_dominated_sorting(tempPF)
    pf_indicator = np.zeros(len(tempPF))==1
    pf_indicator[ndf[0]] = True
    tempPF_pf = tempPF[pf_indicator]
    tempPF_dominated = tempPF[~pf_indicator]
    
    for potentialSolution in tempPF_dominated:
        if potentialSolution not in paretoFrontier:
            penalty = 0    
            logicBool = np.all(tempPF_pf<= potentialSolution, axis=1)
            for j in range(tempPF_pf.shape[0]):
                if logicBool[j]:
                    p = - 1 + np.prod(1 + (potentialSolution-tempPF_pf[j,:]))
                    penalty = max(penalty, p)
            penalties.append(penalty)
    
    myhv = hypervolume(tempPF_pf, ref)
    f = -1*myhv + sum(penalties)
    return f

def pool_job(xStart, criteria_function=None, cons=None, seqFeval=None, seqTol=None):
    opts = {'maxiter':seqFeval, 'tol':seqTol} 
    subMin = optimize.minimize(criteria_function, xStart, constraints=cons, options=opts, method='COBYLA')
    return subMin

def trainSurrogate(kernel, A, GresRescaled, GresPlogRescaled, FresStandardized, FresPlogStandardized, cheapConstr, cheapObj):
    surrogateModel = {'Constraints':{}, 'Objectives':{}, 'type':kernel}
    surrogateModel['Constraints'] = {'PLOGRescaled':[], 'Rescaled':[]}
    surrogateModel['Objectives'] = {'PLOGStandardized':[], 'Standardized':[]}
    
    gi = 0
    for g in GresRescaled:
        if not cheapConstr[gi]:
            surrogateModel['Constraints']['Rescaled'].append(trainRBF(A,g,ptail=True,squares=True,smooth=0.00,rbftype=kernel))
        else:
            surrogateModel['Constraints']['Rescaled'].append(None)
        gi += 1
        
    gi = 0
    for g in GresPlogRescaled:
        if not cheapConstr[gi]:
            surrogateModel['Constraints']['PLOGRescaled'].append(trainRBF(A,g,ptail=True,squares=True,smooth=0.00,rbftype=kernel))
        else:
            surrogateModel['Constraints']['PLOGRescaled'].append(None)
        gi += 1
    
    fi = 0        
    for f in FresStandardized:
        if not cheapObj[fi]:
            surrogateModel['Objectives']['Standardized'].append(trainRBF(A,f,ptail=True,squares=True,smooth=0.00,rbftype=kernel))
        else:
            surrogateModel['Objectives']['Standardized'].append(None)
        fi += 1
    
    fi = 0
    for f in FresPlogStandardized:
        if not cheapObj[fi]:
            surrogateModel['Objectives']['PLOGStandardized'].append(trainRBF(A,f,ptail=True,squares=True,smooth=0.00,rbftype=kernel))
        else:
            surrogateModel['Objectives']['PLOGStandardized'].append(None)
        fi += 1
    return surrogateModel

def cheap_SAMO_COBRA_PhaseII(cobra):
    # print("PHASE II started")
    phase = 'PHASE II'
    pool = cobra['pool']

    if cobra['hypervolumeProgress'] is None:
        raise ValueError("cobraPhaseII: cobra['hypervolumeProgress'] is None! First run smscobraInit")
        
    fn = cobra['fn']
    fnCheap = cobra['fnCheap']
    
    n = len(cobra['A'])
    if n==cobra['initDesPoints']:
        predHV = np.empty(cobra['initDesPoints'])
        predHV[:] = np.nan # structure to store surrogate optimization results
        cobra['optimizerConvergence'] = np.ones(cobra['initDesPoints']) # vector to store optimizer convergence
        feval = np.empty(cobra['initDesPoints'])
        feval[:] = np.nan
        
    if n >= cobra['feval']:
        raise ValueError("ERROR! Number of function evaluations after initialization is larger than total allowed evaluations")
    
    
    def kfold_best_predictor(cobra, pool):
        folds = 10
        surrogateErrors = {}
        for kernel in cobra['RBFmodel']:
            for obji in range(cobra['nObj']):
                surrogateErrors['OBJ'+str(obji)+kernel] = []
                surrogateErrors['OBJ'+str(obji)+'PLOG'+kernel] = []
            for coni in range(cobra['nConstraints']):
                surrogateErrors['CON'+str(coni)+kernel] = []
                surrogateErrors['CON'+str(coni)+'PLOG'+kernel] = []

        for foldi in range(folds):
            trainCobra = {'RBFmodel' : cobra['RBFmodel']}
            trainIndicator = np.array(range(len(cobra['A'])))
            trainCobra['A'] = cobra['A'][trainIndicator%folds!=foldi]
            trainCobra['FresStandardized'] = cobra['FresStandardized'][trainIndicator%folds!=foldi]
            trainCobra['FresPlogStandardized'] = cobra['FresPlogStandardized'][trainIndicator%folds!=foldi]
            trainCobra['GresRescaled'] = cobra['GresRescaled'][trainIndicator%folds!=foldi]
            trainCobra['GresPlogRescaled'] = cobra['GresPlogRescaled'][trainIndicator%folds!=foldi]
            trainCobra['cheapCon'] = cobra['cheapCon']
            trainCobra['cheapObj'] = cobra['cheapObj']
            
            surrogateModels = trainSurrogates(trainCobra, pool)
            
            testCobra = {'RBFmodel' : cobra['RBFmodel']}
            testCobra['A'] = cobra['A'][trainIndicator%folds==foldi]
            testCobra['FresStandardized'] = cobra['FresStandardized'][trainIndicator%folds==foldi]
            testCobra['FresPlogStandardized'] = cobra['FresPlogStandardized'][trainIndicator%folds==foldi]
            testCobra['GresRescaled'] = cobra['GresRescaled'][trainIndicator%folds==foldi]
            testCobra['GresPlogRescaled'] = cobra['GresPlogRescaled'][trainIndicator%folds==foldi]
            testCobra['Fres'] = cobra['Fres'][trainIndicator%folds==foldi]
            testCobra['Gres'] = cobra['Gres'][trainIndicator%folds==foldi]
            
            
            for i in range(len(testCobra['A'])): 
                for obji in range(cobra['nObj']):
                    for kernel in cobra['RBFmodel']:
                        yTrue = testCobra['Fres'][i]
                        conTrue = testCobra['Gres'][i]
                        x = testCobra['A'][i]
                        if not cobra['cheapObj'][obji]:
                            surrogatePlog = surrogateModels[kernel]['Objectives']['PLOGStandardized'][obji]
                            plogSol = interpRBF(x,surrogatePlog)
                            plogSol = plogSol*cobra['FresPlogStandardizedStd'][obji] + cobra['FresPlogStandardizedMean'][obji]
                            plogSol = plogReverse(plogSol)
                            surrogateErrors['OBJ'+str(obji)+'PLOG'+kernel].append((plogSol - yTrue[obji])**2)

                            surrogate = surrogateModels[kernel]['Objectives']['Standardized'][obji]
                            sol = interpRBF(x, surrogate)
                            sol = sol*cobra['FresStandardizedStd'][obji] + cobra['FresStandardizedMean'][obji]
                            surrogateErrors['OBJ'+str(obji)+kernel].append((sol - yTrue[obji])**2)
                        else:
                            surrogateErrors['OBJ'+str(obji)+'PLOG'+kernel].append(0)
                            surrogateErrors['OBJ'+str(obji)+kernel].append(0)
                            
                for coni in range(cobra['nConstraints']):
                    for kernel in cobra['RBFmodel']:
                        if not cobra['cheapCon'][coni]:
                            surrogatePlog = surrogateModels[kernel]['Constraints']['PLOGRescaled'][coni]
                            plogSol = interpRBF(x, surrogatePlog)
                            plogSol = plogReverse(plogSol)
                            plogSol = cobra['GresPlogRescaledDivider'][coni] * plogSol
                            surrogateErrors['CON'+str(coni)+'PLOG'+kernel].append((plogSol - conTrue[coni])**2)
                            
                            surrogate = surrogateModels[kernel]['Constraints']['Rescaled'][coni]
                            sol = interpRBF(x, surrogate)
                            sol = cobra['GresRescaledDivider'][coni] * sol
                            surrogateErrors['CON'+str(coni)+kernel].append((sol - conTrue[coni])**2)
                        else:
                            surrogateErrors['CON'+str(coni)+'PLOG'+kernel].append(0)
                            surrogateErrors['CON'+str(coni)+kernel].append(0)
        
        bestPredictor = {'objKernel':[], 'objLogStr':[], 'conKernel':[], 'conLogStr':[]}
        for obji in range(cobra['nObj']):
            minScore = np.finfo(np.float64).max
            bestPredictor['objKernel'].append(cobra['RBFmodel'][-1])
            bestPredictor['objLogStr'].append('Standardized')
            if not cobra['cheapObj'][obji]:
                for kernel in cobra['RBFmodel']:
                    if np.sqrt(np.mean(np.array(surrogateErrors['OBJ'+str(obji)+kernel]))) < minScore:
                        minScore = np.sqrt(np.mean(np.array(surrogateErrors['OBJ'+str(obji)+kernel])))
                        bestPredictor['objKernel'][obji] = kernel
                        bestPredictor['objLogStr'][obji] = 'Standardized'
                    if np.sqrt(np.mean(np.array(surrogateErrors['OBJ'+str(obji)+'PLOG'+kernel]))) < minScore:
                        minScore = np.sqrt(np.mean(np.array(surrogateErrors['OBJ'+str(obji)+'PLOG'+kernel])))
                        bestPredictor['objKernel'][obji] = kernel
                        bestPredictor['objLogStr'][obji] = 'PLOGStandardized'
                        
        for coni in range(cobra['nConstraints']):
            minScore = np.finfo(np.float64).max
            bestPredictor['conKernel'].append(cobra['RBFmodel'][-1])
            bestPredictor['conLogStr'].append('Rescaled')
            if not cobra['cheapCon'][coni]:
                for kernel in cobra['RBFmodel']:
                    if np.sqrt(np.mean(np.array(surrogateErrors['CON'+str(coni)+kernel]))) < minScore:
                        minScore = np.sqrt(np.mean(np.array(surrogateErrors['CON'+str(coni)+kernel])))
                        bestPredictor['conKernel'][coni] = kernel
                        bestPredictor['conLogStr'][coni] = 'Rescaled'
                    if np.sqrt(np.mean(np.array(surrogateErrors['CON'+str(coni)+'PLOG'+kernel]))) < minScore:
                        minScore = np.sqrt(np.mean(np.array(surrogateErrors['CON'+str(coni)+'PLOG'+kernel])))
                        bestPredictor['conKernel'][coni] = kernel
                        bestPredictor['conLogStr'][coni] = 'PLOGRescaled'    
        return bestPredictor
    
    def define_best_predictor(x, yTrue, conTrue, surrogateModels, cobra):
        if x is not None: # check if dict is empty in first iteration.
            for obji in range(cobra['nObj']):
                for kernel in cobra['RBFmodel']:
                    if not cobra['cheapObj'][obji]:
                        surrogatePlog = surrogateModels[kernel]['Objectives']['PLOGStandardized'][obji]
                        plogSol = interpRBF(x,surrogatePlog)
                        plogSol = plogSol*cobra['FresPlogStandardizedStd'][obji] + cobra['FresPlogStandardizedMean'][obji]
                        plogSol = plogReverse(plogSol)
                        cobra['SurrogateErrors']['OBJ'+str(obji)+'PLOG'+kernel].append((plogSol - yTrue[obji])**2)
                        
                        surrogate = surrogateModels[kernel]['Objectives']['Standardized'][obji]
                        sol = interpRBF(x, surrogate)
                        sol = sol*cobra['FresStandardizedStd'][obji] + cobra['FresStandardizedMean'][obji]
                        cobra['SurrogateErrors']['OBJ'+str(obji)+kernel].append((sol - yTrue[obji])**2)
                    else:
                        cobra['SurrogateErrors']['OBJ'+str(obji)+'PLOG'+kernel].append(0)
                        cobra['SurrogateErrors']['OBJ'+str(obji)+kernel].append(0)
            
            for coni in range(cobra['nConstraints']):
                for kernel in cobra['RBFmodel']:
                    if not cobra['cheapCon'][coni]:
                        surrogatePlog = surrogateModels[kernel]['Constraints']['PLOGRescaled'][coni]
                        plogSol = interpRBF(x, surrogatePlog)
                        plogSol = plogReverse(plogSol)                    
                        plogSol = cobra['GresPlogRescaledDivider'][coni] * plogSol
                        cobra['SurrogateErrors']['CON'+str(coni)+'PLOG'+kernel].append((plogSol - conTrue[coni])**2)
                        
                        surrogate = surrogateModels[kernel]['Constraints']['Rescaled'][coni]
                        sol = interpRBF(x, surrogate)
                        sol = cobra['GresRescaledDivider'][coni] * sol
                        cobra['SurrogateErrors']['CON'+str(coni)+kernel].append((sol - conTrue[coni])**2)
                    else:
                        cobra['SurrogateErrors']['CON'+str(coni)+'PLOG'+kernel].append(0)
                        cobra['SurrogateErrors']['CON'+str(coni)+kernel].append(0)
        
        tempErrors = copy.deepcopy(cobra['SurrogateErrors'])
                
        hvgrowing = np.zeros(len(cobra['hypervolumeProgress']))==1
        for i in range(1,len(cobra['hypervolumeProgress'])):    
            if cobra['hypervolumeProgress'][i] - cobra['hypervolumeProgress'][i-1] > 0:
                hvgrowing[i] = True
        
        hvgrowing[-2*cobra['batch']:] = True # besides the hypervolume improvement iterations, the results from the last two iterations are also taken into account
        
        bestPredictor = {'objKernel':[], 'objLogStr':[], 'conKernel':[], 'conLogStr':[]}
        for obji in range(cobra['nObj']):
            minScore = np.finfo(np.float64).max
            bestPredictor['objKernel'].append(cobra['RBFmodel'][-1])
            bestPredictor['objLogStr'].append('Standardized')
            if not cobra['cheapObj'][obji]:
                for kernel in cobra['RBFmodel']:
                    if np.sqrt(np.mean(np.array(tempErrors['OBJ'+str(obji)+kernel])[hvgrowing])) < minScore:
                        minScore = np.sqrt(np.mean(np.array(tempErrors['OBJ'+str(obji)+kernel])[hvgrowing]))
                        bestPredictor['objKernel'][obji] = kernel
                        bestPredictor['objLogStr'][obji] = 'Standardized'
                    if np.sqrt(np.mean(np.array(tempErrors['OBJ'+str(obji)+'PLOG'+kernel])[hvgrowing])) < minScore:
                        minScore = np.sqrt(np.mean(np.array(tempErrors['OBJ'+str(obji)+'PLOG'+kernel])[hvgrowing]))
                        bestPredictor['objKernel'][obji] = kernel
                        bestPredictor['objLogStr'][obji] = 'PLOGStandardized'
                        
        for coni in range(cobra['nConstraints']):
            minScore = np.finfo(np.float64).max
            bestPredictor['conKernel'].append(cobra['RBFmodel'][-1])
            bestPredictor['conLogStr'].append('Rescaled')
            if not cobra['cheapCon'][coni]:
                for kernel in cobra['RBFmodel']:
                    if np.sqrt(np.mean(np.array(tempErrors['CON'+str(coni)+kernel])[hvgrowing])) < minScore:
                        minScore = np.sqrt(np.mean(np.array(tempErrors['CON'+str(coni)+kernel])[hvgrowing]))
                        bestPredictor['conKernel'][coni] = kernel
                        bestPredictor['conLogStr'][coni] = 'Rescaled'
                    if np.sqrt(np.mean(np.array(tempErrors['CON'+str(coni)+'PLOG'+kernel])[hvgrowing])) < minScore:
                        minScore = np.sqrt(np.mean(np.array(tempErrors['CON'+str(coni)+'PLOG'+kernel])[hvgrowing]))
                        bestPredictor['conKernel'][coni] = kernel
                        bestPredictor['conLogStr'][coni] = 'PLOGRescaled'
        cobra['bestPredictor'].append(bestPredictor)     
        
    def updateInfoAndCounters(cobra, xNew, yNewEval, conNewEval, phase, surrogateModels):
        cobra['A'] = np.vstack((cobra['A'], xNew))
        cobra['lastX'] = xNew
        cobra['Fres'] = np.vstack((cobra['Fres'], yNewEval))
        cobra['Gres'] = np.vstack((cobra['Gres'], conNewEval))
        
        FresStandardized = np.full_like(cobra['Fres'], 0)
        FresStandardizedMean = np.zeros(cobra['nObj'])
        FresStandardizedStd = np.zeros(cobra['nObj'])

        FresPlogStandardized = np.full_like(cobra['Fres'], 0)
        FresPlogStandardizedMean = np.zeros(cobra['nObj'])
        FresPlogStandardizedStd = np.zeros(cobra['nObj'])
        
        for obji in range(cobra['nObj']):
            res, mean, std = standardize_obj(cobra['Fres'][:,obji])        
            FresStandardized[:,obji] = res
            FresStandardizedMean[obji] = mean
            FresStandardizedStd[obji] = std
            
            plogFres = plog(cobra['Fres'][:,obji])
            res, mean, std = standardize_obj(plogFres)
            FresPlogStandardized[:,obji] = res
            FresPlogStandardizedMean[obji] = mean 
            FresPlogStandardizedStd[obji] = std
            
        cobra['FresStandardized'] = FresStandardized
        cobra['FresStandardizedMean'] = FresStandardizedMean
        cobra['FresStandardizedStd'] = FresStandardizedStd  
        
        cobra['FresPlogStandardized'] = FresPlogStandardized
        cobra['FresPlogStandardizedMean'] = FresPlogStandardizedMean
        cobra['FresPlogStandardizedStd'] = FresPlogStandardizedStd        

        GresRescaled = np.full_like(cobra['Gres'], 0)
        GresRescaledDivider = np.zeros(cobra['nConstraints'])
        GresPlogRescaled = np.full_like(cobra['Gres'], 0)
        GresPlogRescaledDivider = np.zeros(cobra['nConstraints'])
        for coni in range(cobra['nConstraints'] ):
            GresRescaled[:,coni], GresRescaledDivider[coni] = rescale_constr(cobra['Gres'][:,coni])
            plogGres, GresPlogRescaledDivider[coni] = rescale_constr(cobra['Gres'][:,coni])
            GresPlogRescaled[:,coni] = plog(plogGres)
            
        cobra['GresRescaled'] = GresRescaled
        cobra['GresRescaledDivider'] = GresRescaledDivider
        
        cobra['GresPlogRescaled'] = GresPlogRescaled
        cobra['GresPlogRescaledDivider'] = GresPlogRescaledDivider
        
        
        pff = paretofrontFeasible(cobra['Fres'], cobra['Gres'])
        pf = cobra['Fres'][pff]
        cobra['paretoFrontier'] = pf
        cobra['paretoFrontierFeasible'] = pff
        
        hv = hypervolume(pf, cobra['ref'])
        cobra['currentHV'] = hv
        
        newNumViol = np.sum(conNewEval > 0)
        newMaxViol = max(0, max(conNewEval))
        
        if newNumViol == 0:
            cobra['hypervolumeProgress'] = np.append(cobra['hypervolumeProgress'], hv)
        else:
            cobra['hypervolumeProgress'] = np.append(cobra['hypervolumeProgress'], cobra['hypervolumeProgress'][-1])
        
        cobra['numViol'] = np.append(cobra['numViol'], newNumViol)
        cobra['maxViol'] = np.append(cobra['maxViol'], newMaxViol)
        cobra['phase'].append(phase)
        
        for ci in range(cobra['nConstraints']):
            if conNewEval[ci] <= 0:
                cobra['EPS'][ci] = np.maximum(cobra['EPS'][ci]*(1-cobra['epsilonLearningRate']), cobra['seqTol']**0.5)
            else:
                cobra['EPS'][ci] = np.minimum((cobra['epsilonLearningRate']+1)*cobra['EPS'][ci],cobra['epsilonMax'][ci])
            
        define_best_predictor(xNew, yNewEval, conNewEval, surrogateModels, cobra)
        return(cobra)

    def trainSurrogates(cobra, pool):
        surrogateModels = {}
        A = cobra['A']
        FresStandardized = cobra['FresStandardized'].T
        FresPlogStandardized = cobra['FresPlogStandardized'].T
        GresRescaled = cobra['GresRescaled'].T
        GresPlogRescaled = cobra['GresPlogRescaled'].T
        kernels = cobra['RBFmodel']
                 
        if pool:
            partialTrainSurrogate = partial(trainSurrogate, A=A, GresRescaled=GresRescaled, GresPlogRescaled=GresPlogRescaled, FresStandardized=FresStandardized, FresPlogStandardized=FresPlogStandardized, cheapConstr=cobra['cheapCon'], cheapObj=cobra['cheapObj'])
            surrogateModelsList = pool.map(partialTrainSurrogate, kernels)
            for model in surrogateModelsList:
                surrogateModels[model['type']] = model
        else:
            for kernel in kernels:
                surrogateModels[kernel] = trainSurrogate(kernel, A, GresRescaled, GresPlogRescaled, FresStandardized, FresPlogStandardized, cobra['cheapCon'], cobra['cheapObj'])
          
        return surrogateModels

    def computeStartPoints(cobra, points):
        np.random.seed(cobra['cobraSeed']+len(cobra['A'])+1)
        lb = cobra['lower']
        ub = cobra['upper']
        strategy = cobra['computeStartPointsStrategy']
        
        if strategy=='random':
            startPoints = lb + np.random.rand(len(lb)) * (ub-lb)
        
        elif strategy=='multirandom':
            startPoints = np.random.rand(len(lb)*points)
            startPoints = startPoints.reshape((points,len(lb)))
            startPoints = lb + startPoints * (ub-lb)
        
        elif strategy=='LHS':
            startPoints = lhs(len(lb), samples=points, criterion="center", iterations=5)
            
        elif strategy=='midle':
            startPoints = (lb + ub)/ 2
        else:
            # do something smart?
            raise ValueError("This strategy does not exist for computeStartPoints")
            
        return startPoints
    
    def combineBest(submins, cobra, surrogateModels, bestPredictor, fnCheap):
        submins_xs = []
        for submin in submins:
            submins_xs.extend(submin['x'])
        submins_xs = np.array(submins_xs)
        nsol = int(len(submins_xs)/cobra['nVar'])
        submins_xs = np.maximum(submins_xs, list(cobra['lower'])*nsol)
        submins_xs = np.minimum(submins_xs, list(cobra['upper'])*nsol)
        submins_xs = np.reshape(submins_xs, (nsol, cobra['nVar']))
        submins_xs = np.unique(submins_xs,axis=0)
        nsol = len(submins_xs)

        
        solutions = np.zeros((nsol,cobra['nObj']))
        constraint_solutions = np.zeros((nsol, cobra['nConstraints']))
        i = 0
        for xi in submins_xs:
            potentialSolution = get_potentialSolution(xi, surrogateModels, bestPredictor, cobra['nObj'], cobra['infillCriteria'], cobra['FresPlogStandardizedStd'], cobra['FresStandardizedStd'], cobra['FresPlogStandardizedMean'], cobra['FresStandardizedMean'], fnCheap, cobra['cheapObj'])
            solutions[i] = potentialSolution
            constraintPrediction = getConstraintPrediction(xi, surrogateModels, bestPredictor, cobra['nConstraints'], cobra['GresPlogRescaledDivider'], cobra['GresRescaledDivider'], fnCheap, cobra['cheapCon'], cobra['EPS'])
            constraint_solutions[i] = constraintPrediction
            i+=1
        
        pfandNew = np.vstack((cobra['paretoFrontier'], solutions))
        pfandNew_constraints = np.vstack((cobra['Gres'][cobra['paretoFrontierFeasible']], constraint_solutions))
        pff = paretofrontFeasible(pfandNew, pfandNew_constraints)
        pff[:len(cobra['Fres'][cobra['paretoFrontierFeasible']])] = False
        paretoDominantSolutions = pfandNew[pff]
        
        if len(paretoDominantSolutions) == cobra['batch']:
            pff = pff[-len(submins_xs):]
            return submins_xs[pff], 'exactBatchSolutions'
        
        elif len(paretoDominantSolutions) > cobra['batch']:
            bestComp = np.array(paretoDominantSolutions[:cobra['batch']])
            newSet = np.vstack((cobra['paretoFrontier'], np.array(bestComp)))
            bestHyp = hypervolume(newSet, cobra['ref'])
            for i in range(10000):
                comb = paretoDominantSolutions[np.random.choice(len(paretoDominantSolutions), size=cobra['batch'], replace=False)]
                newSet = np.vstack((cobra['paretoFrontier'], np.array(comb)))
                hvres = hypervolume(newSet, cobra['ref'])
                if hvres > bestHyp:
                    bestHyp = hvres
                    bestComp = np.array(comb)
            
            indicator = [False]*len(paretoDominantSolutions)
            indicator_i = 0
            for solution in paretoDominantSolutions:
                if any(np.sum(np.equal(solution, bestComp), axis=1)==len(bestComp[0])):
                    indicator[indicator_i] = True
                indicator_i += 1            
            pff = pff[-len(submins_xs):]
            return submins_xs[pff][indicator], 'bestSolutionsFromBatch'
        else:
            constraintViolations = []
            for xi in submins_xs:
                potC = getEquallyImportantConstraintPrediction(xi, surrogateModels, bestPredictor, cobra['nConstraints'], cobra['GresRescaledDivider'], fnCheap, cobra['cheapCon'], EPS=None)
                potC = np.array(potC)
                constraintViolation = potC[potC>0]
                constraintViolations.append(sum(constraintViolation))
            max_c = max(constraintViolations)
            pff = pff[-len(submins_xs):]
            while sum(pff)<cobra['batch']:
                nextMinViolation = np.argmin(constraintViolations)
                pff[nextMinViolation] = True
                constraintViolations[nextMinViolation] = max_c+1            
            x_results = submins_xs[pff]
            return x_results, 'smallestConstraintViolation'
            
    def findSurrogateMinimum(cobra, surrogateModels, bestPredictor, pool):
        submins = []
        # besti = 0
        # bestFun = 0
        # success = []
        cons = []
        
        if cobra['batch'] == 1:
            xStarts = computeStartPoints(cobra, cobra['computeStartingPoints'])
            gCOBRA_partial = partial(gCOBRA, A=cobra['A'], lower=cobra['lower'], upper=cobra['upper'], surrogateModels=surrogateModels, bestPredictor=bestPredictor, nConstraints=cobra['nConstraints'], GresRescaledDivider=cobra['GresRescaledDivider'], fnCheap=fnCheap, cheapCon=cobra['cheapCon'], EPS=cobra['EPS'])
            cons.append({'type':'ineq','fun':gCOBRA_partial})
            compute_infill_criteria_score_partial = partial(compute_infill_criteria_score, surrogateModels=surrogateModels, bestPredictor=bestPredictor, nObj=cobra['nObj'], infillCriteria=cobra['infillCriteria'], FresPlogStandardizedStd=cobra['FresPlogStandardizedStd'], FresStandardizedStd=cobra['FresStandardizedStd'], FresPlogStandardizedMean=cobra['FresPlogStandardizedMean'], FresStandardizedMean=cobra['FresStandardizedMean'], currentHV=cobra['currentHV'], paretoFrontier=cobra['paretoFrontier'], ref=cobra['ref'], fnCheap=fnCheap, cheapObj=cobra['cheapObj'])
        else:
            xStarts = computeStartPoints(cobra, cobra['computeStartingPoints']*cobra['batch'])
            xStarts = xStarts.reshape((cobra['computeStartingPoints'],len(cobra['lower'])*cobra['batch']))
            gCOBRA_partial = partial(batch_gCOBRA, batch=cobra['batch'], A=cobra['A'], lower=cobra['lower'], upper=cobra['upper'], surrogateModels=surrogateModels, bestPredictor=bestPredictor, nConstraints=cobra['nConstraints'], GresRescaledDivider=cobra['GresRescaledDivider'], fnCheap=fnCheap, cheapCon=cobra['cheapCon'], EPS=cobra['EPS'])
            cons.append({'type':'ineq','fun':gCOBRA_partial})
            compute_infill_criteria_score_partial = partial(batch_infill_criteria_score, batch=cobra['batch'], surrogateModels=surrogateModels, bestPredictor=bestPredictor, nObj=cobra['nObj'], infillCriteria=cobra['infillCriteria'], FresPlogStandardizedStd=cobra['FresPlogStandardizedStd'], FresStandardizedStd=cobra['FresStandardizedStd'], FresPlogStandardizedMean=cobra['FresPlogStandardizedMean'], FresStandardizedMean=cobra['FresStandardizedMean'], paretoFrontier=cobra['paretoFrontier'], ref=cobra['ref'], fnCheap=fnCheap, cheapObj=cobra['cheapObj'])
        f = partial(pool_job, criteria_function=compute_infill_criteria_score_partial, cons=cons, seqFeval = cobra['seqFeval'], seqTol=cobra['seqTol'])
        if pool:
            submins = pool.map(f, xStarts)
        else:
            for xStart in xStarts:
                submins.append(f(xStart))

        # for i in range(len(submins)):
        #     subMin = submins[i]
        #     success.append(subMin['success'])
        #     if subMin['fun'] < bestFun and subMin['success']:
        #         bestFun = subMin['fun']
        #         besti = i
        
        xNew, status = combineBest(submins, cobra, surrogateModels, bestPredictor, fnCheap)
        
        # if all(success):
        #     minRequiredEvaluations = (cobra['dimension']+cobra['nConstraints']+cobra['nObj']+cobra['batch'])*20
        #     adjustedAmountEvaluations = int(cobra['seqFeval']*(1-cobra['surrogateUpdateLearningRate']))
        #     cobra['seqFeval'] = max(adjustedAmountEvaluations, minRequiredEvaluations)
            
        #     maxStartingPoints = (cobra['dimension']+cobra['nConstraints']+cobra['nObj'])*10*(1+int(cobra['batch']>1))*(1+int(cobra['oneShot']))
        #     adjustedAmountPoints = int(cobra['computeStartingPoints']*(1+cobra['surrogateUpdateLearningRate']))
        #     cobra['computeStartingPoints'] = min(maxStartingPoints, adjustedAmountPoints)
        # else:
        #     maxRequiredEvaluations = (cobra['dimension']+cobra['nConstraints']+cobra['nObj'])*1000*(1+int(cobra['batch']>1))*(1+int(cobra['oneShot']))
        #     adjustedAmountEvaluations = int(cobra['seqFeval']*(1+cobra['surrogateUpdateLearningRate']))
        #     cobra['seqFeval'] = min(adjustedAmountEvaluations, maxRequiredEvaluations)
            
        #     minRequiredPoints = 2*(cobra['dimension']+cobra['nConstraints']+cobra['nObj']+cobra['batch'])
        #     adjustedAmountPoints = int(cobra['computeStartingPoints']*(1-cobra['surrogateUpdateLearningRate']))
        #     cobra['computeStartingPoints'] = max(adjustedAmountPoints, minRequiredPoints)
                
        cobra['optimizerConvergence'] = np.append(cobra['optimizerConvergence'], status)        
        return xNew   
    
    def get_best_strategy(cobra, pool):
        
        evaluations_stuck_max = 3*cobra['batch']
        a = cobra['hypervolumeProgress'][-1*evaluations_stuck_max:][0]
        b = cobra['hypervolumeProgress'][-1]
        cobra['infillCriteria'] = 'PHV'
        if b-a < cobra['seqTol'] and len(cobra['hypervolumeProgress'])>3*evaluations_stuck_max:
            cobra['infillCriteria'] = 'SMS'
            print('criteria changed to SMS')

        if cobra['oneShot'] or len(cobra['Fres'])>30:
            bestPredictor = kfold_best_predictor(cobra, pool)
            return bestPredictor
        else:
            return cobra['bestPredictor'][-1]
    
    xNew = [None]
    yNewEval = [None]
    conNewEval = [None]
    surrogateModels = {} 

    while n < cobra['feval']:
        ptm = time.time()
        
        bestPredictor = get_best_strategy(cobra, pool)
        
        surrogateModels = trainSurrogates(cobra, pool)        
        
        xNew = findSurrogateMinimum(cobra, surrogateModels, bestPredictor, pool)
        
        yNewEval = []
        conNewEval = []
        if cobra['batch'] == 1:
            xNew = xNew[0]
            stttt = time.time()
            yNew, cNew = fn(xNew)
            cobra['evaluationTime'].append(time.time()-stttt)
            yNewEval.append(yNew)
            conNewEval.append(cNew)
            cobra = updateInfoAndCounters(cobra, xNew, yNew, cNew, phase, surrogateModels)
            xNew = np.array([xNew])
        else:
            # xNew = xNew.reshape(cobra['batch'],int(len(xNew)/cobra['batch']))
            if pool:
                stttt = time.time()
                res = pool.map(fn, xNew)
                cobra['evaluationTime'].append(time.time()-stttt)
            else:
                res = []
                for row in xNew:
                    stttt = time.time()
                    res.append(fn(row))
                    cobra['evaluationTime'].append(time.time()-stttt)
            i = 0
            for result in res:
                yNew, cNew = result
                yNewEval.append(yNew)
                conNewEval.append(cNew)
                cobra = updateInfoAndCounters(cobra, xNew[i], yNew, cNew, phase, surrogateModels)
                i += 1
                
        n = len(cobra['A'])
        
        cobra['optimizationTime'] = np.append(cobra['optimizationTime'], time.time()-ptm)

        if cobra['plot']:
            print(time.time()-ptm, n, cobra['feval'], cobra['hypervolumeProgress'][-1])
            print(cobra['paretoFrontier'])
            visualiseParetoFront(cobra['paretoFrontier'])
    
    functionName = str(fn).split(' ')[6].split('.')[0]
    if cobra['oneShot']:
        outdir = 'oneshot/'+str(functionName)+'/'    
    else:
        outdir = 'batchresults/'+str(functionName)+'/batch'+str(cobra['batch'])+'/'
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    paretoOptimal = paretofrontFeasible(cobra['Fres'],cobra['Gres'])
    paretoFront = cobra['Fres'][paretoOptimal]
    paretoSet = cobra['A'][paretoOptimal]
    paretoConstraints = cobra['Gres'][paretoOptimal]
    runNo = cobra['cobraSeed']

    outputFileParameters = str(outdir)+'par_run'+str(runNo)+'_final.csv'
    outputFileObjectives = str(outdir)+'obj_run'+str(runNo)+'_final.csv'
    outputFileConstraints = str(outdir)+'con_run'+str(runNo)+'_final.csv'
    
    np.savetxt(outputFileParameters, cobra['A'], delimiter=',')
    np.savetxt(outputFileObjectives, cobra['Fres'], delimiter=',')
    np.savetxt(outputFileConstraints, cobra['Gres'], delimiter=',')

    outputFileParameters = str(outdir)+'par_run'+str(runNo)+'_finalPF.csv'
    outputFileObjectives = str(outdir)+'obj_run'+str(runNo)+'_finalPF.csv'
    outputFileConstraints = str(outdir)+'con_run'+str(runNo)+'_finalPF.csv'
    
    np.savetxt(outputFileObjectives, paretoFront, delimiter=',')
    np.savetxt(outputFileParameters, paretoSet, delimiter=',')
    np.savetxt(outputFileConstraints, paretoConstraints, delimiter=',')
    
    if pool:
        pool.close()
        pool.join()

    jsonRes = make_json_serializable(cobra)
    # print(jsonRes.keys())
    with open(str(outdir)+'cobra'+str(runNo)+'.json', 'w') as outfile:
        json.dump(jsonRes, outfile, indent = 3)
    
    return(cobra)
