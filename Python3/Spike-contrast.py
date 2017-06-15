# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 17:59:39 2017

@author: Manuel Ciba
"""

import numpy as np
import time

def generateRandomTestData(num_trains, num_spikes, T):
    spike_trains = np.random.random((num_trains, num_spikes))   # random values between 0 and 1
    spike_trains = np.sort(spike_trains)                        # sort values
    spike_trains = spike_trains * T
    return spike_trains

def generateTestData():                       
    spike_trains = np.array([[1, 1 ,1], [2, 3, 1.5], [5, 8, 6], [9, 10, 9]])
    T = 11
    return spike_trains, T

def get_Theta_and_n_perBin(spike_trains, T, binSize, binStepFactor):
   
    if binStepFactor == 2:
        # init
        binStep = binSize / binStepFactor
        edges = np.arange(0, T+binStep, binStep)
        tmp = spike_trains.shape                                    # number of spike trains
        N = tmp[1]
        hist = np.zeros((len(edges)-1,N));
        # calculate histogram for every spike train
        for i in range(0, N):                                       # for all spike trains
            hist[:,i] = binning_halfOverlap(spike_trains[:,i], T, binSize)
        # calculate parameter over all spike trains
        theta = np.sum(hist, 1)                                      # number of spikes per bin 
        mask = hist != 0
        n = np.sum(mask, 1)                                          # number of active spike trains
        
    return theta, n

def binning_halfOverlap(y, T, binSize, flag_binary=False):
    binStep = binSize/2
    edges = np.arange(0, T+binStep, binStep)
    hist, bin_edges = np.histogram(y, edges)
    hist[0:len(hist)-1] = hist[0:len(hist)-1] + hist[1:len(hist)] 
    return hist

"""
Main Function: SpikeContrast
"""

def SpikeContrast(spike_trains, T):
    # parameter
    binShrinkFactor = 0.9                                       # bin size decreases by factor 0.9 for each iteration
    binStepFactor = 2                                           # bin overlap of 50 %
    binMax = T/2
    isi = np.diff(spike_trains, axis=0)
    isiMin = np.min(isi)
    binMin = np.max([isiMin, 0.01])
    
    tmp = spike_trains.shape                                    # number of spike trains
    N = tmp[1]
    
    # initialization
    numIterations = np.ceil(np.log(binMin / binMax) / np.log(binShrinkFactor))
    numIterations = int(numIterations)
    ActiveST = np.zeros((numIterations, 1))
    Contrast = np.zeros((numIterations, 1))
    C = np.zeros((numIterations, 1))
    
    numAllSpikes = spike_trains.size
    binSize = binMax
    
    for i in range(0, numIterations):                           # for 0, 1, 2, ... numIterations
             
        # calculate Theta and n
        Theta_k, n_k = get_Theta_and_n_perBin(spike_trains, T, binSize, binStepFactor)
        # calcuate C= Contrast * ActiveST
        ActiveST[i] = ((np.sum(n_k*Theta_k))/(np.sum(Theta_k))-1) / (N-1)
        Contrast[i] = (np.sum(np.abs(np.diff(Theta_k)))/ (numAllSpikes*2))                    # Contrast: sum(|derivation|) / (2*#Spikes)
        C[i] = Contrast[i] * ActiveST[i]                                                      # Synchrony curve "C"
        
        binSize *= binShrinkFactor                              # new bin size
    
    # Sync value is maximum of cost function C
    S = np.max(C)
    return S


"""
Main Program
"""

spike_trains, T = generateTestData()
#spike_trains = generateRandomTestData(3, 5, 60)
S = SpikeContrast(spike_trains, T)


T = 100
spike_trains = generateRandomTestData(1000,1000,T)
start = time.time()
S = SpikeContrast(spike_trains, T)
end = time.time()
print(end - start)
