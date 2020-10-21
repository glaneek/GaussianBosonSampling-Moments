#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 20:56:19 2020

@author: frk

This script generates simulated data for a GBS
experiment. The system is defined at the beginning
where the user can define the number of modes (dim),
the beam splitter sequence and their values 
(BeamS, BSargs respectively), the initial phase of 
each mode (PHargs) and the real squeezing parameter "r"
for each mode in SQargs. 

Step 1
The script simulates the system and generate exact
probabilities with the strawberry fields package
for a set of ouput states defined by the parameter
"N": [0,0,0,0],...,[0,0,0,N],....,[N,N,N,N] (N**dim states)
The data are writen by default to a .txt file "DATA.txt". For 
good results aim for a high final cumulative probabilities
i.e 98+%

Step 2
To simulate a real experiment where only a finite sample size
is available, this scripts reconstructs another data file 
"Exp_DATA.txt" where the probabilities for each state are obtained
through random sampling using the exact probabilities in "DATA.txt".
The number of samples is defined by the variable Nruns. 
Default: Nruns=10e4
"""""

import numpy as np
from numpy.linalg import multi_dot
import strawberryfields as sf
from strawberryfields.ops import *
import math
from GBS_functions import *
import random as rd
import datetime


filepath="DATA0134567.txt"
rd.seed(10)
err=[]
for i in range(16):
    r=rd.random()
    err.append(0.04*r-0.02)

#err[0]=0
#err[1]=0
err[2]=0
# err[3]=0

# err[4]=0
#err[5]=0
#err[6]=0
#err[7]=0

err[8]=0
err[9]=0
err[10]=0
err[11]=0
err[12]=0
err[13]=0
err[14]=0
err[15]=0


def generate(err,filepath):

    print("\n--------Simulating experimental data---------\n")
    
    #############################################################
    ########################--STEP 1--###########################
    #############################################################
    
    print("1. Calculating Exact probabilties")
    
    #-----------------------------------#
    #######--DEFINE SYSTEM HERE--########
    
    Pi=np.pi;
    t=datetime.datetime.now().time()
    #filepath="GBS_BS_gauss.txt"
        
    #######Define your system#######
    #How many modes?
    dim=4
    #Beam splitter sequence and arguments
    BeamS=np.array([1,3,2,1,3,2,1,3])
    BeamS=BeamS-1
    #For each beam splitter provide appropriate 
    #angle [(pi/2,0),(pi/3,0)...] etc.
    BSargs = [(Pi/4+err[0], 0),
              (Pi/4+err[1], 0),
              (Pi/4+err[2], 0),
              (Pi/4+err[3], 0),
              (Pi/4+err[4], 0),
              (Pi/4+err[5], 0),
              (Pi/4+err[6], 0),
              (Pi/4+err[7],0)]
    
    #Phase shifters' arguments
    PHargs=[0+err[8],
            0+err[9],
            0+err[10],
            0+err[11]]
    
    #Squeezing arguments
    SQargs=[0.4+err[12],
            0+err[13],
            0+err[14],
            0.4+err[15]]
    #-----------------------------------#
    
    
    #Generate total unitary matrix for BS, PH
    UBS=beam_splitter_args(dim,BSargs,BeamS)
    Uphase = np.diag([np.exp(PHargs[0]*1j),
                      np.exp(PHargs[1]*1j),
                      np.exp(PHargs[2]*1j),
                      np.exp(PHargs[3]*1j)])
    U = multi_dot([UBS,Uphase])
    
    #Open file to write data
    file = open(filepath, "w")
    
    #Initiate program
    gbs = sf.Program(dim)
    
    with gbs.context as q:
        # prepare the input squeezed states
        
        Sgate(SQargs[0]) | q[0]
        Sgate(SQargs[1]) | q[1]
        Sgate(SQargs[2]) | q[2]
        Sgate(SQargs[3]) | q[3]
    
        # linear interferometer
        Interferometer(U) | q
    
    #Define egnine and get results
    eng = sf.Engine(backend="gaussian")
    results = eng.run(gbs)
    
    #States to loop over. Defualt N=6 -> [0,0,0,0]...[5,5,5,5]
    measure_states=[]
    N=13
    print("Generate states from [0,0,0,0] to [N,N,N,N] with N=%d"%N)
    count=0
    for h in range(N):
        for h1 in range(N):
            for h2 in range(N):
                for h3 in range(N):
                    count+=1
                    perc=count/pow(N,4)*100
                    print("\r%.2f%%"%perc,end="")
                    #time.sleep(0.01)
                    if (h+h1+h2+h3)%2==0:
                        measure_states.append([h,h1,h2,h3])
    print("\nStates generated")
    
    print("Fetching exact probabilities")
    
    
    PROB=[]
    p=-1
    non=0
    count=0
    #Loop through states and get results
    for i in measure_states:
        prob = float(results.state.fock_prob(i))
        p+=1
        #print("|{}>: {}".format("".join(str(j) for j in i), prob))
        
        #write to file
        file.write("|{}>: {}\n".format("".join(str(j) for j in i), prob))
    
        
        #Check and count Nan
        if math.isnan(prob):
            non=non+1
            continue;
        
        PROB.append(prob)
        
        count+=1
        perc=count/pow(N,4)*200
        print("\r%.2f%%"%perc,end="")
    
    
    #Check for nan and raise warning
    if non>0:
        print("NaN values encountered\n")
    
    #Print total probability
    cum_prob=np.cumsum(PROB)*100
    print("\nCumsum probability: %.2f%%" %cum_prob[-1])
    
    #Close file
    file.close()
    
    print("Exact probabilities calculated\n")
    
    
    #############################################################
    ########################--STEP 2--###########################
    #############################################################
    
    print("2. Sample from exact probability distribution")
    print("Extract exact probabilities")
    #Extract Exact Data
    #Extract Data
    states=[]
    states=extract_data(filepath)
    Nruns=10e4
    
    
    print("Optain CDF and get Nruns=%d samples"%Nruns)
    #Get CDF of data
    CDF=get_CDF(states)
    X = run_sampling(CDF,Nruns)
    print("--Sampling Done--")
    
    
    print("Bring to right format")
    #Bring sampled PDF to correct format
    ExpStates=np.array(np.copy(states))    
    for i in range(len(states)):
        (ExpStates[i])[1]=X[i]
    
    
    print("Write to data file")
    #Open file to write data
    file = open("Exp_DATA.txt", "w")
    for q in ExpStates:
        i=q[0]
        prob = float(q[1])
        p+=1
        #print("|{}>: {}".format("".join(str(j) for j in i), prob))    
        #write to file
        file.write("|{}>: {}\n".format("".join(str(j) for j in i), prob))
    #Close file
    file.close()
    print("Done: Exp_DATA.txt generated")
    
