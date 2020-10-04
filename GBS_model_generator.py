#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 20:56:19 2020

@author: frk
"""""

import numpy as np
from numpy.linalg import multi_dot
import strawberryfields as sf
from strawberryfields.ops import *
import math
from GBS_functions import beam_splitter_args
import random as rd
import datetime
Pi=np.pi;

t=datetime.datetime.now().time()
#filepath="GBS_BS_gauss.txt"
filepath="DATA.txt"
rd.seed(10)
err=[]
for i in range(16):
    r=rd.random()
    err.append(1.1*r-1)
    
#######Define your system#######
#How many modes?
dim=4
#Beam splitter sequence and arguments
BeamS=np.array([1,3,2,1,3,2,1,3])
BeamS=BeamS-1
#For each beam splitter provide appropriate angle [(pi/2,0),(pi/3,0)...] etc.
BSargs = [(Pi/4+err[0], 0),
          (Pi/4+err[1], 0),
          (Pi/4+err[2], 0),
          (Pi/4+err[3], 0),
          (Pi/4+err[4], 0),
          (Pi/4+err[5], 0),
          (Pi/4+err[6], 0),
          (Pi/4+err[7],0)]

#Phase shifters' arguments
PHargs=[Pi/4+err[8],
        0+err[9],
        Pi/4+err[10],
        0+err[11]]

#Squeezing arguments
SQargs=[0.4+err[12],
        0+err[13],
        0.5+err[14],
        0.4+err[15]]

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
N=12
for h in range(N):
    for h1 in range(N):
        for h2 in range(N):
            for h3 in range(N):
                if (h+h1+h2+h3)%2==0:
                    measure_states.append([h,h1,h2,h3])


PROB=[]
p=-1
non=0
#Loop through states and get results
for i in measure_states:
    prob = float(results.state.fock_prob(i))
    p+=1
    print("|{}>: {}".format("".join(str(j) for j in i), prob))
    
    #write to file
    file.write("|{}>: {}\n".format("".join(str(j) for j in i), prob))

    
    #Check and count Nan
    if math.isnan(prob):
        non=non+1
        continue;
    
    PROB.append(prob)

#Check for nan and raise warning
if non>0:
    print("NaN values encountered\n")

#Print total probability
cum_prob=np.cumsum(PROB)
print("Cumsum probability:  %.4f\n" %cum_prob[-1])

#Close file
file.close()




