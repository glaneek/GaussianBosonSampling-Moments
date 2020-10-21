#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:21:13 2020

@author: frk


"""
import numpy as np
from numpy.linalg import multi_dot
import math
from GBS_functions import *
import random as rd
import datetime
#from GBS_model_generator import err
import matplotlib.pyplot as plt
from scipy.optimize import minimize


Pi=np.pi

#-----------------------------------------------------------#
#Guessed Parameters
#How many modes?
dim=4
#Beam splitter sequence and arguments
BeamS=np.array([1,3,2,1,3,2,1,3])
BeamS=BeamS-1
#For each beam splitter provide appropriate angle [(pi/2,0),(pi/3,0)...] etc.
BSargs = [(Pi/4, 0),(Pi/4,0),(Pi/4, 0),
          (Pi/4, 0),(Pi/4,0),(Pi/4, 0),
          (Pi/4, 0),(Pi/4,0)]

#Phase shifters' arguments
PHargs=[Pi/4,0,Pi/4,0]
PHargs=[0,0,0,0]



#Squeezing arguments
SQargs=[0.4,0,0,0.4]

PHargs=np.array(PHargs)
BSargs=np.array(BSargs)
#-----------------------------------------------------------#

#-----------------------------------------------------------------------
###############FASTER APPROACH##########################################
#-----------------------------------------------------------------------
allmodes2=[[1,1],[1,2],[1,3],[1,4],[2,2],[2,3],[2,4],[3,3],[3,4],[4,4]]
allmodes1=[1,2,3,4]

filepath="DATA_all.txt"
states=[]
states=extract_data(filepath)
##########--COMPARISON--#########
#Get sampled moments and actual moments
momends1, momends2 = get_sampled_moments(allmodes1,allmodes2,states)
#%%
TOL=[allmodes1,allmodes2,momends1,momends2,[],[]]

def all(TOL,PARINDEX):
    rd.seed(10)
    err=[]
    for i in range(16):
        r=rd.random()
        err.append(0.04*r-0.02)
        
    #True values of ALL beam splitters
    BS=np.array(err[:8])+Pi/4
    
    # #Index True (Known), False (To be found)
    # PARINDEX=[False, False,
    #           False, False,
    #           False, True,
    #           True, True]
    
    #Assign initial condition of Pi/4 to unknown parameters
    #as starting points
    x0=[]
    err1=[]
    for i in range(len(PARINDEX)):
        if PARINDEX[i]==False:
            x0.append(Pi/4)
            err1.append(BS[i])
    
    print("Initiate optimization")
    
    BSargs=[]
    for i in range(len(PARINDEX)):
        if PARINDEX[i]:
            BSargs.append((BS[i],0))
        else:
            BSargs.append((0,0))
            
    TOL[4]=BSargs
    TOL[5]=PARINDEX
    
    count=0
    for i in range(len(BSargs)):
        if (BSargs[i])[0]==0:
            BSargs[i]=(x0[count],0)
            count+=1
    
    rss=RSS(x0,TOL)
    print(rss)
    
    res = minimize(RSS,x0,args=TOL)
    errE=res.x
    
    #err1=np.array(BS)
    #err1=np.array(x)+np.hstack([err[0:6],err[7]])
    print("Optimasition done, ready to plot")
    
    print(RSS(res.x,TOL))
    print(RSS(err1,TOL))
    
    return err1, errE, x0

#%%

temp=[False, False, False, False,
      False, False, False, False]

temp=[True, True, True, True,
      True, True, True, True]

all_ind=[]

for k in range(len(temp)):
    temp[k]=False
    all_ind.append(np.copy(temp))
    temp[k]=True
    

plt.figure()

for p in all_ind:
    err1, errE, x0=all(TOL,p)
    #plt.plot(err1,'o-',label="Original")
    plt.plot(errE,'x-',label="Optimization Error")
    #plt.plot(x0, 'x-', label='Starting')
    
plt.show()

#%%

plt.figure()
plt.subplot(2,1,1)
plt.plot(err1,'o-',label="Original")
plt.plot(errE,'x-',label="Optimization Error")
plt.plot(x0, 'x-', label='Starting')
plt.legend()
plt.subplot(2,1,2)
percR=(err1-errE)/err1 *100
plt.plot(percR,'x-',label="Rediduals Percentage")
plt.legend()
plt.show()

#%%
inpu=np.arange(0.76,0.8,0.005)

aoutput=[]

filepath="Exp_DATA.txt"
for i in inpu:
    aoutput.append(RSS(i,filepath))
    print(i)
ama=max(aoutput)
ami=min(aoutput)
amin=inpu[np.argmin(aoutput)]

boutput=[]
filepath="DATA.txt"
for i in inpu:
    boutput.append(RSS(i,filepath))
    print(i)
bma=max(boutput)
bmi=min(boutput)    
bmin=inpu[np.argmin(boutput)]

plt.figure()

plt.plot(inpu,aoutput,label="Exp_DATA.txt")
plt.vlines(amin,ami,ama,label="%.4f"%amin)

plt.plot(inpu,boutput,label="DATA.txt")
plt.vlines(bmin,bmi,bma,label="%.4f"%bmin,colors='r')

plt.vlines(x,ami,bma,label="%.4f"%x[0],colors='g')


plt.legend()
plt.show()

print(mi)
