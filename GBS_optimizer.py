#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 00:58:50 2020

@author: frk
"""
import numpy as np
import matplotlib.pyplot as plt 
from GBS_functions import *
import time
Pi=np.pi

#Moments to calculate
allmodes2=[[1,1],[1,2],[1,3],[1,4],[2,2],[2,3],[2,4],[3,3],[3,4],[4,4]]
allmodes1=[1,2,3,4]

def Vfinal_symplectic_approx(gaussrand):
    #How many modes?
    dim=4
    #Beam splitter sequence and arguments
    BeamS=np.array([1,3,2,1,3,2,1,3])
    BeamS=BeamS-1
    #For each beam splitter provide appropriate angle [(pi/2,0),(pi/3,0)...] etc.
    BSargs = [(Pi/4, 0),(Pi/4,0),(Pi/4, 0),
              ((Pi/4) +gaussrand, 0),(Pi/4,0),(Pi/4, 0),
              (Pi/4, 0),(Pi/4,0)]
    
    #Phase shifters' arguments
    PHargs=[0,0,0,0]
    
    #Squeezing arguments
    SQargs=[0.4,0,0.4,0.4]
        
    #Fetch symplectic covariance for comparisons
    Vfinal=Vfinal_symplectic(dim,SQargs,PHargs,BSargs,BeamS)
    
    return Vfinal

#Prepare data and sample to simulate experimental data!
filepath="GBS_BS_gauss.txt"

#Extract Data
states=[]
states=extract_data(filepath)

#Get CDF of data
CDF=get_CDF(states)

#Sample PDF and get sampled CDF
#np.random.seed(250)    
X = run_sampling(CDF,Nruns=10e4)

#Bring sampled PDF to correct format
ExpStates=np.array(np.copy(states))    
for i in range(len(states)):
    (ExpStates[i])[1]=X[i]

##########--DATA READY--##########

##########--OPTIMISATION--#########
gaussrand=0
threshold=10e-1
diff=1
guesses=[0]
COMP=[gaussrand]
progCOMP=[]
for q in range(80):
    #Get Vfinal matrix from symplectic
    Vfinal=Vfinal_symplectic_approx(gaussrand)
    #Get sampled moments and actual moments
    Vexp, actual, diag=get_sampled_moments(allmodes1,allmodes2,states,Vfinal)
    #Compare
    if q==0:
        comp1=np.sum(comparison(Vexp,actual))
        gaussrand=gaussrand+0.005
        print("HERE\n")
        continue
    else:
        print("Roound 1+")
        comp=np.copy(comp1)
        comp1=np.sum(comparison(Vexp,actual))
        diff=comp1-comp
        
        
        # if (abs(diff)-threshold)<0 and q>1:
        #     print("SHOULDNT BE FUCKING HERE")
        #     print(q)
        #     print(gaussrand)
        #     print(diff)
        #     print(diff+threshold)
        #     print(diff-threshold)
        #     break
        
        if diff<0:
            gaussrand=gaussrand+0.005
        elif diff>0:
            gaussrand=gaussrand-0.005
            
    progCOMP.append(comp1)
    COMP.append(comp1)
    print(q)
    print(gaussrand)
    print(diff)
    guesses.append(gaussrand)
#%%
guesses=guesses[:-1]
x=np.ones(len(guesses))*0.1
plt.figure()
plt.plot(x,"-b",label="True Value",linewidth=4)
plt.plot(guesses,'r',label="Numerical Approximation",linewidth=3)
plt.plot(progCOMP)
plt.ylabel("BS deviation")
plt.xlabel("Iterations")
plt.grid()
plt.legend()
plt.show()

#%%
#print(gaussrand)
# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.plot(x,"-b",label="True Value",linewidth=4)
ax.plot(guesses,'r',label="Numerical Approximation",linewidth=3)
ax.set_ylabel("BS deviation ")
ax.set_xlabel("Iterations")
ax.grid()

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(progCOMP,'-g',label="SSR",linewidth=4)
ax2.set_ylabel("Square Sum of Residuals")
ax.legend()
ax2.legend()
#%%
time = np.arange(10)
temp = np.random.random(10)*30
Swdown = np.random.random(10)*100-10
Rn = np.random.random(10)*100-10

fig = plt.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot( x, '-b', label = "True Error Value",linewidth=3)
lns2 = ax.plot( guesses, '-r', label = "Numerical Approximation",linewidth=3)
ax2 = ax.twinx()
lns3 = ax2.plot(progCOMP, '-g', label = "Squared Sum of Residuals",linewidth=3)

# added these three lines
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=7)

ax.grid()
ax.set_xlabel("Iterations")
ax.set_ylabel("BS deviation ")
ax2.set_ylabel("Deviation from theory")

plt.show()

        
    

