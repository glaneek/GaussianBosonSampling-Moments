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

#Squeezing arguments
SQargs=[0.4,0,0.5,0.4]

PHargs=np.array(PHargs)
BSargs=np.array(BSargs)
#-----------------------------------------------------------#
rd.seed(10)
err=[]
for i in range(16):
    r=rd.random()
    err.append(0.1*r-0.05)
#%%
err[1]=0
err[2]=0
err[3]=0
err[4]=0
err[5]=0
err[6]=0
err[7]=0
err[8]=0
err[9]=0
err[10]=0
err[11]=0
err[12]=0
err[13]=0
err[14]=0
err[15]=0
bo=[[Pi/4,Pi/4+0.01],[Pi/4-0.01,Pi/4+0],
    [Pi/4,Pi/4+0.01],[Pi/4-0.03,Pi/4+0],
    [Pi/4,Pi/4+0.032],[Pi/4-0,Pi/4+0.032],
    [Pi/4-0,Pi/4+0.02],[Pi/4-0.035,Pi/4+0]]

bo=[[Pi/4+0.007, Pi/4+0.0072]]

#Load Data
filepath="Exp_DATA.txt"
states=[]
states=extract_data(filepath)

x=[np.pi/4,np.pi/4,np.pi/4,np.pi/4,
   np.pi/4,np.pi/4,np.pi/4,np.pi/4,
   np.pi/4,0,np.pi/4,0,
   0.4,0,0.5,0.4]

x=[np.pi/4,np.pi/4,np.pi/4,np.pi/4,
   np.pi/4,np.pi/4,np.pi/4,np.pi/4]

x=[np.pi/4+0.007]

argss=[4,1,3,2,1,3,2,1,3,states]

rss=RSS(x)
print(rss)
print("Initiate optimization")
from scipy.optimize import minimize
x0 = np.array(x)
res = minimize(RSS,x0,tol=1e-10,bounds=bo)
errE=res.x

err1=np.array(x)+np.array(err[0])
print("Optimasition done, ready to plot")


plt.figure()
plt.plot(err1,'x-',label="Original")
plt.plot(errE,'x-',label="Optimization Error")
plt.plot(x,'x-',label="Perfect")
plt.legend()
plt.show()
print(RSS(res.x))

#%%
inpu=np.arange(0.78,0.8,0.001)

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

plt.vlines(x,bmi,ama,label="%.4f"%x[0],colors='g')


plt.legend()
plt.show()

print(mi)
