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
from GBS_model_generator import err
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

#Load Data
filepath="Exp_DATA.txt"
states=[]
states=extract_data(filepath)

x=[np.pi/4,np.pi/4,np.pi/4,np.pi/4,
   np.pi/4,np.pi/4,np.pi/4,np.pi/4,
   np.pi/4,0,np.pi/4,0,
   0.4,0,0.5,0.4]
argss=[4,1,3,2,1,3,2,1,3,states]


rss=RSS(x)

print("Initiate optimization")
from scipy.optimize import minimize
x0 = np.array(x)
res = minimize(RSS,x0,method='TNC')
errE=res.x

err1=np.array(x)+np.array(err)
plt.figure()
plt.plot(err1,label="err")
plt.plot(errE,label="errE")
plt.legend()
plt.show()