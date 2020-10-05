#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 23:54:50 2020

@author: frk
"""

import numpy as np
import matplotlib.pyplot as plt 
from GBS_functions import *
import time
Pi=np.pi

#######Define your system#######
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
PHargs=[Pi/4,0,0,0]

#Squeezing arguments
SQargs=[0.4,0,0.5,0.4]

#Moments to calculate
allmodes2=[[1,1],[1,2],[1,3],[1,4],[2,2],[2,3],[2,4],[3,3],[3,4],[4,4]]
allmodes1=[1,2,3,4]


PHargs=np.array(PHargs)
BSargs=np.array(BSargs)

#Data files to loop over
datafiles=['GBS_export.txt','GBS_export1.txt',
                 'GBS_export2.txt','GBS_export3.txt']
#datafiles=['GBS_BS.txt','GBS_BS1.txt','GBS_BS2.txt']

filepath="GBS_sim_export.txt"

#Extract Experimental Data
#Extract Data
states=[]
states=extract_data(filepath)

#Get CDF of data
CDF=get_CDF(states)

    
    
def comp_main(dim,SQargs,PHargs,BSargs,BeamS,allmodes1,allmodes2,ExpStates):
    
    COMP=[]
    Vexp_TOT=np.zeros(len(allmodes2))
    actual_TOT=np.zeros(len(allmodes2))
    
    for ite in range(4):
        #For each case create appropreate symplectic covariance
        #State 1: Beam Splitter is 50/50 and there is no phase
        if ite==0:
            print("State 1\n")
            #Fetch symplectic covariance for comparisons
            Vfinal=Vfinal_symplectic(dim,SQargs,PHargs,BSargs,BeamS)
        #State 2: Beam Splitter has zero reflectivity, no phase
        elif ite==1:
            print("State 2\n")
            (BSargs[3])[0]=0
            #Fetch symplectic covariance for comparisons
            Vfinal=Vfinal_symplectic(dim,SQargs,PHargs,BSargs,BeamS)
        #State 3: Fair beam splitter, phase=Pi/4
        elif ite==2:
            print("State 3\n")
            (BSargs[3])[0]=Pi/4
            PHargs[2]=Pi/4
            #Fetch symplectic covariance for comparisons
            Vfinal=Vfinal_symplectic(dim,SQargs,PHargs,BSargs,BeamS)
        #State 4: 0 reflectivity for beam splitter, phase Pi/4
        elif ite==3:
            print("State 4\n")
            (BSargs[3])[0]=0
            PHargs[2]=Pi/4
            #Fetch symplectic covariance for comparisons
            Vfinal=Vfinal_symplectic(dim,SQargs,PHargs,BSargs,BeamS)
    
            
        #Get sampled moments and actual moments
        Vexp, actual, diag=get_sampled_moments(allmodes1,allmodes2,ExpStates,Vfinal)
        
        
        #Compare
        comp=comparison(Vexp,actual)
    
        COMP.append(np.sum(comp))
        
        # plt.figure()
        # plt.subplot(2,1,1)
        # plt.plot(actual,'ro--',label="actual")
        # plt.plot(Vexp,'--bo',label="experimental")
        # x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
        # plt.xticks(x,labels=allmodes2,rotation="vertical")
        # plt.legend()
        
        # plt.subplot(2,1,2)
        # plt.plot(comp,'k.-',label='resid diag')
        # plt.xticks(x,labels=allmodes2,rotation="vertical")
        # plt.legend()
        # plt.show()
        
        Vexp_TOT=np.vstack([Vexp_TOT,Vexp])
        actual_TOT=np.vstack([actual_TOT,actual])

    Vexp_TOT=Vexp_TOT[1:,:]
    actual_TOT=actual_TOT[1:]
        
    return COMP , Vexp_TOT, actual_TOT

#Vexp_tot=np.zeros(len(allmodes2))
comparison_tot=np.zeros(4)
N=1
test1=np.zeros(10)
test2=np.zeros(10)
test3=np.zeros(10)
test4=np.zeros(10)

for i in range(N):
    print(i)
    #Sample PDF and get sampled CDF
    #np.random.seed(250)    
    X = run_sampling(CDF,Nruns=10e3)
    
    #Bring sampled PDF to correct format
    ExpStates=np.array(np.copy(states))    
    for i in range(len(states)):
        (ExpStates[i])[1]=X[i]
    
    print("--Sampling done--\n")
    
    
    COMP, Vexp, actual=comp_main(dim,SQargs,PHargs,BSargs,BeamS,allmodes1,allmodes2,ExpStates)

    
    (BSargs[3])[0]=Pi/4
    PHargs[2]=0
    comparison_tot=np.vstack([comparison_tot,COMP])
    
    test1=np.vstack([test1,Vexp[0]])
    test2=np.vstack([test2,Vexp[1]])
    test3=np.vstack([test3,Vexp[2]])
    test4=np.vstack([test4,Vexp[3]])
    
        
    #Vexp_tot=np.vstack([Vexp_tot,Vexp])   

comparison_tot=comparison_tot[1:,:]
test1=test1[1:,:]
test2=test2[1:,:]
test3=test3[1:,:]
test4=test4[1:,:]
#Vexp_tot=Vexp_tot[1:,:]
#%%
avCOMP=np.average(comparison_tot,0)

std=np.std(comparison_tot,0)
Conf_int=3.291*std/np.sqrt(N)
uplim=avCOMP[-1]+Conf_int[-1]*np.ones(4)
lowlim=avCOMP[-1]-Conf_int[-1]*np.ones(4)
x=[-0.4,0,2,3.4]
yerrr=Conf_int
#Plot
lab=["State 1","State 2","State 3","State 4"]
plt.figure()
#plt.hlines(uplim,-0.4,3.4,colors='r',zorder=4)#plt.hlines(lowlim,-0.4,3.4,colors='r',zorder=3)
plt.bar(lab,avCOMP,yerr=[Conf_int,Conf_int],color='cornflowerblue',zorder=1
        ,capsize=20,edgecolor="b",label="Res., 99.9% C. Int.")
#plt.errorbar(lab,avCOMP,yerr=yerrr,fmt='o', ecolor='k', capsize=25)
#plt.fill_between(lab,uplim,lowlim,alpha=0.6,color='tab:red',zorder=2)
plt.ylabel("Square Sum of Residuals")
plt.legend(loc=0, framealpha=1)
#plt.ylim(0,0.25)
plt.yscale('log')
plt.grid()
plt.show()


#%%

tests=[test1,test2,test3,test4]
fig, axs=plt.subplots(2,2,sharex=True,sharey=True)
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
for i, ax in enumerate(axs.flat):
    print(i)
    #Test1
    avVexp=np.average(tests[i],0)
    stdVexp=np.std(tests[i],0)
    Conf_int_V=3.291*stdVexp/np.sqrt(N)
    uplim=avVexp+Conf_int_V
    lowlim=avVexp-Conf_int_V
    
    print(axs.flat)
    ax.plot(actual[0],'or--',label="Theory State 1")
    ax.errorbar(x,avVexp,yerr=Conf_int_V,fmt='o',ecolor='k',capsize=7,label="Experimental")
    #plt.fill_between(x,uplim,lowlim,alpha=0.6,color='tab:red')
    x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
    ax.set_xticklabels(labels=allmodes2,rotation="vertical")
    ax.set_yscale('log')
    ax.grid()
    ax.legend()
    

#%%

# plot something
fig, axs = plt.subplots(3,3, figsize=(15, 8), sharex=True, sharey=True)

for i, ax in enumerate(axs.flat):
    print(axs.flat)
    ax.scatter(*np.random.normal(size=(2,200)))
    ax.set_title(f'Title {i}')

# set labels
plt.setp(axs[-1, :], xlabel='x axis label')
plt.setp(axs[:, 0], ylabel='y axis label')


#%%

#Test1
avVexp=np.average(test1,0)
stdVexp=np.std(test1,0)
Conf_int_V=3.291*stdVexp/np.sqrt(N)
uplim=avVexp+Conf_int_V
lowlim=avVexp-Conf_int_V

plt.figure()
plt.subplot(2,2,1)
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
plt.plot(actual[0],'or--',label="Theory State 1")
#plt.plot(Vexp,'--bo',label="experimental")
plt.errorbar(x,avVexp,yerr=Conf_int_V,fmt='o',ecolor='k',capsize=7,label="Experimental")
#plt.fill_between(x,uplim,lowlim,alpha=0.6,color='tab:red')
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
#plt.xticks(x,labels=allmodes2,rotation="vertical")
plt.yscale('log')
#plt.xlabel("Covariance Matrix Elements")
plt.ylabel("Covariance")
plt.grid()
plt.legend()
plt.tight_layout()

#Test2
avVexp=np.average(test2,0)
stdVexp=np.std(test2,0)
Conf_int_V=3.291*stdVexp/np.sqrt(N)
uplim=avVexp+Conf_int_V
lowlim=avVexp-Conf_int_V

plt.subplot(2,2,2)
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
plt.plot(actual[1],'or--',label="Theory State 2")
#plt.plot(Vexp,'--bo',label="experimental")
plt.errorbar(x,avVexp,yerr=Conf_int_V,fmt='o',ecolor='k',capsize=7,label="Experimental")
#plt.fill_between(x,uplim,lowlim,alpha=0.6,color='tab:red')
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
#plt.xticks(x,labels=allmodes2,rotation="vertical")
+#plt.xlabel("Covariance Matrix Elements")
plt.ylabel("Covariance")
plt.grid()
plt.legend()
plt.tight_layout()


#Test3
avVexp=np.average(test3,0)
stdVexp=np.std(test3,0)
Conf_int_V=3.291*stdVexp/np.sqrt(N)
uplim=avVexp+Conf_int_V
lowlim=avVexp-Conf_int_V

plt.subplot(2,2,3)
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
plt.plot(actual[2],'or--',label="Theory State 3")
#plt.plot(Vexp,'--bo',label="experimental")
plt.errorbar(x,avVexp,yerr=Conf_int_V,fmt='o',ecolor='k',capsize=7,label="Experimental")
#plt.fill_between(x,uplim,lowlim,alpha=0.6,color='tab:red')
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
plt.xticks(x,labels=allmodes2,rotation="vertical")
plt.xlabel("Covariance Matrix Elements")
plt.ylabel("Covariance")
plt.grid()
plt.legend()
plt.tight_layout()


#Test 4
avVexp=np.average(test4,0)
stdVexp=np.std(test4,0)
Conf_int_V=3.291*stdVexp/np.sqrt(N)
uplim=avVexp+Conf_int_V
lowlim=avVexp-Conf_int_V

plt.subplot(2,2,4)
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
plt.plot(actual[3],'or--',label="Theory State 4")
#plt.plot(Vexp,'--bo',label="experimental")
plt.errorbar(x,avVexp,yerr=Conf_int_V,fmt='o',ecolor='k',capsize=7,label="Experimental")
#plt.fill_between(x,uplim,lowlim,alpha=0.6,color='tab:red')
x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
plt.xticks(x,labels=allmodes2,rotation="vertical")
plt.xlabel("Covariance Matrix Elements")
plt.ylabel("Covariance")
plt.grid()
plt.legend()
plt.tight_layout()
        



