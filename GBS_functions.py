#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 20:43:59 2020

@author: frk
"""
import numpy as np
from numpy.linalg import multi_dot
from scipy.linalg import block_diag
import random as rd
import time


##############################################################################
#######################--STRAWBERRYFIELDS SECTION--###########################
##############################################################################

def beam_splitter_args(dim,BSargs,B):
    
    #Transmission, reflection arguments
    t_r_amplitudes = [(np.cos(q), np.exp(p*1j)*np.sin(q)) for q,p in BSargs]

    #Generate Unitary transformations of beam spleaters
    BSunitaries = [np.array([[t, -np.conj(r)], [r, t]]) for t,r in t_r_amplitudes]
    print(np.round(BSunitaries,2))
    UBS=np.diag(np.ones(dim))
    p=0
    UBSs=[]
    for i in B:
        if i==0:
            UBSs = block_diag(BSunitaries[p],[[1]],[[1]])
        elif i==1: 
            UBSs = block_diag([[1]], BSunitaries[p],[[1]])
        elif i==2:
            UBSs = block_diag([[1]],[[1]],BSunitaries[p])
        else:
            print("Error\n")
            break
        p=p+1
 
        UBS=multi_dot([UBS,UBSs])
                
    return UBS


######----Experimental Moments----######

def moments_experiment(modes,mstates):
    '''This function returns 1st and 2nd order moments'''
    mom=0
    #loop through outcomes
    for i in range(len(mstates)):
        #for each state and prob
        state=(mstates[i])[0]
        probmom=(mstates[i])[1]
        #calculate moment
        mod=[state[p-1] for p in modes]
        momprod=np.prod(mod)*probmom
        mom=mom+momprod
    return mom


##############################################################################
######################--SYMPLECTIC MATRICES SECTION--#########################
##############################################################################

def basic_sym_SQ(r):
    S_sq=np.array([[np.exp(-r),0],[0,np.exp(r)]])
    
    return S_sq

def sym_squeezer(SQargs,dim):
    '''This function constructs the total squeezing symplectic matrix'''
    S_SQ=[]
    #loop through beam splitters
    for i in range(dim):
        #create diagonal matrices
        front=int(i*2)
        back=int(dim*2-front-2)
        #fetch the basic form
        r=SQargs[i]
        S_sq=basic_sym_SQ(r)
        S_sq_gen=block_diag(np.diag(np.ones(front))
                            ,S_sq,
                            np.diag(np.ones(back)))
        S_SQ.append(S_sq_gen)

    SSQ=np.diag(np.ones(dim*2))
    for nu in S_SQ:
        SSQ=multi_dot([SSQ,nu])

    return SSQ;

def basic_sym_PH(th):
    S_th=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
    
    return S_th

def sym_phase_shifter(PHargs,dim):
    '''This function constructs the 
        phase shifting symplectic matrix'''
    S_PH=[]
    #loop through beam splitters
    for i in range(dim):
        #create diagonal matrices
        front=int(i*2)
        back=int(dim*2-front-2)
        #fetch the basic form
        th=PHargs[i]
        S_th=basic_sym_PH(th)
        S_th_gen=block_diag(np.diag(np.ones(front))
                            ,S_th,
                            np.diag(np.ones(back)))
        S_PH.append(S_th_gen)

    SPH=np.diag(np.ones(dim*2))
    for nu in S_PH:
        SPH=multi_dot([SPH,nu])

    return SPH;

def basic_sym_BS(chi):
    S_bs=np.array([[np.cos(chi),0 ,-np.sin(chi),0],
                   [0,np.cos(chi),0,-np.sin(chi)],
                   [np.sin(chi),0,np.cos(chi),0],
                   [0,np.sin(chi),0,np.cos(chi)]])
    
    return S_bs
    
def sym_beam_splitters(BSargs,BeamS,dim):
    '''This function constructs the total beam
        splitter symplectic matrix'''
    S_BS=[]
    #loop through beam splitters
    for i in range(len(BeamS)):
        #create diagonal matrices
        front=int(BeamS[i]*2)
        back=int(dim*2-front-4)
        #fetch the basic form
        chi=(BSargs[i])[0]
        S_bs=basic_sym_BS(chi)
        S_bs_gen=block_diag(np.diag(np.ones(front))
                            ,S_bs,
                            np.diag(np.ones(back)))
        S_BS.append(S_bs_gen)

    SBS=np.diag(np.ones(dim*2))
    for nu in S_BS:
        SBS=multi_dot([SBS,nu])

    return SBS


def Vfinal_symplectic(dim,SQargs,PHargs,BSargs,BeamS):
    '''This function returns the final covariance from
        symplectic matrices'''
    #Fetch symplectic matrices
    #Squeezing
    SS=sym_squeezer(SQargs,dim)
    SST=np.transpose(SS)
    #Phase shift
    SPH=sym_phase_shifter(PHargs,dim)
    SPHT=np.transpose(SPH)
    #Beam Splitters
    SBS=sym_beam_splitters(BSargs,BeamS,dim)
    SBST=np.transpose(SBS)
    
    
    #Tranformation of original covariance matrix
    V=np.diag(0.5*np.ones(8))
    Vfinal=multi_dot([SBS,SPH,SS,V,SST,SPHT,SBST])
    return Vfinal
    

######----Symplectic Moments----######

def moments_sym_1(l,Vfinal):
    '''Generates first order moment for specified mode, 
        from symplectic matrices'''
    mom1=0.5*(Vfinal[2*(l-1),2*(l-1)]+Vfinal[2*(l-1)+1,2*(l-1)+1]-1)
    return mom1

def moments_sym_2(l,k,Vfinal):
    '''Generates second order moment for specified mode, 
        from symplectic matrices'''
    mom2=0.25*(Vfinal[2*(l-1),2*(l-1)]*Vfinal[2*(k-1),2*(k-1)]
               +2*pow(Vfinal[2*(l-1),2*(k-1)],2)
               +Vfinal[2*(l-1)+1,2*(l-1)+1]*Vfinal[2*(k-1)+1,2*(k-1)+1]
               +2*pow(Vfinal[2*(l-1)+1,2*(k-1)+1],2)
               +Vfinal[2*(l-1),2*(l-1)]*Vfinal[2*(k-1)+1,2*(k-1)+1]
               +2*pow(Vfinal[2*(l-1),2*(k-1)+1],2)
               +Vfinal[2*(l-1)+1,2*(l-1)+1]*Vfinal[2*(k-1),2*(k-1)]
               +2*pow(Vfinal[2*(l-1)+1,2*(k-1)],2)
               -(Vfinal[2*(l-1),2*(l-1)]
                 +Vfinal[2*(l-1)+1,2*(l-1)+1]
                 +Vfinal[2*(k-1),2*(k-1)]
                 +Vfinal[2*(k-1)+1,2*(k-1)+1])
               +1)
    
    if l==k:
        mom2=mom2-0.25

    return mom2


##############################################################################
#################################--GENERAL--##################################
##############################################################################


def extract_data(filepath):
    '''Extracts data text file, from discrete probability distribution
        generated using strawberrie fields'''
    with open(filepath) as fp:
       mstates=[]
       line = fp.readline()
       while line:
           #print("{}".format(line.strip()))
           #Split line
           splitl=line.split(":")
           #Treat states
           ln=list((splitl[0])[1:-1])
           tint=[int(i) for i in ln]
           #Treat associated probabilities
           probs1=float(splitl[1])
           #append to list
           mstates.append([tint,probs1])
           #move to next line
           line = fp.readline()
    print("\n     --Data extracted--\n")
    return mstates;


def run_sampling(CDF1, Nruns=5000):
    '''Generate samples'''
    Xs = np.zeros_like(CDF1,dtype=int)
    for k in np.arange(Nruns):
        a = np.random.uniform(0,1)
        el=np.argmax(CDF1>=a)
        Xs[el] += 1

    return Xs/np.sum(Xs)


##############################################################################
#########################--SAMPLING SIMULATION--##############################
##############################################################################

def get_CDF(Expstates1):
    
    #Get empirical distribution to cumulative 
    PDF=[]
    for i in range(len(Expstates1)):
        probal=(Expstates1[i])[1]
        PDF.append(probal)
    CDF=[]
    CDF=np.cumsum(PDF)
    
    return CDF


def get_sampled_moments(allmodes1,allmodes2,Expstates1,Vfinal):
    
    #MIND THE ORDER
    momends2=[]
    for modes in allmodes2:
        mom=moments_experiment(modes,Expstates1)
        momends2.append(mom)
    
    momends1=[]
    for modes in allmodes1:
        mom=moments_experiment([modes],Expstates1)
        momends1.append(mom)
          
    Vexp=[]
    actual=[]
    diag=[]
    i=0
    for mode in allmodes2:
        n=2*(mode[0]-1)
        m=2*(mode[1]-1)
        tot=momends2[i]-momends1[int(n/2)]*momends1[int(m/2)]
        atot=0.5*(Vfinal[n,m]**2+Vfinal[n,m+1]**2
                  +Vfinal[n+1,m]**2+Vfinal[n+1,m+1]**2)
        adtot=0.5*(Vfinal[n,m]**2+Vfinal[n+1,m+1]**2)
        if mode[0]==mode[1]:
            atot=atot-0.25
            adtot=adtot-0.25
        i=i+1
        Vexp.append(tot)
        actual.append(atot)
        diag.append(adtot)
        
    return Vexp, actual, diag

def comparison(Vexp,actual):
    
    #Compare Vexp and Actual
    comparison=[]
    for j in range(len(Vexp)):
            summation=(Vexp[j]-actual[j])**2
            comparison.append(summation)
    #print(comparison)
    return comparison

# plt.figure()
# plt.subplot(2,1,1)
# plt.plot(actual,'ro--',label="actual")
# plt.plot(Vexp,'--bo',label="experimental")
# x=np.linspace(0,len(allmodes2)-1,len(allmodes2))
# plt.xticks(x,labels=allmodes2,rotation="vertical")
# plt.legend()

# plt.subplot(2,1,2)
# diff1=np.array(actual)-np.array(Vexp)
# plt.plot(diff1,'k.-',label='resid diag')
# plt.xticks(x,labels=allmodes2,rotation="vertical")
# plt.legend()
# plt.show()    
# #Bring Data to correct format
# for i in range(len(Expstates1)):
#     (Expstates1[i])[1]=X[i]
#Perfect data: mstates
#Sampled data: Expstates

# #Write Data
# #Open file to write data
# file = open("GBS_export_sampling.txt", "w")
# for i in range(len(Expstates)):
#     prob=(Expstates[i])[1]
#     state=((Expstates)[i])[0]
#     #print("|{}>: {}".format("".join(str(j) for j in state), prob))
#     #write to file
#     file.write("|{}>: {}".format("".join(str(j) for j in state), prob))

# #Close file
# file.close()

#Calculate covariance measure V_13+V_24
#<n_1 n_2> -<n_1><n_2>=1/2(V_13^2+V_24^2+V_14^2+V_23^2)
#Therefore we need
#<n_1 n_1>
#<n_1 n_2> 
#<n_1 n_3> 
#<n_1 n_4> 
#<n_2 n_2> 
#<n_2 n_3> 
#<n_2 n_4> 
#<n_3 n_3> 
#<n_3 n_4> 
# #<n_4 n_4>  

# #functions to get the necessary output combinations
# #i.e |0011>, |0100> etc.
# def to_base(s, b):
#     '''takes any number as an input (s) and returns it
#     in any basis (b)''' 
    
#     BS="012345678910"
#     res = ""
#     while s:
#         res+=BS[s%b]
#         s//= b
#     return res[::-1] or "0"

# def state_comb(dim, base):
#     '''dim is the dimension of the interferometer and
#     base is used to generate all the combinations given that base
#     i.e dim 2 base 2 would return |00> |01> |10> |11>'''
#     count=0
#     measured=[]
#     for i in range(10000):
#         t=list(map(int,to_base(i,base)))
#         if len(t)<dim:
#             add=dim-len(t)
#             t=list(map(int,np.zeros(add)))+t 
#         elif len(t)>dim:
#             break;
        
#         if sum(t)<=base**dim and sum(t)%2==0:
#             count+=1
#             measured.append(t)
#     return measured;
