#!/usr/bin/env python
# coding: utf-8




import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('pylab', 'inline')
from scipy.integrate import odeint
import sympy as sp
sp.init_printing()
from scipy import linalg as la
import math


#################################################################################################


def M_matrix(beta,mu,A):
    """
    Function for create M matrix for linear SIS epidemic model.
    Parameters:
              beta(float):Infection rate
              mu(float): Recovery rate
              A(numpy array): Adjacency matrix 
    Returns:
            M(numpy array): Matrix M
    """
    M = beta*A - mu*np.eye(len(A))     
    return M
  
#################################################################################################
def create_Ms(beta,mu,As): 
    """
    Function for create M matrices for linear SIS epidemic model over switching temporal network.
    Parameters:
              beta(float):Infection rate
              mu(float): Recovery rate
              As(numpy array): All Adjacency matrices 
    Returns:
            Ms(numpy array): Matrices M
    """
    N = len(As[0]) # number of nodes
    T = len(As) # number of adjacency matrices
    Ms = []
    for t in range(T):
        Ms.append(M_matrix(beta,mu,As[t]))
    return Ms
#################################################################################################
def extract_max_eval_vary_betas(betas,mu,A):
    """
    Function for computing largest eigenvalue for M matrix. 
    Parameters:
              betas(numpy array): Infection rate
              mu(float): Recovery rate
              A(numpy array): Adjacency matrix 
    Returns:
            max_evals(float): Largest eigenvalue of matrix M.
    """
    max_evals = np.zeros(len(betas))
    for t,beta in enumerate(betas):
        M = M_matrix(beta,mu,A)
        e1,e2 = eig(M)
        max_evals[t] = max(e1)
    return max_evals


#################################################################################################


def monodromy_matrix(Ms,durations):
    """
    Function for creating monodromy matrix. 
    Parameters:
              Ms(numpy array): Matrices M
              durations(list): Duartion of the time interval
    Returns:
            Monodromy(numpy array): Monodromy matrix $\hat{M}$
    """

    expMdts = []
    for t,M in enumerate(Ms):
        lam,U = eig(M) # eigenvalue and eigenvector of matrix M
        sigma = np.diag(np.exp(durations[t] * lam))
        expMdt = dot(U,sigma)
        expMdt = dot(expMdt,inv(U))
        expMdts.append(expMdt)
        
    Monodromy = np.eye(len(Ms[0]))   
    for t,expMdt in enumerate(expMdts):
        Monodromy = dot(expMdt,Monodromy)
        
    return Monodromy


#################################################################################################


def get_max_flo_mult(beta):
    """
    Function for creating largest Floquet multiplier. 
    Parameter:
              beta(float): Infection Rate
    Returns:
            max(Floquet_multipliers)(float): Largest Floquet multiplier
    """
    Ms = create_Ms(beta,mu,As)
    mono = monodromy_matrix(Ms,durations)
    floquet_multipliers,floquet_vectors = eig(mono)
    return max(floquet_multipliers)
 #################################################################################################

def get_max_flo_exponent(beta,T):
    """
    Function for computing largest Floquet exponent. 
    Parameter:
              beta(float): Infection Rate
              T(Integer): Period
    Returns:
            lamda_F(float): Largest Floquet exponent
    """
    max_flo_multiplier = zeros(len(betas))
    T=sum(durations) # period of the system
    for i,beta in enumerate(betas):
         max_flo_multiplier[i] = get_max_flo_mult(beta)
    lamda_F=log(max_flo_multiplier)/T
    return lamda_F







