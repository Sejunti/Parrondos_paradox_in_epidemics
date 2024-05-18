#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy.optimize
import sympy as sp
from scipy import linalg as la
import math
import numpy as np
from random import *
from scipy import sparse
import networkx as nx
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix




def power_method(A,tolerance):
    """
    Function for finding the largest eigenvalue of the matrix by Power method. 
    Parameters:
              A(numpy array): Adjacency matrix 
              tolerance(float): Tolerance 
    Returns:
            beta(float): Infection rate
            z(float): Largest eigenvalue
    """
    x0 = ones(shape(A)[0])
    z_old = 0
    x=x0
    step_change = 10
    while step_change > tol:
        z = x/norm(x,2)
        step_change = norm(z-z_old,inf)
        z_old = z
        x = A@z
        beta = z@x

    return beta,z 


###############################################################



def M_matrix(beta,mu,A):
    """
    Create M matrix for linear SIS epidemic model.
    Parameters:
              beta(float):Infection rate
              mu(float): Recovery rate
              A(numpy array): Adjacency matrix 
    Returns:
            M(numpy array): Matrix M
    """
    M = beta*A - mu*np.eye(len(A))     
    return M
  
###############################################################

def create_Ms(beta,mu,As): 
    """
    Create M matrices for linear SIS epidemic model over switching temporal network.
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

###############################################################



def T_matrix(US,invUS,lams,beta,durations):
    """
    Function for find matrix T by singular value decomposition for the SIS model over switching temporal network.
    Parameters:
              US: Eigenvectors
              invUS: Inverse eigenvectors
              lams: Eigenvalue
              beta(float):Infection rate
              durations(list): Duartion of the time interval
              As(numpy array): All Adjacency matrices 
    Returns:
           Matrix_T(numpy array): Matrix T
    """

    expMdts = []
    for t,U in enumerate(US):
        sigma = np.diag(np.exp(durations[t] * (beta*lams[t]-mu)))
        expMdt = dot(U,sigma)
        expMdt = dot(expMdt,invUS[t])
        expMdts.append(expMdt)
    Matrix_T = np.eye(len(US[0]))   
    for t,expMdt in enumerate(expMdts):
        Matrix_T = dot(expMdt,Monodromy)
        
    return Matrix_T


###############################################################



def get_x_intercept(function,tolerance):
    """
    Function for finding the root by root finding algorithm.
    Parameters:
              function: Given function
              tolerance(float): Tolerance value
    Returns:
           sol.root(float): Root of the function
    """
    sol = scipy.optimize.root_scalar(fun, bracket=[0.0, 0.1], method='brentq',xtol=tolerance)
    return sol.root


###############################################################



def compute_epidemic_threshold(As_list):
    """
    Function for computing the epidemic threshold for both static and switching network. 
    Parameter:
              As_list(numpy array): List of all adjacency matrices
              
    Returns:
           sol.root(float): Root of the function
           M1_list: Epidemic threshold for first static network 
           M2_list: Epidemic threshold for second static network 
           Fl_list: Epidemic threshold for switching network
    """
    M1_list = [] # Epidemic threshold for first static network 
    M2_list = [] # Epidemic threshold for second static network 
    Fl_list = [] # Epidemic threshold for switching network

    for As in As_list:
        US=[]
        invUS=[]
        lams=[]
    for A in As:
        lam,U=eig(A)
        lams.append(lam)
        US.append(U)
        invUS.append(inv(U))
   
    def lam_3(beta):
        M_beta=monodromy_matrix_2(US,invUS,lams,[beta],dts)
        e,_ = power_method(M_beta)
        lamda=log(e)/T
        return lamda 
    e1,_ = power_method(As[0])
    e11,_ = power_method(As[1])
    M1_list.append(mu/(e1))
    M2_list.append(mu/(e11))
    Fl_list.append(get_x_intercept(lam_3))
    return M1_list,M2_list,Fl_list







