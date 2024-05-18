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
from scipy.linalg import eig
from SIS_Model_Over_Perioic_Temporal_Networks import *  




def approximate_matrix(beta,mu,T,barA):
    """
    Function for finding the first-order approximation of  matrix T. 
    Parameters:
              beta(float):Infection rate
              mu(float):Recovery rate
              T(integer): Period 
              barA(numpy array): Average adjacency matrix
    Returns:
           barT(numpy array): First order approximation of  matrix T
    """ 
    c=matrix_exponetial(-T*mu*I)
    h=I+(beta*T*barA)
    barT=np.dot(h,c)
    return barT


###################################################################

def approximate_largest_eigenvalue(beta,mu,T,barA):
    """
    Function for finding the first-order approximation of the largest eigenvalue of T. 
    Parameters:
              beta(float):Infection rate
              mu(float):Recovery rate
              T(integer): Period 
              barA(numpy array): Average adjacency matrix
    Returns:
           approximate_max_eigen_value(float): First order approximation of largest eigenvalue of T
    """ 
    a=np.exp(-T*mu)
    eigen_value,eigen_vectors = eig(barA)
    max_eigen_value=max(eigen_value)
    approximate_max_eigen_value=a+beta*max_eigen_value*T*a
    return approximate_max_eigen_value


###################################################################



def approximate_lamda_F(beta,mu,T,barA):
    """
    Function for finding the first-order approximation of the Floquet exponent. 
    Parameters:
              beta(float):Infection rate
              mu(float):Recovery rate
              barA(numpy array): Average adjacency matrix
    Returns:
           approximate_lamda_F(float): First order approximation of the Floquet exponent
    """ 
    eigen_value,eigen_vectors = eig(barA)
    max_eigen_value=max(eigen_value)
    approximate_lamda_F=-mu+beta*max_eigen_value
    return approximate_lamda_F





