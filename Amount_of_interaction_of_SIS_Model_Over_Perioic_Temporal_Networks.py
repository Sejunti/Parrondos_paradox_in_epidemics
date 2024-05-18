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
import matplotlib.ticker as ticker
from SIS_Model_Over_Perioic_Temporal_Networks import *  




def amount_of_interaction_of_static_network(A):
    """
    Function for computing the amount of interaction of the static network. 
    Parameter:
              A(numpy array): Adjacency matrix
    Returns:
            s1(float): Amount of interaction of the static network.
    """
    upper_triangular_entry = np.triu(A)
    s1=np.sum(upper_triangular_entry)
    return s1


###############################################################


def amount_of_interaction_of_switching_network(As, durations,T):
    """
    Function for computing the amount of interaction of the static network. 
    Parameter:
              As(numpy array): All Adjacency matrices
              durations(list): Duartion of the time interval
              T(interger): Period 
              
    Returns:
            s2(float): Amount of interaction of the switching network.
    """
    T=sum(durations)
    interaction = [amount_of_interaction_of_static_network(A) for A in As]
    s2= np.mean(interaction) * T / len(As)
    return s2





