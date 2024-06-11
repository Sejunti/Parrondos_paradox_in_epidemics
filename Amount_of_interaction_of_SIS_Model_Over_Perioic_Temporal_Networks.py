import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sp
from scipy import linalg as la
import math
import matplotlib.ticker as ticker
from SIS_Model_Over_Perioic_Temporal_Networks import *  


def amount_of_interaction_of_static_network(A):
    """
    Function gives the amount of interaction of the static network.
    Parameter:
              A(numpy array): Adjacency matrix
    Returns:
            s1(float): Amount of interaction of the static network
    """
    s1=np.sum(A)
    return s1


###############################################################

def amount_of_interaction_of_switching_network(As, durations,T):
    """
    Function gives the amount of interaction for the switching network.
    Parameter:
              As(numpy array): All Adjacency matrices
              durations(list): Duration of each static network in one period
              T(integer): Period           
    Returns:
            s2(float): Amount of interaction of the switching network.
    """
    T=sum(durations)
    interaction = durations[0]*np.sum(As[0])+durations[1]*np.sum(As[1])
    s2= interaction/T
    return s2






