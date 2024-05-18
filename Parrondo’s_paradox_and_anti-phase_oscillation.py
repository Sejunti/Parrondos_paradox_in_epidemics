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
from Simulation_of_SIS_Model_Over_Perioic_Temporal_Network import * 


#########################################################################


def fraction_anti_phase(z, times, t_start, t_end):
    """
    Function for measure  the anti-phase oscillation.
    Parameters:
              z(numpy array): State of the system at each time step
              times(numpy array): Time step
              t_start(integer): Start of time interval
              t_end(integer): End of time interval
    Returns:
            q(float): Anti-phase oscillation 
    """

    time_interval = (times >= t_start) & (times <= t_end)
    times, z = get_trajectory(beta,mu,As,periods,durations,z0)
    z_interval = z[time_interval]
    q=(np.sum(np.diff(z_interval[:, 0]) * np.diff(z_interval[:, 1]) < 0)) / len(z_interval[:, 0])
    return q







