import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sp
from scipy import linalg as la
import math
from Simulation_of_SIS_Model_Over_Perioic_Temporal_Network import * 

def fraction_anti_phase(z):
    """
    Function measures the extent of anti-phase oscillation by the fraction of time.
    Parameter:
              z(numpy array): State of the system at each time step
    Returns:
            q(float): Fraction of time  
    """

    times, z = get_trajectory(beta,mu,As,periods,durations,z0)
    time_interval = (times >= t_start) & (times <= t_end)
    z_interval = z[time_interval]
    q=(np.sum(np.diff(z_interval[:, 0]) * np.diff(z_interval[:, 1]) < 0)) / len(z_interval[:, 0])
    return q







