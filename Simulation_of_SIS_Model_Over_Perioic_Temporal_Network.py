#!/usr/bin/env python
# coding: utf-8



import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('pylab', 'inline')
from random import *
from scipy import sparse
from scipy.integrate import odeint
from scipy import linalg as la
import math
from random import sample
from numpy.random import Generator, PCG64


#################################################################################################


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
  
#################################################################################################
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


#################################################################################################


def matrix_exponetial(a):
    """
    Function for finding the exponential of a matrix.  
    Parameter:
              a(numpy array): Matrix 
    Returns:
             s(float): Matrix exponential 
    """
    lam,u=eig(a)
    s1=np.diag(np.exp(lam))
    s2=np.dot(u,s1)
    s=np.dot(s2,inv(u)) 
    return s


#################################################################################################


# Linear Simulation of SIS model 
def integrate_interval(M,times,start): 
    """
    Function for integrating the linear SIS model over the interval.  
    Parameters:
              M(numpy array): Matrix 
              times(numpy array): Time interval
              start(numpy array):Initial Condiiton
    Returns:
             u(numpy array): State of the system at each time step
    """


    u = [ list(np.dot(matrix_exponetial( (time-times[0]) *M), start)) for time in times]
    return u
#################################################################################################

def get_trajectory(beta,mu,As,periods,durations,z0):
    """
    Function for generating the trajectory of the system over time.  
    Parameters:
              beta(float):Infection rate
              mu(float): Recovery rate
              As(numpy array): All Adjacency matrices 
              periods(Integer): Number of period 
              durations(list): Duartion of the time interval
              z0(numpy array): Initial condition
    Returns:
             times(numpy array): Time step
             z(numpy array): State of the system at each time step
    """
    
    As_temp = As*periods
    durations = durations*periods
    transition_times = list(np.cumsum(np.array(durations)))
    transition_times.insert(0,0)

    # initial condition
    u0 = z0
        
    Ms = create_Ms(beta,mu,As_temp)

    
    for s in range(len(Ms)):

        start_time = transition_times[s]
        end_time = transition_times[s+1]
        n = 1000*(end_time-start_time)
        interval_times = np.linspace(start_time,end_time,n)

        u = integrate_interval(Ms[s],interval_times,u0)
        u0 = u[-1]
        if s==0:
            times = interval_times
            z = u
        if s>0:
            times = np.concatenate((times,interval_times))
            z = np.concatenate((z,u))
        
    return times,z


#################################################################################################


#Nonlinear Simulation of the SIS model 

def integrate_interval_2(M,As, times, start):
    """
    Function for integrating the nonlinear SIS model over the interval.  
    Parameters:
              M(numpy array): Matrix 
              As(numpy array): All Adjacency matrices 
              times(numpy array): Time interval
              start(numpy array):Initial Condiiton
    Returns:
             u(numpy array): State of the system at each time step
    """
    if mod(t,7)<5:
        A=As[0]
    else:
        A=As[1]
    def diff_eq(u, t): #Function for defining the nonlinear SIS model.  
        return np.dot(M, u) - beta * np.dot(np.dot(np.diag(u), A), u)

    u = odeint(diff_eq, start, times)
    return u
#################################################################################################
def get_trajectory_2(beta, mu, As, periods, durations, z0):
    """
    Function for generating the trajectory of the system over time.  
    Parameters:
              beta(float):Infection rate
              mu(float): Recovery rate
              As(numpy array): All Adjacency matrices 
              periods(Integer): Number of period 
              durations(list): Duartion of the time interval
              z0(numpy array): Initial condition
    Returns:
             times(numpy array): Time step
             z(numpy array): State of the system at each time step
    """
    As_temp = As * periods
    durations = durations * periods
    transition_times = list(np.cumsum(np.array(durations)))
    transition_times.insert(0, 0)

    # initial condition
    u0 = z0

    Ns = create_Ns(beta, mu, As_temp)
#     Us = create_Us(u0, As_temp)

    for s in range(len(Ns)):
        start_time = transition_times[s]
        end_time = transition_times[s + 1]
        n = 1000 * (end_time - start_time)
        interval_times = np.linspace(start_time, end_time, n)

        u = integrate_interval_2(Ms[s], As_temp[s], interval_times, u0)  # Pass As_temp[s] as A
        u0 = u[-1]
        if s == 0:
            times = interval_times
            z = u
        if s > 0:
            times = np.concatenate((times, interval_times))
            z = np.concatenate((z, u))

    return times, z

#################################################################################################


# Agent based model simulation by Gillespie Algorithm
#--- Set up PRNG: ---
seed= 42                     # Set seed of PRNG state 
rg = Generator(PCG64(seed))  # Initialize bit generator (here PCG64) with seed


#################################################################################################


def calculate_node_propensities(A, states, beta, mu):
    """
    Function for generating the node propensities.  
    Parameters:
              A(numpy array): Adjacency matrix 
              states:  Node state 
              beta(float):Infection rate
              mu(float): Recovery rate
    Returns:
            Lambda_list: Total node rate 
    """       
    N = len(states)
    Lambda_list = np.zeros(N)
    for i in range(N): 
        if states[i] == 0:  # If node i is susceptible calculate total propensity for infection
            neighbor_states = states[A[i]>0] 
            lambda_i = beta * np.sum(neighbor_states) # Total node rate from number of neighbors that are infectious
            Lambda_list[i] = lambda_i 
        elif states[i]==1:  # If node is infectious set propensity to mu
            Lambda_list[i] = mu
    return Lambda_list


#################################################################################################


def draw_next_event_direct(Lambda_list, states):
    '''
    Function for selecting the reaction channel and the waiting time until the event.
    Parameters:
              Lambda_list: Total node rate  
              states:  Node state      
    Returns:
             i:Selected reaction channel
             tau: The waiting time 
    '''
    Lambda = np.sum(Lambda_list)
    u1, u2 = np.random.rand(2) # Draw two uniform random variates from (0,1):

    tau = -np.log(1. - u1) / Lambda # Draw waiting time

    target_sum = u1 * Lambda # Select reaction by linear search
    sum_i = 0
    for i in range(len(states)):
        
        sum_i += Lambda_list[i]
        
        if sum_i >= target_sum:
            break

    return i,tau


#################################################################################################


def update_states(A, states, X1, X2, Lambda_list, i):
    '''
    Function for update the states of nodes and propensities after an event
    Parameters:
              A(numpy array): Adjacency matrix 
              states:  Node state 
              X1= First community or first block in SBM
              X2= Second community or second block in SBM
              Lambda_list: Total node rate  
              i:Selected reaction channel
               
    Returns:
            X1= Updated first community or first block in SBM
            X2= Updated second community or second block in SBM
            Lambda_list: Updated total node rate 
            states:  Updated node state 
    '''
    state_before = states[i] 
    N = len(states)

    if state_before == 0:  # If state_before was S, update to I
        if i < int(N/2): # Update state counts
            X1[0] -= 1   
            X1[1] += 1
        else:
            X2[0] -= 1
            X2[1] += 1

        states[i] = 1  # Update node state to infectious    
        Lambda_list[i] =mu  # Update i's propensity
        susceptible_neighbors = np.where((A[i] > 0) & (states == 0))[0] # Update propensities of i's neighbors
        if len(susceptible_neighbors) > 0:
            Lambda_list[susceptible_neighbors] = Lambda_list[susceptible_neighbors]+beta


    else:  # Else  state_before was I, update to S
        if i < int(N/2):
            X1[1] -= 1
            X1[0] += 1
        else:
            X2[1] -= 1
            X2[0] += 1
        susceptible_neighbors = np.where((A[i] > 0) & (states == 0))[0]# Update propensities of i's neighbors
        if len(susceptible_neighbors) > 0:
            Lambda_list[susceptible_neighbors] = Lambda_list[susceptible_neighbors] - beta
        
        states[i] = 0  # Update node state to susceptible
        Lambda_list[i] = beta * len(np.where((A[i] > 0) & (states == 1))[0])# Update i's propensity

    return X1, X2, Lambda_list, states


#################################################################################################


def direct_method_SIS_graph(A1,A2, beta, mu, T, states):
    '''
    Function for direct method by Gillespie algorithm over SIS model on network.
    Parameters:
              A1(numpy array): First adjacency  matrix 
              A2(numpy array): Second adjacency  matrix 
              beta(float):Infection rate
              mu(float): Recovery rate
              T(integer)=Simulation duration
              states:  Node state 
               
    Returns:
            X_t(numpy array)= States at time t 
    '''
    N = len(A1)  
    I1 = np.sum(states[:int(N/2)]) # Total number of infectious in first block
    S1 = len(states)/2-I1 # Total number of susceptible in first block
    I2 = np.sum(states[int(N/2):]) # Total number of infectious in second block
    S2 = len(states)/2-I2 # Total number of susceptible in second block
    t = 0 # initial time 
    

    X_t = []# Vector to save temporal evolution of state numbers over time
    X_t.append([t, S1, I1, S2, I2])
    
    while t < T: # Keep drawing events until t >= T
        if mod(t,7)<5:
            A=A1
        else:
            A=A2
        Lambda_list = calculate_node_propensities(A, states, beta, mu)
    
        if np.isclose(np.sum(Lambda_list), 0.):# Check for Lambda == 0 (no more reactions can happen)
            X_t.append([T, S1, I1, S2, I2])
            break

        i, tau = draw_next_event_direct(Lambda_list,  states) #Direct event sampling step

        t += tau

        if mod(t,7)<5:
            A=A1
        else:
            A=A2
        [S1, I1], [S2, I2], Lambda_list, states = update_states(A, states, [S1, I1], [S2, I2], Lambda_list, i)# Update states
    
        X_t.append([t, S1, I1, S2, I2]) #Save current state numbers to X_t

    return np.array(X_t).transpose()


#################################################################################################


def initialize(A, p): 
    '''
    Function for the initialization of the network.
    Parameters:
              A(numpy array): Adjacency  matrix 
              p(float): Fraction of initially infected individual 
    Returns:
            states: Node state  
    '''

    N = len(A)
    infec_ind1 = sample(range(int(N/2)), int(p*N/2)) 
    infec_ind2 = sample(range(int(N/2), N), int(p*N/2)) 
    states = np.zeros(N)
    states[infec_ind1 + infec_ind2] = 1
    return states


#################################################################################################


def fraction_of_infected_individual(number_of_simulations,p,A1,A2,beta, mu,T,states):
    '''
    Function for finding the fraction of infected individual in each block. 
    Parameters:
              number_of_simulations(integer): Number of simulation 
              p(float): Fraction of initially infected individual 
              A1(numpy array): First adjacency  matrix 
              A2(numpy array): Second adjacency  matrix 
              beta(float):Infection rate
              mu(float): Recovery rate
              T(integer): Simulation duration
              states:  Node state   
    Returns:
            X_array(numpy array): Current Node state  
            
    '''
    X_array = []
    for q in range(number_of_simulations):
        states = initialize(A1, p)
        X_t = direct_method_SIS_graph_linear(A1,A2, beta, mu, T, states)    
        X_array.append(X_t)
    return X_array






