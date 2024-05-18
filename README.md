# Parrondos_paradox_in_epidemics
# Parrondos_paradox_in_epidemics

When you use the code, please cite the following  paper:
[Parrondo’s paradox in susceptible-infectious-susceptible dynamics over periodic temporal networks]


This repository provides the code written in Python for the Parrondo's paradox ananlysis in SIS model over periodic temporal network.


# How to use the code?
We divided the whole analysis into a few parts and for each part provided a Python (.py) file containing relevant functions and required libraries. 

- `SIS_Model_Over_Perioic_Temporal_Networks.py` contains  some functions. These functions can be used to determine the the largest eigenvalue of a matrix and the largest Floquet exponent. 

    - Function 'extract_max_eval_vary_betas(betas,mu,A)'  use for computing the largest eigenvalue of the static network which takes infection rate ($\beta$), recovery rate ($\mu$) and the adjacency matrix, $\textbf{A}$ as input and gives the largest eigenvalue as output. 
    - Function 'get_max_flo_exponent(beta,T)'  computes  the largest Floquet exponent, $\lambda_{\text{F}}$. The input of this function is infection rate, $\beta$ and the period, $T$ and gives $\lambda_{\text{F}}$ as output. 
    
    
- `Simulation_of_SIS_Model_Over_Perioic_Temporal_Network.py` has functions  for simulating the linear, nonlinear and agent-based simulation of SIS model on switching networks. 
    - Function `get_trajectory(beta,mu,As,periods,durations,z0)` simulates the linear SIS model on switching networks. The inputs are:
        - "$\beta$": infection rate,
        - "$\mu$": recovery rate,
        - "$\textbf{As}$": list of all adjacency matrices,
        - "periods": period of the system,
        - "durations": list of the time interval,
        - "$z0$": initial condition 
    - The function gives state of the system at each time step of linear SIS model as output. 
    - Function 'get_trajectory_2(beta, mu, As, periods, durations, z0)' simulates the nonlinear SIS model on switching networks which takes same input as function `get_trajectory(beta,mu,As,periods,durations,z0)' and the output of the function is the state of the system at each time step of nonlinear SIS model.
    - Function  'fraction_of_infected_individual(number_of_simulations,p,A1,A2,beta, mu,T,states)' is used to find the fraction of infected individual at each community for agent-based simulation of SIS model simulation on switching networks. The input functions are:
              -'number_of_simulations': number of simulation, 
              -'$p$' Fraction of initially infected individual in each block,
              -'$\textbf{A}^{(1)}$': first adjacency  matrix, 
              -'$\textbf{A}^{(2)}$': second adjacency  matrix,
              -'$\beta$':infection rate,
              -'$\mu$': recovery rate,
              -'$T$': simulation duration
              -'states':  node state   
  
- `Amount_of_interaction_of_SIS_Model_Over_Perioic_Temporal_Networks.py` has functions for finding the relationship between amount of interaction and epidemic spreading. 
    - Function `amount_of_interaction_of_static_network(A)` gives the amount of interaction of the static network and takes the adjacency matrix, $\textbf{A}$ as input.
    - Function `amount_of_interaction_of_switching_network(As, durations,T)` computes the amount of interaction for  switching  network and the list of all adjacency matrices, $\textbf{As}$, list of the time interval, durations and period, $T$ are the input of this function. 
    
- 'Generality_of_the_Parrondo’s_paradox.py' is used to measure the epidemic threshold for both the static and switching network. 
    - Function 'power_method(A,tolerance)' is used for finding the largest eigenvalue of the matrix by Power method whose input is $\textbf{A}$ the adjacency matrix and the tolerance.
    - Function 'get_x_intercept(function,tolerance)' computes the root of the system by root finding algorithm which takes function and tolerance as input.
    - Function 'compute_epidemic_threshold(As_list)' measures the epidemic threshold for both static and switching network which takes the list of all adjacency matrices as input and returns the epidemic threshold for static networks and switching networks.
    
-'Perturbation_theory_for_the_largest_Floquet_exponent.py'  has function for finding the first order approximation of largest Floquet exponent. 
  - Function 'approximate_matrix(beta,mu,T,barA)' computes the first-order approximation of  matrix $\mathcal{T}$ and the input of this function are the infection rate, recovery rate, period and $\overline{\textbf{A}}$ is the time-averaged adjacency matrix.
  - Function 'approximate_largest_eigenvalue(beta,mu,T,barA)'measures the first-order approximation of the largest eigenvalue of $\mathcal{T}$ and input of this function same as input of the function  'approximate_matrix(beta,mu,T,barA)'.
  - Function 'approximate_lamda_F(beta,mu,T,barA)' computes the first order approximate largest Floquet exponent,  $\lambda_{\text{F}}$ and the input of this function are the infection rate, recovery rate, period and $\overline{\textbf{A}}$ is the time-averaged adjacency matrix.

-'Parrondo’s_paradox_and_anti-phase_oscillation.py'  computes the anti-phase oscillation.
   - Function 'fraction_anti_phase(z, times, t_start, t_end)' measures the anti-phase oscillation in one period as output and  the input functions are:
      - '$z$': state of the system at each time step of linear SIS model
      - 'times': time step
      - 't_start': period start value
      - 't_end': period end value
   
       

               
    
    
    
