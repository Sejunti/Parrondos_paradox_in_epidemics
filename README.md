# Parrondos_paradox_in_the susceptible-infectious_susceptible_(SIS)_model_on_temporal_networks


When you use the code, please cite the following  paper:
[Parrondo’s paradox in susceptible-infectious-susceptible dynamics over periodic temporal networks]


This repository provides Python code for analysis of the Parrondo's paradox in the susceptible-infectious-susceptile (SIS) dynamics over periodically switching temporal networks.


# How to use the code?
We divided the whole analysis into a few parts and for each part provided a Python (.py) file containing relevant functions and required libraries. 

- `SIS_Model_Over_Perioic_Temporal_Networks.py` contains some functions. These functions can be used to compute the largest eigenvalue of a matrix and the largest Floquet exponent. 

    - Function 'extract_max_eval_vary_betas(betas,mu,A)' computes the largest eigenvalue of the static network. Its inputs are the infection rate (betas: $\beta$), recovery rate (mu: $\mu$), and adjacency matrix, $\textbf{A}$. Its output is the largest eigenvalue, $\lambda_{\text{max}}$. 
    - Function 'get_max_flo_exponent(beta,T)' computes the largest Floquet exponent, $\lambda_{\text{F}}$. Its inputs are the infection rate and the period, $T$. Its output is $\lambda_{\text{F}}$. 
    
- `Simulation_of_SIS_Model_Over_Perioic_Temporal_Network.py` has functions for running the linear deterministic SIS dynamics (i.e., linearized individual-based approximation), nonlinear deterministic SIS dynamics (i.e., individual-based approximation), and stochastic agent-based SIS dynamics on switching networks. 
    - Function `get_trajectory(beta,mu,As,periods,durations,z0)` simulates the linear deterministic SIS dynamics on periodic switching networks. Its output is the state of the system at each time step and the inputs are:
        - $\beta$: infection rate (beta),
        - $\mu$: recovery rate (mu),
        - $\textbf{As}$: list of all adjacency matrices,
        - 'periods': period of the dynamical system,
        - 'durations': list of the time interval,
        - $z0$: initial condition.  
    - Function `get_trajectory_2(beta, mu, As, periods, durations, z0)` simulates the nonlinear deterministic SIS dynamics on periodic switching networks. Its inputs and outputs are the same as those of `get_trajectory(beta,mu,As,periods,durations,z0)`.
    - Function  `fraction_of_infected_individual(number_of_simulations,p,A1,A2,beta, mu,T)` finds the fraction of infectious nodes in each block for stochastic agent-based simulations on periodic switching networks. Its inputs are:
        - 'number_of_simulations': number of simulation, 
       -  $p$: Fraction of initially infectious nodes in each block,
       -  $\textbf{A}^{(1)}$: first adjacency matrix, 
       -  $\textbf{A}^{(2)}$: second adjacency matrix,
        - $\beta$: infection rate (beta),
       - $\mu$: recovery rate (mu),
       - $T$: simulation duration
      
- `Amount_of_interaction_of_SIS_Model_Over_Perioic_Temporal_Networks.py` has functions for finding the relationship between amount of interaction and epidemic spreading. 
    - Function `amount_of_interaction_of_static_network(A)` gives the amount of interaction of the static network. Its input is the adjacency matrix, $\textbf{A}$.
    - Function `amount_of_interaction_of_switching_network(As, durations,T)` gives the amount of interaction for the periodic switching network. Its inputs are the list of all adjacency matrices, $\textbf{As}$, duration of each static network in one period, and period, $T$. 
    
- `Generality_of_the_Parrondo’s_paradox.py` computes the epidemic threshold for both the static and switching networks. 
    - Function `power_method(A,tolerance)` finds the largest eigenvalue of the matrix, $\textbf{A}$ by the power method. Its inputs are the adjacency matrix, $\textbf{A}$, and the tolerance of the power method.
    - Function `get_x_intercept(function,tolerance)` computes the root of the function [NM: What do you mean by 'system'? 'system' is almost always a vague word, so I'd either avoid or specify. Do you mean 'the root of the function given as input'?] [Maisha: System means the function.I changed it as function now.]by a root finding algorithm.
    - Function `compute_epidemic_threshold(As_list)` measures the epidemic threshold for both static and switching network which takes the list of all adjacency matrices as input and returns the epidemic threshold for static networks and switching networks.
      
- `Perturbation_theory_for_the_largest_Floquet_exponent.py`  has function for finding the first order approximation of largest Floquet exponent. 
    - Function `approximate_matrix(beta,mu,T,barA)` computes the first-order approximation of matrix $\mathcal{T}$. Its inputs are the infection rate, recovery rate, period (i.e., $T$), and the time-averaged adjacency matrix (i.e., $\overline{\textbf{A}}$).
    - Function `approximate_largest_eigenvalue(beta,mu,T,barA)` computes the first-order approximation of the largest eigenvalue of $\mathcal{T}$. Its inputs are the same as those of `approximate_matrix(beta,mu,T,barA)`.
    - Function `approximate_lamda_F(beta,mu,T,barA)` computes the first-order approximation of the largest Floquet exponent, $\lambda_{\text{F}}$. Its inputs are the same as those of `approximate_matrix(beta,mu,T,barA)`. 

- `Parrondo’s_paradox_and_anti-phase_oscillation.py` computes the anti-phase oscillation.
    - Function `fraction_anti_phase(z, t_start, t_end)` measures the extent of anti-phase oscillation by the fraction of time during one cycle or period. Its inputs are:
       - $z$ : state of the SIS model at each time step,
       
