#!/usr/bin/env python
# coding: utf-8

# command line arguments are as follows:
# python bay_opt.py <value_of_alpha> <num_of_files_with_weights> <directory_with_graph_weights>

# Importling the relevant libraries

import numpy as np
import matplotlib.pyplot as plt
import time
import glob
import sys
from skopt import gp_minimize
from skopt.learning import GaussianProcessRegressor
from skopt.learning.gaussian_process.kernels import ConstantKernel, Matern
import warnings
warnings.filterwarnings("ignore")

# Mapping the command line arguments to variables
alpha = float(sys.argv[1])
num_of_files = int(sys.argv[2])
directory = sys.argv[3]                                               # ensure that directory includes *.rud in the end
    
results_filepath = "results_bo.txt"    
results_file = open(results_filepath, 'w')
results_file.write("N \t Time (s) \t MAX-CUT \t Filename \n")


# Definitions of several functions

def coupling_matrix(N, weights, J1=None):                             # J1 is to be used if there is a necessity to explicitly define a coupling matrix
    if J1 is not None:                                   
        return J1
    J = np.zeros((N, N))
    k = 0
    for i in range(N-1):
        for j in range(N-1, i, -1):
            J[i,j] = weights[k]                                       
            J[j,i] = weights[k]                                       # Ensuring that J is symmetric
            k = k + 1
    return J

def voltage(alpha, beta, N, Tmax, ret='x', pi_multiplier=1):          # the variable pi_multiplier can be used to run the simulations for different biasing points of the IM
    x = np.random.normal(-0.4, 0.4, size=(N, Tmax))
    f = np.zeros((N, Tmax))                                           # the feedback function f
    #g = np.random.normal(0, 0.01, size=(N, Tmax))                     # gaussian noise to be added to each spin at each iteration
    
    J = coupling_matrix(N, weights)
    
    for k in range(Tmax-1):
        for n in range(N):
            sum1 = 0
            if beta != 0:
                for j in range(N):
                    sum1 += J[j, n]*(x[j, k])
            f[n, k] = alpha*(x[n, k]) + beta*sum1                     # equation 3.2 in the paper
            #f[n, k] *= 0.9                                           # since f = V/V_pi, f can be scaled to see the impact of different voltages that can be applied and if an amplifier would be needed
            x[n, k+1] = ((np.cos(f[n, k] - pi_multiplier*np.pi/4))**2 - 0.5 )    # equation 3.1 in the paper
    
    if ret=='x':
        return x
    if ret=='f':                       
        return f                                                      # Returning f if we want to observe the behaviour of f for different parameters

def energy(Tmax, N, alpha, beta, array=0):
    x = voltage(alpha, beta, N, Tmax)
    E = np.zeros(Tmax)
    J = coupling_matrix(N, weights=weights)
    for k in range(Tmax):
        for i in range(N):
            for j in range(i, N):
                E[k] += J[i, j]*np.sign(x[i, k])*np.sign(x[j, k])
    if array == 0:
        return -E[-1]
    if array == 1:                                                   #array is needed if visualization of convergence is needed
        return -E

    
# Below functions plots the surrogate model and the sampling points used to train the model
def plot_approximation(gpr, X, X_sample, Y_sample, show_legend=False):
    mu, std = gpr.predict(X, return_std=True)
    plt.fill_between(X.ravel(), 
                     mu.ravel() + 1.96 * std, 
                     mu.ravel() - 1.96 * std, 
                     alpha=0.1) 
    plt.plot(X, mu, 'b-', lw=1, label='Surrogate function')
    plt.plot(X_sample, Y_sample, 'kx', mew=3, label='samples')
    
    if show_legend:
        plt.legend()
        
def initialize_weights_variables(filename, density=50, min_weight=1, max_weight=10):
    global weights
    weights = []

    if density != 50:
        for i in range(len(weights)):
            if np.random.randint(0,100) >= density:
                weights[i] = 0
            else:
                weights[i] = -np.random.randint(min_weights, max_weights)
    else:

        with open(filename, 'r') as file:
            for line in file:
                words = line.split()
                if len(words) == 3:
                    weight = int(words[2])
                    weights.append((weight-1)/2)

        file.close()

    print("density is ", abs(100*sum(weights)/len(weights)), '%')

    E_init = np.zeros(2).reshape(-1, 1)

    p = [0.5, -0.5, -len(weights)]
    roots = np.roots(p)

    N = int(np.round(roots[roots>0][0]))           # Finding the positive root 

    return weights, N, E_init


# Bayesian Optimization function to find the optimal values of $\alpha$ and $\beta$ that achieve the ground state configuration.

def bay_opt(alpha, optimal_energy, N, Tmax, E_init, noise=0.005, beta_max=0, num=10, plot_graph=1,n_calls=4,
           length_scale=0.1):
    start = time.time()
    beta_min = N*(1-alpha)/(abs(sum(weights)))

    if beta_max == 0:
        beta_max = num*beta_min
    else:
        beta_max = beta_max

    bounds = np.array([[beta_min, beta_max]])

    X_init = np.array([[beta_min], [beta_max]])
    E_init[0] = energy(Tmax, N, alpha, X_init[0])
    E_init[1] = energy(Tmax, N, alpha, X_init[1])

    # In case n_calls needs an adaptive definition
    # n_calls = int((9*N*(1-alpha)/(abs(sum(weights))))/0.003) - 2
    n_calls = n_calls

    n52 = ConstantKernel(1.0) * Matern(length_scale=length_scale, length_scale_bounds=(1e-10, 10), nu=2.5)
    gpr = GaussianProcessRegressor(kernel=n52, alpha=noise)

    r = gp_minimize(lambda x: energy(Tmax, N, alpha, np.array(x)), 
                    bounds.tolist(),
                    base_estimator=gpr, 
                    acq_func='EI',      # expected improvement
                    xi=0.1,             # exploitation-exploration trade-off
                    n_calls=n_calls,    # number of iterations
                    n_random_starts=0,  # initial samples are provided thus n_random_starts=0
                    x0=X_init.tolist(), # initial samples
                    y0=E_init.ravel())

    gpr.fit(r.x_iters, r.func_vals)

    optimal_energy[0].append(r.x[0])
    optimal_energy[1].append(r.fun)
    optimal_energy[2].append(alpha)

    #The first row stores the value of beta
    #The second row stores the value of the accuracy
    #The third row stores the value of alpha

    # Plot the fitted model and the noisy samples
    if plot_graph == 1:
        X = np.linspace(bounds[:, 0], bounds[:, 1], 1000).reshape(-1, 1) #beta values for plotting
        plt.title("alpha = {}".format(alpha))
        plt.xlabel("beta")
        plt.ylabel("Energy")
        plot_approximation(gpr, X, r.x_iters, r.func_vals, show_legend=True)
        plt.show()

    end = time.time()
    
    answers = [optimal_energy, r, end-start, n_calls, beta_min, beta_max, end-start]
    
    return answers

if num_of_files > 1:

    # Bayesian Optimization done for multiple files at once

    file_paths = glob.glob(directory)

    answers_multiple = []
    optimal_energy_multiple = []
    alpha_g_multiple = []
    beta_g_multiple = []
    i = 0
    noise = 0.01

    for i in range(num_of_files):
        
        file_name = file_paths[i]
        weights, N, E_init = initialize_weights_variables(file_name)
        Tmax = 50
        print("N = ", N)
        answers_multiple.append(bay_opt(alpha=alpha, optimal_energy=[[], [], []], N=N, Tmax=Tmax, E_init=E_init, 
                             noise=noise, num=8, n_calls=2, plot_graph=0))

        optimal_energy_multiple.append(answers_multiple[-1][0])
        beta_g_multiple.append(optimal_energy_multiple[-1][0][0])
        alpha_g_multiple.append(optimal_energy_multiple[-1][2][0])
        C = 0.5*(sum(weights)+optimal_energy_multiple[-1][1][0])

        print("Time taken = ", answers_multiple[-1][2])
        print("Cut value = ", 0.5*(sum(weights)+optimal_energy_multiple[-1][1][0]))
        print(optimal_energy_multiple[-1],"\n\n")
        
        results_file.write("{} \t {} \t {} \t {} \n".format(N, np.round(answers_multiple[-1][2],2), abs(C), file_name))

elif num_of_files == 1:

    # One file at a time

    start = time.time()
    filename = directory
    weights, N, E_init = initialize_weights_variables(filename, density=50)
    Tmax = 50
    end = time.time()
    print("time taken for defining the system = ", end-start, "s")


    answers = bay_opt(alpha=alpha, optimal_energy=[[], [], []], N=N, 
                        Tmax=Tmax, E_init=E_init, 
                        #beta_max=5*N*(1-alpha)/(abs(sum(weights))),
                        num=8, n_calls=1, plot_graph=0)


    optimal_energy = answers[0]
    print("Optimal value is obtained when alpha = {}, beta = {}".format(
        optimal_energy[2][0], optimal_energy[0][0]))
    print("Time taken to find the optimal gains is {}s".format(answers[6]))
    C = 0.5*(sum(weights)+optimal_energy[1][0])
    print("\n\n Max cut value after ", Tmax ," iterations = ", C)

    print("Do you want to plot the convergence of the energy for these optimal values for 150 iterations to obtain a potentially better value of MAX CUT? [y/n]")
    choice = input("\n")
    
    results_file.write("{} \t {} \t {} \t {} \n".format(N, np.round(answers[6],2), abs(C), filename))
    
    if choice == 'y':

        print("\n\n Plotting the convergence plot for the optimal values of alpha and beta for 150 iterations \n\n ")

        alpha_g = answers[0][2][0]
        beta_g = answers[0][0][0]

        y = energy(150, N, alpha_g, beta_g, 1)
        plt.title("Convergence of energy. N = {}, Density = {}".format(N, 50))
        plt.plot(y)
        plt.ylabel("Energy")
        plt.xlabel("Number of iterations")
        plt.show()
        print("alpha={}, beta={}".format(alpha_g, beta_g))


        C = 0.5*(sum(weights)+optimal_energy[1][0])
        print("\n\n Max cut value after ", Tmax ," iterations = ", C)
        C = 0.5*(min(y)+sum(weights))
        print("\n\n Max cut value after 150 iterations = ", C)
