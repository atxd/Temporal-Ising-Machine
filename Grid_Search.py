# Command line arguments are in the following syntax
# python Grid_Search.py <alpha_value> <directory> <num_of_files>

import numpy as np
import sys
import time
import glob


def coupling_matrix(N, weights, J1=None):
    if J1 is not None:
        return J1
    J = np.zeros((N, N))
    k = 0
    for i in range(N-1):
        for j in range(N-1, i, -1):
            J[i,j] = weights[k]
            k = k + 1
    J = J + np.transpose(J)
    # ensuring that J is symmetric
    
    return J

def voltage(alpha, beta, N, Tmax, ret='x', pi=1):  
    x = np.random.normal(-0.4, 0.4, size=(N, Tmax))
    f = np.zeros((N, Tmax))   #the function f
    #g = np.random.normal(0, 0.01, size=(N, Tmax)) #gaussian noise
    
    J = coupling_matrix(N, weights)

    for k in range(Tmax-1):
        for n in range(N):
            sum1 = 0
            if beta != 0:
                for j in range(N):
                    sum1 += J[j, n]*(x[j, k])
            f[n, k] = alpha*(x[n, k]) + beta*sum1                              #equation 3.2 in the paper with beta=0
            #f[n, k] *= 0.9
            x[n, k+1] = ((np.cos(f[n, k] - pi*np.pi/4))**2 - 0.5 )                #equation 3.1 in the paper
    
    if ret=='x':
        return x
    if ret=='f':
        return f

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
    if array == 1:
        return -E
        
def assign_weights(filename, density=50):
    global weights
    weights = []

    with open(filename, 'r') as file:
        for line in file:
            words = line.split()
            if len(words) == 3:
                weight = int(words[2])
                weights.append(weight)

    file.close()

    if density != 50:
        for i in range(len(weights)):
            if np.random.randint(0,100) >= density:
                weights[i] = 0
            else:
                weights[i] = -abs(weights[i])

    for i in range(len(weights)):
        if weights[i] == 1:                         # weights now belong to {-1, 0}
            weights[i] = 0

    print("density is ", abs(100*sum(weights)/len(weights)), '%')

    p = [0.5, -0.5, -len(weights)]
    roots = np.roots(p)

    N = int(np.round(roots[roots>0][0]))
    
    return weights, N


# Mapping the command line arguments to the respective variables
for i in range(0, len(sys.argv)-1):
    if sys.argv[i]=='ALPHA':
        alpha = float(sys.argv[i+1])
    elif sys.argv[i]=='NUM_OF_FILES':
        num_of_files = int(sys.argv[i+1])
    elif sys.argv[i]=='DIRECTORY':
        directory = sys.argv[i+1]


start = time.time()
file_paths = glob.glob(directory)
results_filepath = "results_gs.txt"
results_file = open(results_filepath, 'w')
results_file.write("N \t Time (s) \t MAX-CUT \t Filename \n")

for i in range(num_of_files):
    filename = file_paths[i]
    weights, N = assign_weights(filename, density=50)
    Tmax = 50
    end = time.time()
    print("N = {}".format(N))
    print("time taken for defining the system = ", end-start, "s")

    values = []
    beta_min = N*(1-alpha)/(abs(sum(weights)))
    beta_max = 9*beta_min
    beta_space = np.linspace(beta_min, beta_max, 5)

    start = time.time()
    for j in beta_space:
        values.append(energy(Tmax, N, alpha, j))
    end = time.time()
    print("Time: {}s".format(end-start))

    C = 0.5*(min(values) + sum(weights))
    print("The best value of Cut is ", C)

    results_file.write("{} \t {} \t {} \t {} \n".format(N, np.round(end-start,2), abs(C), filename))
    
results_file.close()