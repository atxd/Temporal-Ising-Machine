# Employing Bayesian Optimization in Temporal Ising Machine
An Ising spin system is simulated by creating a set of artifical spins that are characterized in a temporal basis i.e. each spin corresponds to a photovoltage for a specified time range. An intensity modulator (MZM) is used to determine the spins. The transfer functions are given by,

```math
x_i[k+1] = \text{cos}^2\left(f_i[k] + \zeta_i[k] - \frac{\pi}{4}\right) - \frac{1}{2}
```

```math

f_i[k] = \alpha x_i[k] + \beta \sum_j J_{i,j} x_j[k]
```

where $x_i[k]$ denotes the photovoltage of the $i$-th spin at the $k$-th instant of time. $\zeta_i[k]$ corresponds to gaussian noise in the system. $f_i[k]$ is the feedback RF bias provided to the intensity modulator. Coefficients $\alpha$ and $\beta$ are the linear gain and coupling gain respectively. The $i$-th spin is now given by $\sigma_i$ = sgn($x_i$).

Clearly, the evolution of this spin system depends on the gain parameters $\alpha$ and $\beta$. The Ising Hamiltonian can then be written as,
```math \label{eq: Hamiltonian}
H = \sum_{i, j} J_{i,j} \sigma_i \sigma_j
```

It can be shown that the minimum energy (the ground state energy) of this Hamiltonian can be mapped to the max-cut value of a graph with $N$ nodes, where $N$ is the same as the number of spins in the Ising Hamiltonian, and the edges connecting the nodes in a graph corresponds to the interaction energies of 2 spins in the Ising Hamiltonina.

The choice of the gain parameters $\alpha$ and $\beta$ is crucial since it determines the final stable value of each spin (if there is a stable configuration). Previously, a brute force approach (Grid Search algorithm) was used to find the optimal values of the gain parameters $\alpha$ and $\beta$ that ensure that the system evolves into its ground state configuration. However, in this project, an attempt has been made to employ Bayesian Optimization to search for the optimal values of the gain parameters. The improvements in solution accuracy and in the time taken to find the optimal solution was observed in many of the cases.

# Using the Graph Generator (rudy)
Directions to use the graph generator are given in the README.md file in the Graph Generator directory. Additionally, some graphs have already been generated and are present in the same directory.

# Using Bay_opt.py
Run the following commands on the terminal and specify the actual values of the arguments wherever necessary.
```
python Bay_opt.py ALPHA <value_of_alpha> NUM_OF_FILES <num_of_files> DIRECTORY <directory_with_graph_weights> DENSITY <density_of_edges>
```
There is no constraint on the order of the arguments, only that the parameter and its value be adjacent entries. Here is an example of a command that one could give
```
python Bay_opt.py ALPHA 0.9 NUM_OF_FILES 8 DENSITY 50 DIRECTORY "Graph Generator\\*.rud"
```
A new file results_bo_density_<density_of_edges>.txt will be created which will save the results (Number of Nodes (N), Time taken to find the solution (T), MAX-CUT value and filename. By default, the weights used will have 50% density and all the weights are either 0 or -1. The choice of nodes whose edges are connected is random.

# Using Grid_Search.py
Run the following commands on the terminal and specify the actual values of the arguments wherever necessary.
```
python Grid_Search.py ALPHA <alpha_value> DIRECTORY <directory> NUM_OF_FILES <num_of_files> DENSITY <density_of_edges>
```
There is no constraint on the order of the arguments, only that the parameter and its value be adjacent entries. Here is an example of a command that one could give
```
python Grid_Search.py DENSITY 50 ALPHA 0.9 DIRECTORY "Graph Generator\\*.rud" NUM_OF_FILES 16
```
A new file results_gs_density_<density_of_edges>.txt will be created which will save the results (Number of Nodes (N), Time taken to find the solution (T), MAX-CUT value and filename. By default, the weights used will have 50% density and all the weights are either 0 or -1. The choice of nodes whose edges are connected is random.

# Interpreting the results files
Both the results files consists of data that specifies the size of the graph, the time taken by the simulation to calculate the cut value and the file that contained the weights. The results are shown in a tabular form. Additionally, the "Comparison_Results.py" file can be run to get one table that provides a ratio of the times and cut values obtained using Bayesian Optimization wrt Grid Search. There are no command line arguments required to run this Python file. The comparison results are stored in a text file named "comparison_density_<density_of_edges>.txt"
