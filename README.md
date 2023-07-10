# Emplying Bayesian Optimization in Temporal Ising Machine
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

The choice of the gain parameters $\alpha$ and $\beta$ is crucial since it determines the final stable value of each spin (if there is a stable configuration). Previously, a brute force approach (Grid Search algorithm) was used to find the optimal values of the gain parameters $\alpha$ and $\beta$ that ensure that the system evolves into its ground state configuration. However, in this project, an attempt has been made to employ Bayesian Optimization to search for the optimal values of the gain parameters. The corresponding improvements in the search for the optimal state have also been noticed.

# Using the Graph Generator (rudy)
Directions to use the graph generator are given in the README.md file in the Graph Generator directory.

# Using Bay_opt.py
Run the following commands on the terminal and specify the actual values of the arguments wherever necessary.
```
python Bay_opt.py <value_of_alpha> <num_of_files> <directory_with_graph_weights>
```
Here is an example of a command that one could give
```
python Bay_opt.py 0.9 8 "Graph Generator\\*.rud"
```
A new file results_bo.txt will be created which will save the results (Number of Nodes (N), Time taken to find the solution (T), MAX-CUT value and filename. By default, the weights used will have 50% density and all the weights are either 0 or -1. The choice of nodes whose edges are connected is random.

# Using Grid_Search.py
Run the following commands on the terminal and specify the actual values of the arguments wherever necessary.
```
python Grid_Search.py <alpha_value> <directory> <num_of_files>
```
Here is an example of a command that one could give
```
python Grid_Search.py 0.9 "Graph Generator\\*.rud" 16
```
A new file results_gs.txt will be created which will save the results (Number of Nodes (N), Time taken to find the solution (T), MAX-CUT value and filename. By default, the weights used will have 50% density and all the weights are either 0 or -1. The choice of nodes whose edges are connected is random.

# Interpreting the results files
Both the results files consists of data that specifies the size of the graph, the time taken by the simulation to calculate the cut value and the file that contained the weights. The results are shown in a tabular form.
