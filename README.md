# Emplying Bayesian Optimization in Temporal Ising Machine
An Ising spin system is simulated by creating a set of artifical spins that are characterized in a temporal basis i.e. each spin corresponds to a photovoltage for a specified time range. An intensity modulator (MZM) is used to determine the spins. The transfer functions are given by,
```math \label{eq: transfer equation}
x_i[k+1] = \text{cos}^2\left(f_i[k] + \zeta_i[k] - \frac{\pi}{4}\right) - \frac{1}{2}
\label{eq: feedback term}
f_i[k] = \alpha x_i[k] + \beta \sum_j J_{i,j} x_j[k]
```

where $x_i[k]$ denotes the photovoltage of the $i$-th spin at the $k$-th instant of time. $\zeta_i[k]$ corresponds to gaussian noise in the system. $f_i[k]$ is the feedback RF bias provided to the intensity modulator. Coefficients $\alpha$ and $\beta$ are the linear gain and coupling gain respectively. The $i$-th spin is now given by $\sigma_i$ = sgn($x_i$).

Clearly, the evolution of this spin system depends on the gain parameters $\alpha$ and $\beta$. The Ising Hamiltonian can then be written as,
```math \label{eq: Hamiltonian}
H = \sum_{i, j} J_{i,j} \sigma_i \sigma_j
```

It can be shown that the minimum energy (the ground state energy) of this Hamiltonian can be mapped to the max-cut value of a graph with $N$ nodes, where $N$ is the same as the number of spins in the Ising Hamiltonian, and the edges connecting the nodes in a graph corresponds to the interaction energies of 2 spins in the Ising Hamiltonina.

The choice of the gain parameters $\alpha$ and $\beta$ is crucial since it determines the final stable value of each spin (if there is a stable configuration). Previously, a brute force approach (Grid Search algorithm) was used to find the optimal values of the gain parameters $\alpha$ and $\beta$ that ensure that the system evolves into its ground state configuration. However, in this project, we have tried to employ Bayesian Optimization to search for the optimal values of the gain parameters. The corresponding improvements in the search for the optimal state have also been noticed.
