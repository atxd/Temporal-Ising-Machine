#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

density=float(50)

bo_filename = "results_bo_density_" + str(density) + ".txt"
gs_filename = "results_gs_density_" + str(density) + ".txt"

bo_nodes = []
bo_times = []
bo_cuts = []

with open(bo_filename, 'r') as bo_file:
    for line in bo_file:
        words = line.split()
        if words[0] !='N':
            bo_nodes.append(int(words[0]))
            bo_times.append(float(words[1]))
            bo_cuts.append(float(words[2]))

gs_nodes = []
gs_times = []
gs_cuts = []
with open(gs_filename, 'r') as gs_file:
    for line in gs_file:
        words = line.split()
        if words[0] != 'N':
            gs_nodes.append(int(words[0]))
            gs_times.append(float(words[1]))
            gs_cuts.append(float(words[2]))

print("\nPlotting the times of Bayesian Optimization and Grid Search approaches\n")


plt.title("Comparison of times")
plt.xlabel("Number of Nodes")
plt.ylabel("Times (s)")
plt.plot(bo_nodes, bo_times, marker='x', ls='', color='red', label='Bayesian Optimization')
plt.plot(gs_nodes, gs_times, marker='x', ls='', color='blue', label='Grid Search')
plt.legend()
plt.show()


# Solution Accuracy Comparison

bo_times.sort()
gs_times.sort()
gs_nodes.sort()
bo_nodes.sort()
gs_cuts.sort()
bo_cuts.sort()


#Checks if there is any mismatch in the nodes of the two files. If there is a mismatch it means that the same files have not been used for BO and GS

i_mismatch = [None]
for i in range(len(gs_nodes)):
    if bo_nodes[i] - gs_nodes[i] != 0:
        i_mismatch.append(i)
        print("Nodes at index", i+1, "are different")


time_ratios = []
cut_ratios = []

i = 0
while(i < len(gs_nodes) and i not in i_mismatch):
    time_ratios.append(bo_times[i]/gs_times[i])
    cut_ratios.append(bo_cuts[i]/gs_cuts[i])
    i += 1



comparison_filename = "comparison_density_" + str(density) + ".txt"
comparison_file = open(comparison_filename, 'w')

comparison_file.write("Nodes \t\t Time \t\t CutValue \n")

print("\nTimes and Cut values of BO relative to GS\n")

print("Nodes\t\t Time \t\t Cutvalue")
for i in range(len(time_ratios)):
    print(bo_nodes[i], '\t\t', np.round(time_ratios[i],3), '\t\t', np.round(cut_ratios[i],3))
    comparison_file.write("{} \t\t {} \t\t {} \n".format(bo_nodes[i], np.round(time_ratios[i],3),
                                                         np.round(cut_ratios[i],3)))


# Viewing the actual values if needed
'''
print("\nResults of both approaches listed separately\n")

print("\tBay Opt \t\t\t\t Grid Srch")
print("Nodes\t\tCut value\t\tNodes\t\tCut value")
for i in range(len(gs_nodes)):
    print(bo_nodes[i],'\t\t',bo_cuts[i],'\t\t',gs_nodes[i],'\t\t',gs_cuts[i])
'''
