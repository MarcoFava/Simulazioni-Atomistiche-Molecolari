import matplotlib.pyplot as plt
import Montecarlo as mc
import system as sys
import numpy as np
import random
import copy
from SAW_module import setup, mc_move, simulation, plot_val_over_NMC, unfold, visualize_chain, overlap

N = 100
L = 1000
snapshots = 100
# steps = 12*N
burnin = 3*N
steps = snapshots*burnin


# job(N=N,L=L,steps=steps,burnin=burnin)
data = simulation(N=N,L=L,steps=steps,burnin=burnin)

# visualize_chain(data['trajectory'][0])
# visualize_chain(data['trajectory'][1])
# visualize_chain(data['trajectory'][2])


my_sys = sys.System(N=N,L=L)
montecarlo = mc.Montecarlo(system=my_sys)

print(len(data['trajectory']))
no_overlaps = 0
tot_overlap = 0
for i in range(len(data['trajectory'])):
    #if i%1 == 0: print(i)
    for j in range(len(data['trajectory'])):
        if i<=j : pass
        tot_overlap += 1

        my_sys.positions = data['trajectory'][i]
        cm = my_sys.centre_of_mass

        my_sys.positions = copy.deepcopy(data['trajectory'][j])
        my_sys.centre_of_mass = cm

        if overlap(data['trajectory'][i],my_sys.positions) is None: 
            no_overlaps += 1
            print(f'rate overlaps = {no_overlaps/tot_overlap}, {i}, {j}')
            # print(overlap(data['trajectory'][i],my_sys.positions))
            # break
print(no_overlaps/tot_overlap)
