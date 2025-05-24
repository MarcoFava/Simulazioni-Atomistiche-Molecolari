import matplotlib.pyplot as plt
import Montecarlo as mc
import system as sys
import numpy as np
import random
import copy
import pantarei as rei
from SAW_module import setup, mc_move, simulation, plot_val_over_NMC, unfold, visualize_chain, overlap, prob_overlap

#scheduler = rei.Scheduler.ThreadScheduler(backend='process')
task = rei.Task(simulation) #,clear_first=True)
#job = rei.Thread(task, scheduler=scheduler)

N = 100
L = 1000
snapshots = 50
snapshots = 10
# steps = 12*N
burnin = 3*N
steps = snapshots*burnin



data = task(N=N,L=L,steps=steps,burnin=burnin)
#scheduler.wait()
#data = simulation(N=N,L=L,steps=steps,burnin=burnin)
# visualize_chain(data['trajectory'][2])

len_traj = len(data['trajectory'])
for i in range(len_traj):
    my_sys.positions = copy.deepcopy(data['trajectory'][i])
    #        for deg in [90, 180, 270]:
    for _ in range(3):
        data['trajectory'].append(my_sys.rotate_pos(deg=90))
print(f"Total number of configurations is: {len(data['trajectory'])}")


bins = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
bins = [0, 4, 8, 12, 16, 20]
drs = [bins[i]+0.5 for i in range(len(bins)-1)]
counts = []
for dr in drs:
    counts.append(prob_overlap(dr=dr,data=data))

plt.stairs(counts,bins,fill=True)
plt.savefig('prob_hist.png')


#my_sys = sys.System(N=N,L=L)
#montecarlo = mc.Montecarlo(system=my_sys)

#print(len(data['trajectory']))
#no_overlaps = 0
#tot_overlap = 0
#for i in range(len(data['trajectory'])):
#   #if i%1 == 0: print(i)
#    for j in range(len(data['trajectory'])):
#        if i<=j : pass
#        tot_overlap += 1
#
#        my_sys.positions = data['trajectory'][i]
#        cm = my_sys.centre_of_mass
#
#        my_sys.positions = copy.deepcopy(data['trajectory'][j])
#        my_sys.centre_of_mass = cm + dr
#
#        if overlap(data['trajectory'][i],my_sys.positions) is None: 
#            no_overlaps += 1
#            print(f'rate overlaps = {no_overlaps/tot_overlap:.3}, {i}, {j}')
#            # print(overlap(data['trajectory'][i],my_sys.positions))
#            # break
#print(no_overlaps/tot_overlap)
