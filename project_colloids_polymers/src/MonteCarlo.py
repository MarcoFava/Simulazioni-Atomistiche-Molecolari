# import numpy as np
import random as rnd
from numba import njit
from math import exp



class MonteCarlo:
    def __init__(self, system, interaction, seed=None):
        self.system = system
        self.positions = system.positions
        self.box = system.box
        self.interaction = interaction
        self.n_coll = system.n_coll
        self.n_part = system.n_part

        self.seed = seed
        if seed is not None: rnd.seed(seed)

        # self.beta = 1./system.temperature

        self.tot_mc_moves = 0
        self.tot_coll_moves = 0
        self.accepted_moves = 0
        self.accepted_moves_coll = 0

        self.mc_param = 0.4 * system.box[0]
        self.mc_param_coll = 2*system.r_part
        self.mc_param_part = system.r_coll


    @property
    def accept_ratio(self):
        if self.tot_mc_moves > 0 : return self.accepted_moves/self.tot_mc_moves
        else: return 0.
    
    @property
    def accept_ratio_coll(self):
        # if self.tot_mc_moves > 0: return self.accepted_moves_coll/(self.tot_mc_moves /(self.n_coll+self.n_part)*self.n_coll)
        if self.tot_coll_moves > 0: return self.accepted_moves_coll/self.tot_coll_moves
        else: return 0.

    def run(self):
        self.mc_param = 0.4 * self.box[0]
        self.tot_mc_moves += self.system.current_N
        for _ in range(self.system.current_N):
            self.displacement_move_hs()
        # if self.interaction.check_overlap_all(self.positions,self.box): 
        #     raise KeyError('Error: Overlap found')


    def displacement_move_hs(self):
        """
        Select a particle i at random
        Propose a random displacement dr of particle i (ex. over a cube of side dl)
        Check for any overlaps using the hard spheres potential
        If the move is rejected because of any overlap, keep the old configuration
        """
        pos = self.positions
        index = rnd.randint(0,self.system.current_N-1)
        if index < self.n_coll: 
            self.tot_coll_moves += 1
            # mc_param = 0.1*self.box[0]
            mc_param = self.mc_param_coll
        else:
            # mc_param = 0.4*self.box[0]
            mc_param = self.mc_param_part
        temp_pos = pos[:,index].copy()

        # pos[:,index] += [(rnd.random()-0.5)*self.mc_param for _ in range(pos.shape[0])]
        # mc_param=0.1
        displace_particle(pos,index,mc_param)
        self.interaction.apply_pbc(pos[:,index],self.box)

        if not self.interaction.check_single_overlap(pos,self.box,self.n_coll,index):
            if index < self.n_coll: self.accepted_moves_coll += 1
            self.accepted_moves += 1
            # print('accepted')
        else:
            pos[:,index] = temp_pos
            # print('rejected')


    def displacement_move(self):
        """
        Select a particle i at random
        Propose a random displacement dr of particle i (ex. over a cube of side dl)
        Calculate the energy change dU between the new and old configurations
        Accept the proposed move with probability min{1,exp(-beta*dU)}. 
        If the move is rejected, keep the old configuration
        """
        pos = self.positions[:,:self.system.current_N]
        index = rnd.randint(0,self.system.current_N-1)
        temp_pos = pos[:,index].copy()

        U0 = self.interaction.compute_single_energy(pos,self.box,index)

        # pos[:,index] += [(rnd.random()-0.5)*self.mc_param for _ in range(pos.shape[0])]
        displace_particle(pos,index,self.mc_param)
        self.interaction.apply_pbc(pos[:,index],self.box)

        U1 = self.interaction.compute_single_energy(pos,self.box,index)
        energy_diff = U1 - U0

        if (energy_diff < 0):
            # print(f"accepted: {energy_diff:2}")
            self.accepted_moves += 1
            self.system._energy += energy_diff
        else:
            p = exponential(-self.beta*energy_diff)
            # p = np.exp(-self.beta*energy_diff)
            r = rnd.random()
            if r < p:
                # print(f"accepted: {r:.3}, {p:.2}, {energy_diff:2}")
                self.accepted_moves += 1
                self.system._energy += energy_diff
            else:
                pos[:,index] = temp_pos
                # print(f"rejected: {r:.3}, {p:.2}, {energy_diff:2}")



@njit
def displace_particle(pos, index, mc_param):
    # dr = (np.random.rand(pos.shape[0]) - 0.5) * mc_param
    # pos[:, index] += dr
    pos[0,index] += (rnd.random()-0.5)*mc_param
    pos[1,index] += (rnd.random()-0.5)*mc_param
    pos[2,index] += (rnd.random()-0.5)*mc_param

@njit
def exponential(x):
    return exp(x)
    # return np.exp(x)
