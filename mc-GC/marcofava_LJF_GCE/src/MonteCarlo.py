import numpy as np
import random as rnd
from numba import njit
from math import exp



class MonteCarlo:
    def __init__(self, system, interaction, eta=None, seed=None):
        self.system = system
        self.positions = system.positions
        self.box = system.box
        self.interaction = interaction

        self.seed = seed
        if seed is not None: rnd.seed(seed)

        self.beta = 1./system.temperature
        self.z = system.z
        # self.chem_potential = system.chem_potential
        # self.Lambda = 1./np.sqrt(3*system.temperature)
        # self.reduced_Lambda_3 = self.Lambda**3/self.system.volume
        # self.reduced_Lambda_3_inv = 1./self.reduced_Lambda_3
        
        if eta is not None:
            self.weight = np.array([np.exp(-self.beta*val)  for val in eta])
        else: self.weight = np.ones(1000)

        self.tot_mc_moves = 0
        self.accepted_moves_d = 0
        self.accepted_moves_i = 0
        self.accepted_moves_r = 0

        self.mc_param = 0.4 * system.box[0]

        self.weight_d = 1
        self.weight_i = 1
        self.weight_r = 1


    @property
    def accept_ratio(self):
        if self.tot_mc_moves > 0 : 
            return (self.accepted_moves_d + self.accepted_moves_i + self.accepted_moves_r)/self.tot_mc_moves
        else: return 0.

    def run(self):
        self.mc_param = 0.4 * self.box[0]
        for _ in range(self.system.current_N):
            self.step()

    def step(self):
        self.tot_mc_moves += 1
        # choice = rnd.choice(['d','i','r'])
        # if choice == 'd':
        #     self.displacement_move()
        # elif choice == 'i':
        #     self.insertion_move()
        # elif choice == 'r':
        #     self.removal_move()

        W = self.weight_d + self.weight_i + self.weight_r
        W_d = self.weight_d
        W_i = W_d + self.weight_i
        r = rnd.random()*W
        if r < W_d:
            self.displacement_move()
        elif r< W_i:
            self.insertion_move()
        else:
            self.removal_move()




    def displacement_move(self):
        """
        Move a random particle and accept the move if 
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
            self.accepted_moves_d += 1
            self.system._energy += energy_diff
        else:
            p = exponential(-self.beta*energy_diff)
            # p = np.exp(-self.beta*energy_diff)
            r = rnd.random()
            if r < p:
                # print(f"accepted: {r:.3}, {p:.2}, {energy_diff:2}")
                self.accepted_moves_d += 1
                self.system._energy += energy_diff
            else:
                pos[:,index] = temp_pos
                # print(f"rejected: {r:.3}, {p:.2}, {energy_diff:2}")


    def insertion_move(self):
        """
        Place a new particle at a random position in the cell
        Calculate the energy change dU between the new and old configurations
        Accept the proposed move with probability min{1,V/Lambda**3/(N+1)*exp(-beta*(dU-mu)) = exp(-beta*dU)*z/(N+1)}. 
        If the move is rejected, keep the old configuration
        """
        _N = self.system.current_N
        _new_N = _N + 1
        pos = self.positions[:,:_new_N]
        # new_pos = np.array([(2*rnd.random()-1.)*self.box[k] for k in range(pos.shape[0])]).T
        new_pos = generate_new_pos(self.box)
        
        index = _new_N - 1
        pos[:,index] = new_pos

        dU = self.interaction.compute_single_energy(pos,self.box,index)
        # energy_diff = dU - self.chem_potential

        # p = self.reduced_Lambda_3_inv / _new_N * np.exp(-self.beta*energy_diff)
        # p = np.exp(-self.beta*dU) * self.z /_new_N * self.system.volume
        _prefactor = self.z * self.system.volume
        p = exponential(-self.beta*dU)/_new_N * _prefactor * self.weight[_new_N]/self.weight[_N]
        r = rnd.random()
        if r < p:
            self.accepted_moves_i += 1
            self.system.current_N = _new_N
            self.system._energy += dU
            # print(f"accepted: {r:5.3}, {p:5.2}, {self.eta[_N]/self.eta[_new_N]:5.2}, {_new_N} --A")
        # else:
            # print(f"rejected: {r:5.3}, {p:5.2}, {self.eta[_N]/self.eta[_new_N]:5.2}, {_new_N} R")



    def removal_move(self):
        """
        Remove an existing particle chosen at random
        Calculate the energy change dU between the new and old configurations
        Accept the proposed move with probability min{1,Lambda**3*N/V*exp(-beta*(dU+mu)) = N*exp(-beta*dU)/z}. 
        If the move is rejected, keep the old configuration
        """
        _N = self.system.current_N
        _new_N = _N-1
        pos_0 = self.positions[:,:_N].copy()

        index = rnd.randint(0,_new_N-1)
        
        dU = -self.interaction.compute_single_energy(pos_0,self.box,index)
        # energy_diff = dU + self.chem_potential

        # p = self.reduced_Lambda_3 * _N * np.exp(-self.beta*energy_diff)
        _prefactor_inv = 1. / self.z / self.system.volume
        p = exponential(-self.beta*dU) * _N * _prefactor_inv * self.weight[_new_N]/self.weight[_N]

        r = rnd.random()
        if r < p:
            self.accepted_moves_r += 1
            self.system._energy += dU
            self.system.current_N = _new_N
            self.positions[:,index] = pos_0[:,-1]
            
        #     print(f"accepted: {r:5.3}, {p:5.2}, {self.eta[_N]/self.eta[_new_N]:5.2}, {_new_N} --A")
        # else:
        #     print(f"rejected: {r:5.3}, {p:5.2}, {self.eta[_N]/self.eta[_new_N]:5.2}, {_new_N} R")


@njit
def displace_particle(pos, index, mc_param):
    # dr = (np.random.rand(pos.shape[0]) - 0.5) * mc_param
    # pos[:, index] += dr
    pos[0,index] += (rnd.random()-0.5)*mc_param
    pos[1,index] += (rnd.random()-0.5)*mc_param
    pos[2,index] += (rnd.random()-0.5)*mc_param

@njit
def generate_new_pos(box):
    dim = box.shape[0]
    new_pos = np.empty(dim)
    for k in range(dim):
        new_pos[k] = (2.0 * np.random.rand() - 1.0) * box[k]
    return new_pos
    
@njit
def exponential(x):
    return exp(x)
    # return np.exp(x)
