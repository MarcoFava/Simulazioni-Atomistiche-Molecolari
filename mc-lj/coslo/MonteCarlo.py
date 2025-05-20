import numpy as np
import random as rnd



class MonteCarlo:
    def __init__(self, system, interaction):
       self.system = system
       self.positions = system.positions
       self.box = system.box
       self.interaction = interaction

       self.beta = 1./system.temperature

       self.n_mc_moves = 0
       self.n_accepted_moves = 0

    @property
    def accept_ratio(self):
        if self.n_mc_moves > 0 : return self.n_accepted_moves/self.n_mc_moves
        else: return 0.
    
    def run(self,mc_param=0.4):
        for _ in range(self.positions.shape[0]):
            self.move(mc_param)
        # return self.positions


    def move(self,mc_param=0.4):
        self.n_mc_moves += 1

        index = rnd.randint(0,self.positions.shape[1]-1)
        temp_pos = self.positions[:,index].copy()
        U0 = self.interaction.compute_single_energy(self.positions,self.box,index)

        self.positions[:,index] += [(rnd.random()-0.5)*mc_param for _ in range(self.positions.shape[0])]

        U1 = self.interaction.compute_single_energy(self.positions,self.box,index)
        energy_diff = U1 - U0

        if (energy_diff < 0):
            # COSLO: one should be careful with this idea of accumulating energies.
            # Keep in mind that roundoff errors will accumulate (albeit VERY slowly).
            # I suggest you do not and rather recompute the full energy when needed
            # in System. System should call Interaction.compute_energy() when requested
            
            # print(f"accepted: {energy_diff:2}")
            self.n_accepted_moves += 1
            self.system._energy += energy_diff
            # return self.positions
        else:
            p = np.exp(-self.beta*energy_diff)
            r = rnd.random()
            # COSLO: parenthesis are not needed
            # if (r < p):
            if r < p:
                # print(f"accepted: {r:.3}, {p:.2}, {energy_diff:2}")
                self.n_accepted_moves += 1
                self.system._energy += energy_diff
                # return self.positions
            else:
                self.positions[:,index] = temp_pos
                # print(f"rejected: {r:.3}, {p:.2}, {energy_diff:2}")
                # return self.positions

