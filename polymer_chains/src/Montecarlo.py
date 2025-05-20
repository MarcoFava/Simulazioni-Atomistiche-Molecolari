import numpy as np
import random

class Montecarlo:
    def __init__(self,system):
        self.system = system
        self.positions = system.positions
        self.lattice = system.lattice
    
    def move(self):
        # Either take the first particle and add it to the tail (drop=0, insert=-1)
        # or the opposite (drop=-1, insert=0)
        drop, insert = random.choice([[0, -1], [-1, 0]])

        # Attempt to put the particle at position `pos`
        delta = random.choice([[1, 0], [-1, 0], [0, 1], [0, -1]])

        pos = self.positions[insert][:] + np.array(delta)

        # If the new site is empty, put the particle there
        if self.lattice[pos[0], pos[1]] == 0:
            new = self.positions.pop(drop)
            self.lattice[new[0], new[1]] = 0
            self.lattice[pos[0], pos[1]] = 1
            new[:] = pos

            # Either add to the head or the tail
            if insert == 0:
                self.positions.insert(0, new)
            else:
                self.positions.append(new)
    
    def run(self):
        for _ in range(np.shape(self.positions)[0]):
            self.move()