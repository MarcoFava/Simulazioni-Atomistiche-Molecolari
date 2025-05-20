import numpy as np
import copy
from f2py_jit import jit

potential = 'lj'
potential_file = f'src/potential_module_{potential}.f90'
interaction_file = 'src/interaction_module.f90'
files = [potential_file, interaction_file]

f90 = jit(files, inline=True, flags='-fbounds-check')

N = 10
positions = np.ones((3,N),dtype=float,order='F')
box = np.ones(3)*0.75

# index = 0
# for index in range(positions.shape[1]):
#     pos = copy.deepcopy(positions[:,index])
#     f90.interaction_module.pbc(pos,box)
#     positions[:,index] = pos

def pbc(positions,box):
    f90.interaction_module.full_pbc(positions,box)

print(positions)
pbc(positions,box)

print(positions)

