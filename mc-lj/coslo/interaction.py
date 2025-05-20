import sys
from f2py_jit import jit

class Interaction:
    # COSLO: I suggest to rename which_potential -> potential?
    def __init__(self, parameters, potential=None, files=None,
                 flags='-O3 -ffast-math'):

        self.which_potential = which_potential

        # COSLO: use standard way to check this 
        # if (files!=None):
        if files is not None:
            # COSLO: there was a typo here!
            self.files = files
        else:
            # COSLO: standard way to interpolate strings is with f-strings
            # potential_file = 'src/potential_module_'+which_potential+'.f90'
            potential_file = f'src/potential_module_{which_potential}.f90'
            interaction_file = 'src/interaction_module.f90'
            self.files = [potential_file, interaction_file]
            try:
                open(potential_file)
                open(interaction_file)
            except:
                print('ERROR: Unable to open the files ' + potential_file + ' and/or ' + interaction_file)

        self.flags = flags
        
        self.f90 = jit(self.files, inline=True,
                flags=self.flags)

        # initialize potential
        self.f90.potential_module.initialize(**parameters)
        
    def compute_single_energy(self, positions, box, index):
        return self.f90.interaction_module.compute_single_interaction(positions,box,index+1)

    # COSLO: rename it compute_energy(), because total_energy suggests E_pot+E_kin
    def compute_total_energy(self, positions, box):
        return self.f90.interaction_module.compute_interaction(positions,box)

    def compute_virial(self, positions, box):
        return self.f90.interaction_module.compute_virial(positions,box)
    
    def compute_pressure(self, positions, box, density, kB_T, volume):
        W = self.compute_virial(positions, box)
        return density*kB_T + W/volume
