from f2py_jit import jit

class Interaction:
    def __init__(self, parameters, potential=None, files=None,
                 flags='-O3 -ffast-math'):

        self.potential = potential
        
        if files is not None:
            self.fiels = files
        else:
            potential_file = f'src/potential_module_{potential}.f90'
            interaction_file = 'src/interaction_module.f90'
            self.files = [potential_file, interaction_file]
            try:
                open(potential_file)
                open(interaction_file)
            except:
                print('ERROR: Unable to open the files ' + potential_file + ' and/or ' + interaction_file)

        self.flags = flags
        
        self.f90 = jit(self.files, inline=True, flags=self.flags)

        # initialize potential
        self.f90.potential_module.initialize(**parameters)


    def apply_pbc(self,position, box):
        self.f90.interaction_module.pbc(position,box)

    def apply_full_pbc(self,positions, box):
        self.f90.interaction_module.full_pbc(positions,box)
        
    def compute_single_energy(self, positions, box, index):
        return self.f90.interaction_module.compute_single_interaction(positions,box,index+1)
    
    def compute_energy(self, positions, box):
        return self.f90.interaction_module.compute_interaction(positions,box)

    def compute_virial(self, positions, box):
        return self.f90.interaction_module.compute_virial(positions,box)
    
    def compute_pressure(self, positions, box, density, kB_T, volume):
        W = self.compute_virial(positions, box)
        return density*kB_T + W/volume
    
    def check_overlap(self,positions,box,n_coll):
        return self.f90.interaction_module.check_overlap(positions,box,n_coll)
        # return self.f90.interaction_module.check_overlap(positions,box)
    
    def check_single_overlap(self,positions,box,n_coll,index):
        return self.f90.interaction_module.check_single_overlap(positions,box,n_coll,index+1)
        # return self.f90.interaction_module.check_single_overlap(positions,box,index+1)
    
    def check_overlap_all(self,positions,box):
        return self.f90.interaction_module.check_overlap_all(positions,box)
    
