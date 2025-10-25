import numpy as np

class System:

    def __init__(self, interaction, parameters, dim=3, config='rnd-coll', positions=None, chem_spec=None):
        self.temperature = None
        self.interaction = interaction
        self.dim = dim

        if positions is None:
            self.initialize_lattice(config,parameters)
        else:
            self.positions = positions
        
        self.current_N = self.positions.shape[1]

        # self._density = self.current_N/self.volume
        # self._energy = self.interaction.compute_energy(self.positions[:,:self.current_N],self.box)


################################################################
    @property
    def volume(self):
        return np.prod(self.box)

    @property
    def packing_fract_coll(self):
        if self.chem_spec is not None: return 4/3*np.pi * self.r_coll**3 * self.n_coll  / self.volume

    @property
    def packing_fract_part(self):
        if self.chem_spec is not None: return 4/3*np.pi * self.r_part**3 * self.n_part  / self.volume

    # @property
    # def density(self):
    #     self._density = self.current_N/self.volume
    #     return self._density
    
    # @density.setter
    # def density(self, value):
    #     # re-scale the box and positions to match the new density, compute again also the energy
    #     factor = 1/(value/self._density)**(1/3) # k = rho'/rho -> factor = 1/k**(1/3)
    #     self.box *= factor                      # box' = box/k**(1/3)
    #     self.positions *= factor                # pos' = pos/k**(1/3)
    #     self._density = value
    #     self._energy = self.interaction.compute_energy(self.positions[:,:self.current_N],self.box)
    #     return self._density

    # @property
    # def pressure(self):
    #     return self.density*self.temperature + self.virial

    # @property
    # def energy(self):
    #     return self._energy/self.current_N
    
    # @property
    # def virial(self):
        return self.interaction.compute_virial(self.positions,self.box)/self.volume
        # return self.interaction.compute_virial(self.positions[:,:self.current_N],self.box)/self.volume


################################################################
    def initialize_lattice(self,config,parameters):

        if (config=='cubic'):
            self.M = parameters['M']
            self.a = parameters['a']
            self.box = np.ones(self.dim)*self.M*self.a
            pos = []
            for i in range(self.M):
                for j in range(self.M):
                    for k in range(self.M):
                        pos.append([i*self.a - self.box[0]/2, j*self.a - self.box[1]/2, k*self.a - self.box[2]/2])
            self.positions =  np.array(pos).T

        elif (config=='fcc'):
            self.box = np.ones(self.dim)*parameters['box']
            self.positions = self.fcc_positions(parameters['n'],parameters['box']).T

        elif (config=='rnd-coll'):
            self.box = np.ones(self.dim)*parameters['box']
            self.n_coll = parameters['n_coll']
            self.n_part = parameters['n_part']
            self.r_coll = parameters['r_coll']
            self.r_part = parameters['r_part']
            self.chem_spec = np.array([0]*self.n_coll + [1]*self.n_part)
            self.positions = self.rnd_coll_part_positions(parameters['n_coll'],parameters['n_part'],self.box)

        else:
            raise ValueError("Not supported initial configuration \nPossible supported configurations: \'cubic\', \'fcc\'.\nPass otherwise the positions directly.")

################################################################
    def check_pbc(self):
        for i in range(len(self.box)):
            check = self.positions[i,:self.current_N] > self.box[i]
            if(check.sum()>0):
                print(f"There are {check.sum()} elements not respecting PBCs along the axes {i}")

################################################################
    def fcc_positions (self, n, box ):
        """Sets up the fcc lattice: four molecules per unit cell."""

        from itertools import product

        # Arguments are the number of particles, box length,

        nc = np.rint ( (n/4)**(1.0/3.0) ).astype(np.int_)
        assert n==4*nc**3, "{}{:d} -> {:d}".format('n, nc mismatch ',n,4*nc**3)
        cell = box / nc  # Unit cell
        box2 = box / 2.0 # Half box length
        r = np.empty((n,3),dtype=np.float64)

        r_fcc = np.array ( [ [0.25,0.25,0.25],[0.25,0.75,0.75],[0.75,0.75,0.25],[0.75,0.25,0.75] ], dtype=np.float64 )

        i = 0
        
        for ix, iy, iz in product(range(nc),repeat=3): # triple loop over unit cells
            for a in range(4): # loop over atoms in unit cell
                r[i,:] = r_fcc[a,:] + np.array ( [ix,iy,iz] ).astype(np.float64) # in range 0..nc
                r[i,:] = r[i,:] * cell                                          # in range 0..box
                r[i,:] = r[i,:] - box2                                          # in range -box2..box2
                # if not soft:
                #     assert not overlap ( r[i,:], e[i,:], r[:i,:], e[:i,:], box, length ), "Density too high"
                i = i + 1

        return r

################################################################
    def rnd_coll_part_positions(self, n_coll, n_part, box):
        import random as rnd
        iter = 0
        positions = np.ones((box.shape[0],n_coll+n_part),order='F')*2*box[0] # position all out of the box

        while True:
            # positions[:,0:n_coll] = (np.random.random((box.shape[0],n_coll))-0.5)*box[0]
            # positions[:,0:n_coll] = np.array([[0,0],[0,0],[0,2*self.r_coll]])
            # distance between colloids set to 2*(1 + 0.1*q)*r_coll 
            k = 2/np.sqrt(3)*(self.r_coll+0.1*self.r_part)
            positions[:,0:n_coll] = np.array([[0,0,0],[k,k,k]]).T + np.ones((3,2))*(-box[0]/2+self.r_coll)
            
            if not self.interaction.check_overlap(positions,box,n_coll): break
            iter += 1
            if iter >= 1000: 
                raise KeyError('Error: too many interations when assigning random positions to the colloids')

        i=n_coll
        while i < n_part+n_coll:
            positions[:,i] = [(rnd.random()-0.5)*box[i] for i in range(box.shape[0])]
            if not self.interaction.check_single_overlap(positions,box,n_coll,i): i += 1
            iter += 1
            if iter >= n_part*10: raise KeyError('Error: too many interations when assigning random positions to the particles')

        # check pbc
        for i in range(len(box)):
            check = positions[i,:] > box[i]
            if(check.sum()>0):
                raise KeyError(f"There are {check.sum()} elements not respecting PBCs along the axes {i}")

        return positions

################################################################
    # def save_pos(self,filename='./save_data/temp_pos.txt',coll=True,part=False):
    #     with open(filename,'a+') as f:
    #         if coll: 
    #             for pos in self.positions[:,:self.n_coll].T:
    #                 f.write(f'{pos[0]:22}{pos[1]:22}{pos[2]:22}\n')
    #         if part:
    #             for pos in self.positions[:,self.n_coll+1:]:
    #                 f.write(f'{pos}\n')
    #         if coll or part: 
    #             f.write('##############################################\n\n')
