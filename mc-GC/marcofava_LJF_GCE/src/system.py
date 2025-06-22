import numpy as np

class System:

    def __init__(self, temperature, z, interaction, parameters, dim=3, config='cubic', positions=None):
        self.temperature = temperature
        # self.chem_potential = chem_potential
        self.z = z
        self.interaction = interaction
        self.dim = dim

        if positions is None:
            self.initialize_lattice(config,parameters)
        else:
            self.positions = positions
        
        self.current_N = self.positions.shape[1]
        self.tot_N = parameters['tot_N']
        self.positions = np.hstack((self.positions,np.zeros([3,self.tot_N-self.current_N])),dtype=np.float64)
        self.positions = np.asfortranarray(self.positions)
        
        self._density = self.current_N/self.volume
        self._energy = self.interaction.compute_energy(self.positions[:,:self.current_N],self.box)


################################################################
    @property
    def volume(self):
        V = 1
        for i in range(self.box.shape[0]):
            V *= self.box[i]
        return V

    @property
    def density(self):
        self._density = self.current_N/self.volume
        return self._density
    
    @density.setter
    def density(self, value):
        # re-scale the box and positions to match the new density, compute again also the energy
        factor = 1/(value/self._density)**(1/3) # k = rho'/rho -> factor = 1/k**(1/3)
        self.box *= factor                      # box' = box/k**(1/3)
        self.positions *= factor                # pos' = pos/k**(1/3)
        self._density = value
        self._energy = self.interaction.compute_energy(self.positions[:,:self.current_N],self.box)
        return self._density

    @property
    def pressure(self):
        return self.density*self.temperature + self.virial

    @property
    def energy(self):
        return self._energy/self.current_N
    
    @property
    def virial(self):
        return self.interaction.compute_virial(self.positions[:,:self.current_N],self.box)/self.volume


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

        else:
            raise ValueError("Not supported initial configuration \nPossible supported configurations: \'cubic\', \'fcc\'")

################################################################
    def check_pbc(self):
        # pbc implemented when computing the energy of interaction in energy_module.f90 
        for i in range(len(self.box)):
            check = self.positions[i,:self.current_N] > self.box[i]
            if(check.sum()>0):
                print(f"There are {check.sum()} elements not respecting PBCs along the axes {i}")

################################################################
    def update_values(self):
        self._energy = self.interaction.compute_energy(self.positions[:,:self.current_N],self.box)

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
