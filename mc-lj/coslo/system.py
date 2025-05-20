import numpy as np

class System:

    def __init__(self, temperature, interaction, M=2, a=1.0, dim=3, config='crystal', positions=None):
        self.temperature = temperature
        self.interaction = interaction
        self.M = M
        self.a = a
        self.box = np.array([M*a,M*a,M*a])

        if positions is None:
            self.positions = self.initialize_lattice(config=config)
        else:
            self.positions = positions
        
        self._density = self.positions.shape[1]/self.volume
        self._energy = self.interaction.compute_total_energy(self.positions,self.box)

        self.track_pressure = []
        self.track_energy = []
        self.track_values()

    @property
    def volume(self):
        V = 1
        for i in range(self.box.shape[0]):
            V *= self.box[i]
        return V

    @property
    def density(self):
        return self._density
    
    @density.setter
    def density(self, value):
        # re-scale the box and positions to match the new density, compute again also the energy
        factor = 1/(value/self._density)**(1/3) # k = rho'/rho -> factor = 1/k**(1/3)
        self.box *= factor                      # box' = box/k**(1/3)
        self.positions *= factor                # pos' = pos/k**(1/3)
        self.a *= factor
        self._density = value
        self._energy = self.interaction.compute_total_energy(self.positions,self.box)
        return self._density

    @property
    def pressure(self):
        W = self.interaction.compute_virial(self.positions,self.box)
        return self.density*self.temperature + W/self.volume

    @property
    def energy(self):
        return self._energy/np.shape(self.positions)[1]


    ################################################################
    def initialize_lattice(self,config='crystal'):
        "Build a cubic crystal"
        pos = []
        if (config=='crystal'):
            for i in range(self.M):
                for j in range(self.M):
                    for k in range(self.M):
                        pos.append([i*self.a - self.box[0]/2, j*self.a - self.box[1]/2, k*self.a - self.box[2]/2])
            return np.array(pos).T
        else:
            raise ValueError("Not supported initial configuration \nPossible supported configurations: \'crystal\'")

    ################################################################
    # COSLO: I do not think it is the system's responsibility to track its energy/pressure.
    # It is up to the simulation to do this kind of things. The System has no notion
    # of "time evolution" - the simulation or the trajectory has it.
    def track_values(self, is_pressure=True, is_energy=True):
        if is_pressure: self.track_pressure.append(self.pressure)
        if is_energy: self.track_energy.append(self.energy)

    ################################################################
    def check_pbc(self):
        # pbc implemented when computing the energy of interaction in energy_module.f90 
        for i in range(len(self.box)):
            check = self.positions[i,:] > self.box[i]
            # COSLO: again () not needed
            if(check.sum()>0):
                print(f"There are {check.sum()} elements not respecting PBCs along the axes {i}")

