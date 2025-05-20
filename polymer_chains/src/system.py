import numpy as np

class System:
    def __init__(self,N,L):
        self.N = N
        self.L = L
        self.positions, self.lattice = self.setup(N,L)

    @property
    def centre_of_mass(self):
        # return sum(np.array(self.positions)/self.N)
        return np.array(self.positions).mean(axis=0)
        # return self.positions.mean()
    
    @centre_of_mass.setter
    def centre_of_mass(self,cm):
        dr = cm-self.centre_of_mass
        self.shift_centre_of_mass(dr=dr)
        return cm

    def shift_centre_of_mass(self,dr):
        for i in range(np.shape(self.positions)[1]):
            for j in range(np.shape(self.positions)[0]):
                self.positions[j][i] += int(dr[i])
            



    def setup(self, N, L):
        """
        Setup the initial conformation as a linear chain wrapped around
        the cell boundaries
        """
        lattice = np.zeros((L+1, L+1), dtype=int)
        positions = []
        for x in range(L):
            for y in range(L):
                if x % 2 != 0:
                    y = L - y
                positions.append([x, y])
                lattice[x, y] = 1
                if len(positions) >= N:
                    break
            if len(positions) >= N:
                break
        return positions, lattice


    def unfold(self):
        import copy
        _positions = copy.deepcopy(self.positions)

        for i in range(1, len(_positions)):
            for alpha in range(len(self.positions[0])):
                # If the distance between successive monomers is larger than 1
                # we have crossed the boundary of the central cell.
                # In that case we unfold all successive monomers along the chain
                if _positions[i][alpha] - _positions[i-1][alpha] > 1:
                    for j in range(i, len(_positions)):
                        _positions[j][alpha] -= self.L
                if _positions[i][alpha] - _positions[i-1][alpha] < -1:
                    for j in range(i, len(_positions)):
                        _positions[j][alpha] += self.L
        return _positions