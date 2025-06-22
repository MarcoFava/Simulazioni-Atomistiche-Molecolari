import numpy as np
import random
# from numba import njit

def setup(N, L):
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

# @njit
def mc_move(positions, lattice):
    # Either take the first particle and add it to the tail (drop=0, insert=-1)
    # or the opposite (drop=-1, insert=0)
    drop, insert = random.choice([[0, -1], [-1, 0]])
    # drop_insert = [[0, -1], [-1, 0]]
    # choice = random.randint(0,1)
    # drop, insert = drop_insert[choice]

    # Attempt to put the particle at position `pos`
    delta = random.choice([[1, 0], [-1, 0], [0, 1], [0, -1]])
    # delta = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    # choice = random.randint(0,3)
    # delta = delta[choice]

    pos = positions[insert][:] + np.array(delta)

    # If the new site is empty, put the particle there
    if lattice[pos[0], pos[1]] == 0:
        new = positions.pop(drop)
        lattice[new[0], new[1]] = 0
        lattice[pos[0], pos[1]] = 1
        new[:] = pos

        # Either add to the head or the tail
        if insert == 0:
            positions.insert(0, new)
        else:
            positions.append(new)

def unfold(positions, L):
    import copy
    _positions = copy.deepcopy(positions)

    for i in range(1, len(_positions)):
        for alpha in range(len(positions[0])):
            # If the distance between successive monomers is larger than 1
            # we have crossed the boundary of the central cell.
            # In that case we unfold all successive monomers along the chain
            if _positions[i][alpha] - _positions[i-1][alpha] > 1:
                for j in range(i, len(_positions)):
                    _positions[j][alpha] -= L
            if _positions[i][alpha] - _positions[i-1][alpha] < -1:
                for j in range(i, len(_positions)):
                    _positions[j][alpha] += L
    return _positions

# def simulation(N,L,steps,burnin):
#     positions, lattice = setup(N, L)
#     R_values = []
#     for step in range(steps):
#         for attempt in range(len(positions)):
#             mc_move(positions, lattice)
#         if (step > burnin):
#             R_values.append(R(positions))
    
#     return positions, lattice, R_values

def simulation(N,L,steps,burnin):
    import src.system as sys
    import src.Montecarlo as mc
    # import system as sys
    # import Montecarlo as mc
    import copy

    my_sys = sys.System(N=N,L=L)
    montecarlo = mc.Montecarlo(system=my_sys)

    R_values, trajectory = [], []
    for step in range(steps):
        montecarlo.run()

        # if step >= burnin:
        if (step>0) & (step%burnin == 0):
            my_sys.centre_of_mass = [0,0]
            trajectory.append(copy.deepcopy(my_sys.positions))
            # R_values.append(R(my_sys.positions))
    
    return {'trajectory': trajectory
            # , 'R_values': R_values
            }

def prob_overlap(dr,data,clear_first=False):
    import src.system as sys
    import src.Montecarlo as mc
    # import system as sys
    # import Montecarlo as mc
    import copy
    # import pantarei as rei
    
    # task = rei.Task(overlap,clear_first=clear_first)

    N = len(data['trajectory'][0])
    L = 1000
    my_sys = sys.System(N=N,L=L)
    montecarlo = mc.Montecarlo(system=my_sys)

    
    no_overlaps = 0
    tot_overlap = 0
    for i in range(len(data['trajectory'])):
        #if i%1 == 0: print(i)
        for j in range(len(data['trajectory'])):
#            for deg in [90, 180, 270]:
            if i<=j : pass
            
            tot_overlap += 1

            my_sys.positions = data['trajectory'][i]
            cm = my_sys.centre_of_mass
            
            my_sys.positions = copy.deepcopy(data['trajectory'][j])
            my_sys.centre_of_mass = cm + dr
            
            if overlap(data['trajectory'][i],my_sys.positions) is None:
            # if task(positions_0=data['trajectory'][i],positions_1=my_sys.positions) is None:
                no_overlaps += 1
                # print(f'rate overlaps = {no_overlaps/tot_overlap:.3}, {i}, {j}')
                # print(overlap(data['trajectory'][i],my_sys.positions))
                # break
    print(no_overlaps/tot_overlap)
    return (tot_overlap - no_overlaps)/tot_overlap




def overlap(positions_0, positions_1):
    """
    Checks for overlaps between positions_0 and positions_1
    Returns the first overlap found (if found) or None (if not found)
    """
    where_overlap = None
    for i in range(np.shape(positions_0)[0]):
        is_overlap = False
        for j in range(np.shape(positions_1)[0]):
            if positions_0[i] == positions_1[j]:
                where_overlap = positions_0[i]
                is_overlap = True
                break
        if is_overlap: break
    return where_overlap


def R(positions):
    return np.sqrt(np.sum((np.array(positions[-1])-np.array(positions[0]))**2))

def visualize_chain(positions):
    import matplotlib.pyplot as plt
    cm = sum(np.array(positions)/np.shape(positions)[0])
    plt.plot([r[0] for r in positions], [r[1] for r in positions], '-o')
    plt.plot(cm[0],cm[1],'o')

def plot_val_over_NMC(val,title=None, label=None,figname=None, x_multipl=1):
    import matplotlib.pyplot as plt

    x = np.linspace(0,len(val)*x_multipl,len(val))

    plt.figure()

    if label is None:
        # plt.plot(x,val,label=f'average={np.average(val):.3}')
        plt.plot(x,val)
    else:
        plt.plot(x,val,label=label)
        plt.legend()

    plt.grid()

    if title is not None: plt.title(title)
    if figname is not None: plt.savefig(figname)
