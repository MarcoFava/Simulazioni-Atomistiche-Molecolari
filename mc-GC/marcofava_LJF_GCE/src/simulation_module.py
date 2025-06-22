import numpy as np
import matplotlib.pyplot as plt
import random as rnd
# import MonteCarlo as mc
# import interaction
# import system as sys




def simulation(temperature, density, steps, burnin_steps, param_interact=None, param_sys=None, mc_par_over_box=0.1, configuration='fcc'):
    import copy
    import src.MonteCarlo as mc
    import src.interaction as interaction
    import src.system as sys

#########################################################
# setup default parameters
    parameters_interaction = dict()
    if param_interact is None:
        parameters_interaction['sigma'] = 1.0       # Unità di misura di lunghezza
        parameters_interaction['epsilon'] = 1.0     # Unità di misura di energia

        parameters_interaction['rcut'] = 2.5*1.0    # Distanza di taglio per il potenziale cut and shift
    else:
        parameters_interaction = param_interact

    parameters_system = dict()
    # configuration = 'fcc'
    if param_sys is None:
        M = 6                                               # Cells per side
        a = 1.0                                             # Lattice spacing
        box = 4
        n = box**3*4                                        # 4, 32, 108, 256, 500, 864, 1372
        if configuration == 'cubic':
            parameters_system['M'] = M
            parameters_system['a'] = a
        elif configuration == 'fcc':
            parameters_system['n'] = n
            parameters_system['box'] = box
    else:
        parameters_system = param_sys

    # temperature = [].append(temperatures)
    # energies = dict()
    # energies[f'{temperature:.2}'], virials[f'{temperature:.2}'], slopes[f'{temperature:.2}'], intecepts[f'{temperature:.2}'], rvalues[f'{temperature:.2}'] = [], [], [], [], []


#########################################################
# create classes
    interact = interaction.Interaction(parameters_interaction, potential='lj')#,flags='-fbounds-check')

    my_sys = sys.System(temperature, interact, parameters_system, config=configuration)
    my_sys.density = density
    montecarlo = mc.MonteCarlo(my_sys, interact)
    rnd.seed(0)
    # mc_parameter = 0.6/11*my_sys.box[0]
    mc_parameter = mc_par_over_box*my_sys.box[0]


#########################################################
# run montecarlo until burnin without saving values
    P, E, W = [], [], []
    save_pos = []#.append(my_sys.positions)
    save_pos.append(copy.deepcopy(my_sys.positions))
    for i in range(int(burnin_steps)):
        montecarlo.run(mc_parameter)
        # if i % 100 == 0:
            # E.append(my_sys.energy)
            # P.append(my_sys.pressure)
            # W.append(my_sys.virial)
    # reset acceptance ratio 
    temp_accep_ratio = montecarlo.accept_ratio
    montecarlo.n_accepted_moves = 0
    montecarlo.tot_mc_moves = 0

    # run montecarlo until end of simulation and save every 100 mc cycles
    for i in range(int(steps - burnin_steps)):
        montecarlo.run(mc_parameter)
        if i % 100 == 0:
            E.append(my_sys.energy)
            # P.append(my_sys.pressure)
            W.append(my_sys.virial)


        if i % ((steps - burnin_steps)/10) == 0:
            save_pos.append(copy.deepcopy(my_sys.positions))
            my_sys.check_pbc()
            # print(f'Cycle {i}/{int(steps)}')
            # print(f'{np.abs((E[-1]-interact.compute_energy(my_sys.positions,my_sys.box)/my_sys.current_N)/E[-1])*100:.2e}')


#########################################################
# return data
    # print(f'temperature: {temperature:.2}, density: {density:5.3}, accept ratio: {montecarlo.accept_ratio:.2}')
    
    # return a dict of arrays
    P, E, W = np.array(P), np.array(E), np.array(W)
    P = temperature*density + W
    return {'pressure': P, 'energy': E, 'virial': W,
            'mean_pressure': P.mean(), 'mean_energy': E.mean(),
            'accept_ratio': montecarlo.accept_ratio, 'box': my_sys.box,
            'positions': save_pos}

    

def plot_val_over_NMC(val,title=None, label=None,figname=None, x_multipl=1, ylabel=None):
    x = np.linspace(0,len(val)*x_multipl,len(val))

    plt.figure()

    if label is None:
        # plt.plot(x,val,label=f'average={np.average(val):.3}')
        plt.plot(x,val)
    else:
        plt.plot(x,val,label=label)
        plt.legend()

    if ylabel is not None:
        plt.ylabel(ylabel)
        plt.xlabel('N_MC')
        
    plt.grid()

    if title is not None: plt.title(title)
    if figname is not None: plt.savefig(figname)

def visualize_3dmol(positions, cell=None, colors=None, radii=None, chain=False,
                    center=None, color='white', radius=0.5):
    """
    Visualize a particle configuration using 3dmol http://3dmol.csb.pitt.edu/
    """
    import py3Dmol
    if center is None:
        center = [0.5] * len(cell)
    if colors is None:
        colors = [color] * len(positions)
    if radii is None:
        radii = [radius] * len(positions)
    view = py3Dmol.view()
    view.setBackgroundColor('white')
    for i in range(len(positions)):
        view.addSphere({'center': {'x': positions[i][0],
                                   'y': positions[i][1],
                                   'z': positions[i][2]},
                        'radius': radii[i], 'color': colors[i]})
    if chain:
        for i in range(1, len(positions)):
            view.addCylinder({'start': {'x': positions[i-1][0],
                                        'y': positions[i-1][1],
                                        'z': positions[i-1][2]},
                              'end': {'x': positions[i][0],
                                      'y': positions[i][1],
                                      'z': positions[i][2]},
                              'color': 'white'})
    if cell is not None:
        view.addBox({'center': {'x': center[0], 'y': center[1], 'z': center[2]},
                     'dimensions': {'w': cell[0], 'h': cell[1], 'd': cell[2]},
                     'wireframe': True, 'color': "#000000"})
    return view
