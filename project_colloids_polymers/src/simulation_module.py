import numpy as np
import matplotlib.pyplot as plt
import random as rnd
# import MonteCarlo as mc
# import interaction
# import system as sys



def initial_parameters(**kwargs):
    param_interact = kwargs.get('param_interact',None)
    param_sys = kwargs.get('param_sys',None)
    potential = kwargs.get('potential','hs_coll')
    seed = kwargs.get('seed',None)
    
    # setup default parameters
    parameters_interaction = dict()
    if param_interact is None:
        if potential == 'hs_coll':
            parameters_interaction['r_coll'] = kwargs.get('r_coll',0.5)         # Raggio colloidi
            parameters_interaction['r_part'] = kwargs.get('r_part',0.1)         # Raggio particelle
            parameters_interaction['n_coll'] = kwargs.get('n_coll',2)           # Numero colloidi
            parameters_interaction['n_part'] = kwargs.get('n_part',100)         # Numero particelle
        else:
            parameters_interaction['sigma'] = kwargs.get('sigma',1.0)           # Unità di misura di lunghezza
            parameters_interaction['epsilon'] = kwargs.get('epsilon',1.0)       # Unità di misura di energia

            default_rcut = 2.5*parameters_interaction['sigma']
            parameters_interaction['rcut'] = kwargs.get('rcut',default_rcut)    # Distanza di taglio per il potenziale cut and shift
    else:
        parameters_interaction = param_interact

    parameters_system = dict()
    configuration = kwargs.get('configuration','rnd-coll')
    if param_sys is None:
        M = 6                                               # Cells per side
        a = 1.0                                             # Lattice spacing
        box = 4
        n = box**3*4                                        # 4, 32, 108, 256, 500, 864, 1372
        if configuration == 'cubic':
            parameters_system['M'] = kwargs.get('M',M)
            parameters_system['a'] = kwargs.get('a',a)
        elif configuration == 'fcc':
            parameters_system['n'] = kwargs.get('n',n)
            parameters_system['box'] = kwargs.get('box',box)
        elif configuration == 'rnd-coll':
            parameters_system['n_coll'] = kwargs.get('n_coll',2)
            parameters_system['n_part'] = kwargs.get('n_part',100)
            parameters_system['r_coll'] = kwargs.get('r_coll',0.5)
            parameters_system['r_part'] = kwargs.get('r_part',0.1)
            parameters_system['box'] = kwargs.get('box',10)
    else:
        parameters_system = param_sys


    return {
            'parameters_interaction': parameters_interaction, 
            'parameters_system': parameters_system,
            'configuration': configuration,
            'potential': potential,
            'seed' : seed
            }



def simulation(steps, burnin=0,**kwars):
    import copy
    import time
    from src.interaction    import Interaction
    from src.system         import System
    from src.MonteCarlo     import MonteCarlo
    from src.trajectory     import Trajectory


    parameters                = initial_parameters(**kwars)
    parameters_interaction    = parameters['parameters_interaction']
    parameters_system         = parameters['parameters_system']
    configuration             = parameters['configuration']
    potential                 = parameters['potential']
    seed                      = parameters['seed']
    clear_first               = kwars.get('clear_first',True)
    # which_vals = kwars.get('which_vals', {'energy': False, 'pressure': False, 'density': False, 'N': False})
    
    # Create classes 
    # Possible flags: #,flags='-fbounds-check')#,flags='-Wall -Wextra -fcheck=all')
    interact      = Interaction(parameters_interaction, potential=potential)
    my_sys        = System(interact, parameters_system, config=configuration)
    montecarlo    = MonteCarlo(my_sys, interact, seed=seed)
    trajectory    = Trajectory(sim_params=parameters,mode='w',clear_first=clear_first)

    trajectory_dir = trajectory.trajectory_dir
    print(f'seed: {seed}, trajectory_dir:\'{trajectory_dir}\'\n')

    # Double check that there are no overlaps before starting the simulation
    if interact.check_overlap_all(my_sys.positions,my_sys.box): 
        raise KeyError('Error: Overlap found after initialization')


    #########################################################
    perc_prog = int(steps/20) # Print output every 5% of the total simulation
    # save_pos = []
    # dr_coll = []

    # keep track of the time spent in the simulation
    tot_time = 0
    t0 = time.time()                ## timer start ##

    for i in range(int(steps)):
        montecarlo.run()

        # Every 100 steps, keep track of the positions in the trajectory
        if i % 100 == 0 and i >= burnin:
            trajectory.write(my_sys.positions,i)
            # save_pos.append(my_sys.positions.copy())
            # dr_coll.append(np.linalg.norm(my_sys.positions[0] - my_sys.positions[1]))

        # Every 5% of the total run, print some progress
        if i % perc_prog == 0: 
            if i >= burnin:
                tot_time += time.time()-t0
                t0 = time.time()
                print(f'Progress s-{montecarlo.seed}: {int(i/steps*100)}%-{tot_time/60:.1f}\', ' +
                        f'packing frac: coll={my_sys.packing_fract_coll:.3}, part={my_sys.packing_fract_part:.3}'+ 
                        f', accept ratio colloids: {montecarlo.accept_ratio_coll:.2}')
            else:
                print(f'Progress {montecarlo.seed}: burnin->{int(i/steps*100)}%-{tot_time/60:.1f}\'')

    tot_time += time.time()-t0
    print(f'Progress s-{montecarlo.seed}: ENDED, ' +
          f'tot_time={tot_time:.1f}s, effective time={tot_time/steps/my_sys.current_N:.1}s, ' +
          f'accept ratio: {montecarlo.accept_ratio:.2}, accept ratio colloids: {montecarlo.accept_ratio_coll:.2}'
        #   f'packing frac: coll={my_sys.packing_fract_coll:.3}, part={my_sys.packing_fract_part:.3}'
          )
    

    # Read again the positons to be returned
    save_pos = []
    for pos in trajectory:
        save_pos.append(pos)


    return {
            'accept_ratio': montecarlo.accept_ratio,'accept_ratio_coll': montecarlo.accept_ratio_coll, 'box': my_sys.box
            ,'tot_time': tot_time, 'effective time':tot_time/steps/my_sys.current_N
            ,'trajectory': save_pos
            ,'packing_fract_coll':my_sys.packing_fract_coll,'packing_fract_part':my_sys.packing_fract_part
            ,'initial parameters' : parameters
            ,'trajectory_dir' : trajectory_dir
            }


# def compute_gdr(positions,box,dr,n_bins):
#     rmax = dr*n_bins
#     hbox = box/2
#     hist = 0

#     for i in range(size)



# def gdr(positions,box,n_bins,last_N):
#     dr = box/n_bins
#     counts = np.zeros(n_bins)
    
#     for pos in positions:
#         for i in range(last_N):
#             for j in range(last_N):
#                 if j>=i+1:
#                     r = min(np.prod(np.square(pos[:,j]-pos[:,i])),abs(np.prod(np.square(pos[:,j]-pos[:,i]))-box))
#                     print(int(r/dr),r/dr,r//dr)
#                     counts[int(r/dr)] += 1
#     return counts/counts.sum()


def plot_val_over_NMC(val,title=None, label=None,figname=None, x_multipl=1, ylabel=None,blocks=None):
    x = np.linspace(0,len(val)*x_multipl,len(val))

    plt.figure()

    if label is None:
        # plt.plot(x,val,label=f'average={np.average(val):.3}')
        plt.plot(x,val)
    else:
        plt.plot(x,val,label=label)
        plt.legend()
    
    if blocks is not None:
        y_min = min(val)
        y_max = max(val)
        print(y_min,y_max)
        for x_block in blocks:
            plt.axvline(x_block,ymin=y_min,ymax=y_max,color='r',linestyle='--')

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
        if center is None:
            center = [0.5] * len(cell)
            center = [0.] * len(cell)
        view.addBox({'center': {'x': center[0], 'y': center[1], 'z': center[2]},
                     'dimensions': {'w': cell[0], 'h': cell[1], 'd': cell[2]},
                     'wireframe': True, 'color': "#000000"})
    return view
