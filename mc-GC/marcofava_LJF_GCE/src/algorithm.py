import numpy as np
import random as rnd
import time

import src.system as sys
import src.interaction as interaction
import src.MonteCarlo as mc

import pantarei as rei
scheduler = rei.scheduler.ThreadScheduler(backend='process')

class Algorithm:
    def __init__(self):
        self.initialize_parameters(potential='lj_c')


    def initialize_parameters(self,potential='lj_cs',mc_par_over_box=0.1,
                              param_interact=None,param_sys=None,configuration='fcc'):
        # setup default parameters
        self.potential = potential
        self.mc_par_over_box=mc_par_over_box
        self.configuration = configuration

        self.parameters_interaction = dict()
        if param_interact is None:
            self.parameters_interaction['sigma'] = 1.0       # Unità di misura di lunghezza
            self.parameters_interaction['epsilon'] = 1.0     # Unità di misura di energia

            self.parameters_interaction['rcut'] = 2.5*1.0    # Distanza di taglio per il potenziale cut and shift
        else:
            self.parameters_interaction = param_interact

        self.parameters_system = dict()
        if param_sys is None:
            self.parameters_system['tot_N'] = 1000
            M = 5                                               # Cells per side
            a = 1.0                                             # Lattice spacing
            box = 4
            n = box**3*4                                        # 4, 32, 108, 256, 500, 864, 1372
            if configuration == 'cubic':
                self.parameters_system['M'] = M
                self.parameters_system['a'] = a
                self.N = M**3
            elif configuration == 'fcc':
                self.parameters_system['n'] = n
                self.parameters_system['box'] = box
                self.N  = n
            # else:
            #     raise KeyError('Possible default configuration: fcc, cubic')
        else:
            self.parameters_system = param_sys

    
    # def simulation_GC(self,temperature,volume,z,steps,seed=None, which_vals=None):
    def simulation_GC(self,kwargs):
        temperature     = kwargs.get('temperature')
        volume          = kwargs.get('volume')
        z               = kwargs.get('z')
        eta             = kwargs.get('eta',None)
        # p_N             = kwargs.get('p_N')
        steps           = kwargs.get('steps')
        burnin           = kwargs.get('burnin',0)
        seed            = kwargs.get('seed',None)
        which_vals      = kwargs.get('which_vals', {'energy': True, 'pressure': True, 'density': True, 'N': False})

        interact = interaction.Interaction(self.parameters_interaction, potential='lj_c')#,flags='-fbounds-check')
        my_sys = sys.System(temperature, z, interact, self.parameters_system, config=self.configuration)
        my_sys.density = self.N/volume
        montecarlo = mc.MonteCarlo(my_sys, interact, eta=eta, seed=seed)
        # self.mc_param = self.mc_par_over_box * my_sys.box[0]
        # montecarlo.mc_param = 0.15
        # if eta is not None: montecarlo.weight_d = 0


        tot_time = 0
        t0 = time.time()                ## timer start ##
        P,E,d,Ns = [], [], [], []
        perc_prog = int(steps/20)

        for i in range(int(steps)):
            montecarlo.run()
            if i % perc_prog == 0: 
                if i >= burnin:
                    print(f'Progress {montecarlo.seed}: {int(i/steps*100)}%, ' +
                          f'N={my_sys.current_N}, d={my_sys.density:.3}')
                else:
                    print(f'Progress {montecarlo.seed}: burnin-{int(i/steps*100)}%, ' +
                          f'N={my_sys.current_N}, d={my_sys.density:.3}')
            if i % 100 == 0 and i >= burnin:
                if which_vals['energy']: E.append(my_sys.energy)
                if which_vals['pressure']: P.append(my_sys.pressure)
                if which_vals['density']: d.append(my_sys.density)
                if which_vals['N']: Ns.append(my_sys.current_N.real)

                t1 = time.time()        ## timer lap ##
                tot_time += (t1-t0)/my_sys.current_N
                t0 = t1
        
        print(f'Progress {montecarlo.seed}: ENDED, ' +
                f'N={my_sys.current_N}, d={my_sys.density:.3}, accept:{montecarlo.accept_ratio:.2}')
        
        acceptance = dict()
        acceptance['accept_ratio'] = montecarlo.accept_ratio
        acceptance['accepted_moves_d'] = montecarlo.accepted_moves_d
        acceptance['accepted_moves_i'] = montecarlo.accepted_moves_i
        acceptance['accepted_moves_r'] = montecarlo.accepted_moves_r
        acceptance['tot_mc_moves'] = montecarlo.tot_mc_moves

        

        return {'pressure':P, 'energy': E, 'density': d, 'N': Ns,
                'tot_time': tot_time/steps,
                'acceptance': acceptance}
                # 'system': my_sys, 'montecarlo': montecarlo}


    # def run_simulation(self,parallel=True,clear_first=False, wait=False,**kwargs):
    #     task = rei.Task(self.simulation_GC,clear_first=clear_first)
    #     if parallel: job = rei.Thread(task, scheduler=scheduler)
    #     else: job = task

    #     data = job(kwargs=kwargs)
    #     if wait: return 0
    #     scheduler.wait()
    #     return data

    # def plot_hist(self,data):
    #     import matplotlib.pyplot as plt
    #     Ns = data['N']

    #     plt.hist(Ns)