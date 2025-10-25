from src.simulation_module import initial_parameters, simulation
import pantarei as rei

job = rei.Parallel(simulation)

##################
r_coll = 1.8    ##
n_part = 800    ##
#n_part = 600   ## ????
n_part = 400    ##
##################
r_coll = 1.2    ##
n_part = 2200   ##
##################


steps = 1e6
steps = 6e6


##########################
q = 0.2                 ##
r_part = r_coll*q       ##
box = 7                 ##
V = box**3              ##
n_coll = 2              ##
##########################


# r_coll = 1.8, n=800, steps=1e6
# seeds=[101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119] 

# r_coll = 1.8, n=800-600??-400, steps=6e6
# seeds=[120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139]

# r_coll = 1.8, n=800, steps=6e6, change initial dist from 1+0.1q to 1+q
#seeds=[220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239] 

# r_coll = 1.2, n=2200, steps=6e6, initial dist back to 1+0.1q
# seeds = [int(320+k) for k in range(20)]
seeds=[320,321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339]

for seed in seeds:
    data = job(steps=steps,r_coll=r_coll,r_part=r_part,n_coll=n_coll,n_part=n_part,box=box,seed=seed)
    print(data)






# import time
# from src.interaction    import Interaction
# from src.system         import System
# from src.MonteCarlo     import MonteCarlo

# kwars = dict(steps=steps,r_coll=r_coll,r_part=r_part,n_coll=n_coll,n_part=n_part,box=box)

# parameters                = initial_parameters(**kwars)
# parameters_interaction    = parameters['parameters_interaction']
# parameters_system         = parameters['parameters_system']
# configuration             = parameters['configuration']
# potential                 = parameters['potential']
# seed                      = parameters['seed']
# clear_first               = kwars.get('clear_first',True)

# # which_vals = kwars.get('which_vals', {'energy': False, 'pressure': False, 'density': False, 'N': False})

# # Create classes 
# # Possible flags: #,flags='-fbounds-check')#,flags='-Wall -Wextra -fcheck=all')
# interact      = Interaction(parameters_interaction, potential=potential,flags='-Wall -Wextra -fcheck=all')
# my_sys        = System(interact, parameters_system, config=configuration)
# montecarlo    = MonteCarlo(my_sys, interact, seed=seed)



# # Double check that there are no overlaps before starting the simulation
# if interact.check_overlap_all(my_sys.positions,my_sys.box): 
#     raise KeyError('Error: Overlap found after initialization')


# for _ in range(int(steps)):
#     montecarlo.run()

# print(my_sys.r_part,my_sys.volume/V,montecarlo.mc_param,montecarlo.accept_ratio)
