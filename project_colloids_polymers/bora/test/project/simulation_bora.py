from src.simulation_module import initial_parameters, simulation
import pantarei as rei

job = rei.Parallel(simulation)

r_coll = 1.8
n_part = 800


steps = 1e6


##########################
q = 0.2                 ##
r_part = r_coll*q       ##
box = 7                 ##
V = box**3              ##
n_coll = 2              ##
##########################



seeds = [101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119]

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
