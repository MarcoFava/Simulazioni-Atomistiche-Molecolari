def mc(density, temperature, steps, burnin_steps):
    interaction = src.interaction.Interaction(parameters, which_potential='lj')#,flags='-fbounds-check')

    P, E = [], []
    my_sys = sys.System(temperature, interaction, M=M, a=a, config='crystal')
    my_sys.density = density
    montecarlo = mc.MonteCarlo(my_sys, interaction)
    rnd.seed(0)
    mc_parameter = 0.4/11*my_sys.box[0]

    # run montecarlo until burnin without saving values
    for _ in range(burnin):
        montecarlo.run(mc_parameter)
        # P.append(my_sys.pressure)
        # E.append(my_sys.energy)

    # reset acceptance ratio 
    temp_accep_ratio = montecarlo.accept_ratio
    montecarlo.n_accepted_moves = 0
    montecarlo.n_mc_moves = 0

    # run montecarlo until end of simulation and save every 100 mc cycles
    for i in range(int(N_MC - burnin)):
        montecarlo.run(mc_parameter)
        # if i//100 == i/100:
        # COSLO: more explicit
        if i % 100 == 0:
            P.append(my_sys.pressure)
            E.append(my_sys.energy)

    # COSLO: here we return a dict of arrays
    P, E = numpy.array(P), numpy.(E)
    return {'pressure': P, 'energy': E,
            'mean_pressure': P.mean(), 'mean_energy': E.mean()}
