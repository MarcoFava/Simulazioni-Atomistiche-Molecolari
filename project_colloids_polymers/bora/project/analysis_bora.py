from src.simulation_module import initial_parameters, simulation
import pantarei as rei
import numpy as np

job = rei.Parallel(simulation)

r_coll = 1.8
r_coll = 1.2

n_part = 800
n_part = 600
n_part = 400
n_part = 2200


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

# r_coll = 1.8, n=800-400, steps=6e6
# seeds=[120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139]
# r_coll = 1.8, n=600, steps=6e6
# seeds=[120,123,124,125,126,127,128,129,130,132,133,134,135,136,137,138,139]

# r_coll = 1.8, n=800, change initial dist from 1+0.1q to 1+q
#seeds=[220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239] 

# r_coll = 1.2, n=2200, steps=6e6, initial dist back to 1+0.1q
seeds=[320,321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339]

# pic_n = '0' # n=800
# pic_n = '1' # n=600
# pic_n = '2' # n=400
# pic_n = '3' # change initial dist from 1+0.1q to 1+q, n=800
# pic_n = '4' # pic_n 0+3, n=800
pic_n = '5' # r_coll = 1.2, n=2200

# ['accept_ratio', 'accept_ratio_coll', 'box', 'tot_time', 'effective time', 'trajectory', 'packing_fract_coll', 'packing_fract_part', 'initial parameters']

frames = int(len(seeds)*steps/1000)
trajectory = np.ndarray((frames,3,n_coll+n_part))
tot_time = []
eff_time = []
print(trajectory.shape)
for seed in seeds:
    try:
        data = job(steps=steps,r_coll=r_coll,r_part=r_part,n_coll=n_coll,n_part=n_part,box=box,seed=seed)
    except:
        print(f'Problems with seed {seed}')
        continue
    tot_time.append(data['tot_time'])
    eff_time.append(data['effective time'])

    index = int(seeds.index(seed)*steps/1000)
    trajectory[index:index+int(steps/1000),:,:] = data['trajectory']
    print(f'Done seed {seed}')
print(np.mean(tot_time)/60,np.mean(eff_time))



from src.simulation_module import plot_val_over_NMC,pack_frac,effective_potential
import matplotlib.pyplot as plt

#########################################################
packing_fract_coll = pack_frac(r=r_coll,N=n_coll,V=V)  ##
packing_fract_part = pack_frac(r=r_part,N=n_part,V=V)  ##
#########################################################

# Cumulate the distances between colloids to plot the histogram
dist_colloids = []
# for pos in trajectory:
for i in range(0,len(trajectory),1):
    pos = trajectory[i]
    rij = pos[:,1] - pos[:,0]
    # rij = pos[0,1] - pos[0,0]                  # AWA
    # rij -= box * np.round(rij / box) # minimal image convention
    # for k in range(3):
    #     if rij[k] > box/2: rij[k]-= box
    #     elif rij[k] < -box/2: rij[k]+= box
    dist = np.linalg.norm(rij)
    dist_colloids.append(dist/2/r_coll) # normalize by 2*r_coll


# First take a look at the time evolution of the distances between colloids
blocks = [steps*k for k in range(len(seeds))]
plot_val_over_NMC(val=dist_colloids,x_multipl=1000,blocks=blocks,ylabel=f"$r/2r_c$",figname='analysis/time_evol-'+pic_n+'.png')
# plot_val_over_NMC(val=dist_colloids,x_multipl=1000,blocks=blocks,ylabel=f"$r/2r_c$",title='Evolution of colloids distance during simulation',figname='analysis/time_evol-'+pic_n+'.png')


# Now plot the histogram
plt.figure()
bins = 40                                                   # change n of bins here
bins = 50                                                   # change n of bins here
d_bin = (max(dist_colloids) - min(dist_colloids))/bins
ttl = f"Histogram of the distances distances between colloids\nfor the values $\eta_p$={packing_fract_part:.2}, $\eta_c$={packing_fract_coll:.2}, q={q}"
pdr, dr, ax = plt.hist(dist_colloids,bins,density=True)
# plt.title(ttl)
plt.ylabel(f"$p(r)$")
plt.xlabel(f"$r/2r_c$")
plt.ylabel(f"$p(x)$")                  # AWA
plt.xlabel(f"$x$")                  # AWA
plt.savefig('analysis/Histogram-'+pic_n+'.png')
# plt.savefig('analysis/Histogram-x.png')



# Compute now the g(r)
# g(r) prop to p(r) / (4 pi r^2 dr)
gdr = [pdr[i]/(dr[i]**2 *4*np.pi*d_bin) for i in range(len(pdr))]
# gdr /= np.mean(gdr[len(gdr)//4*3:]) # make the last points to be around 1
gdr /= np.mean(gdr[6:12]) # make the last points to be around 1
# gdr /= np.mean(gdr[5]) # make the last points to be around 1
plt.figure()
plt.plot(dr[:-1],gdr,'o-',markersize=3,label='Simulation data')
plt.ylabel('g(r)')
plt.xlabel(f"$r/2r_c$")
plt.grid()
# plt.gca().set_yscale("log", base=10)

plt.savefig('analysis/gdr-'+pic_n+'.png')


# Compute and plot the effective potential
Ueff = -np.log(gdr)
plt.figure()
plt.plot(dr[:-1],Ueff,'o-',markersize=3,label='Simulation data')


# Draw the 2 theoretical lines
r_coll_th = 0.5                     # the data was normalized by 2/r_coll
r_part_th = q*0.5                   # the data was normalized by 2/r_coll
x = np.linspace(0,box/r_coll,1000)
V_reduced = V / (2*r_coll)**3       # the data was normalized by 2/r_coll
density = (n_coll+n_part)/V_reduced

pot_full = [effective_potential(r,r_coll=r_coll_th,r_part=r_part_th,
                                N=(n_coll+n_part),V=V_reduced) for r in x]
pot_AO = [effective_potential(r,r_coll=r_coll_th,r_part=r_part_th,
                              N=(n_coll+n_part),V=V_reduced,which_pot='AO') for r in x]

# plt.plot(x,pot_full,label='Depletion potential')
plt.plot(x,pot_AO,label='AO potential')

ddata = max(Ueff) - min(Ueff)
ddr = max(dr[:-1]) - min(dr[:-1])
plt.ylim([min(Ueff)-ddata/20,max(Ueff)+ddata/20])
plt.xlim([min(dr[:-1])-ddr/20,max(dr[:-1])+ddr/20])

plt.grid(True,which="both")
plt.legend()
ttl = f"Effective potential with theoretical predictions\nfor the values $\eta_p$={packing_fract_part:.2}, $\eta_c$={packing_fract_coll:.2}, q={q}"
plt.title(ttl)
plt.ylabel(f"$V'(r)$")
plt.xlabel(f"$r/2r_c$")

plt.savefig('analysis/Effective_pot-'+pic_n+'.png')
# plt.show()
