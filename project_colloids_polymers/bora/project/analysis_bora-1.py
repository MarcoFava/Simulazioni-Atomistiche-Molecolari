from src.simulation_module import initial_parameters, simulation,plot_val_over_NMC,pack_frac,effective_potential
import pantarei as rei
import numpy as np
import matplotlib.pyplot as plt

job = rei.Parallel(simulation)


##########################
q = 0.2                 ##
box = 7                 ##
V = box**3              ##
n_coll = 2              ##
##########################


def import_dist_colloids(r_coll,n_part,steps,seeds):
    ##########################
    r_part = r_coll*q       ##
    ##########################

    frames = int(len(seeds)*steps/1000)
    trajectory = np.ndarray((frames,3,n_coll+n_part))
    print(trajectory.shape)
    for seed in seeds:
        try:
            data = job(steps=steps,r_coll=r_coll,r_part=r_part,n_coll=n_coll,n_part=n_part,box=box,seed=seed)
        except:
            print(f'Problems with seed {seed}')
            continue

        index = int(seeds.index(seed)*steps/1000)
        trajectory[index:index+int(steps/1000),:,:] = data['trajectory']


    dist_colloids = []
    for i in range(0,len(trajectory),1):
        pos = trajectory[i]
        rij = pos[:,1] - pos[:,0]
        rij -= box * np.round(rij / box) # minimal image convention
        dist = np.linalg.norm(rij)
        if dist <= 2*r_coll+2.2*r_part: dist_colloids.append(dist/2/r_coll) # normalize by 2*r_coll
        # if dist <= 2*r_coll+3*r_part: dist_colloids.append(dist/2/r_coll) # normalize by 2*r_coll
        # dist_colloids.append(dist/2/r_coll) # normalize by 2*r_coll
    
    return np.array(dist_colloids)


def plot_potential_hist(r_coll,n_part,steps,seeds,bins,axis,figure,color=None,which_pot='AO',corr=1):
    if color is None:
        col_dat = '#004aad'
        col_pot = 'red'
    else:
        col_dat = color
        col_pot = color

    dist_colloids = import_dist_colloids(r_coll,n_part,steps,seeds)

    r_part = q*r_coll
    packing_fract_coll = pack_frac(r=r_coll,N=n_coll,V=V)
    packing_fract_part = pack_frac(r=r_part,N=n_part,V=V)


    # Draw the 2 theoretical lines
    r_coll_th = 0.5                     # the data was normalized by 2/r_coll
    r_part_th = q*0.5                   # the data was normalized by 2/r_coll
    x = np.linspace(1,box/r_coll,1000)
    V_reduced = V / (2*r_coll)**3       # the data was normalized by 2/r_coll
    pot_full =  [effective_potential(r,r_coll=r_coll_th,r_part=r_part_th,
                                    N=(n_coll+n_part)*corr,V=V_reduced) for r in x]
    pot_AO =    [effective_potential(r,r_coll=r_coll_th,r_part=r_part_th,
                                  N=(n_coll+n_part)*corr,V=V_reduced,which_pot='AO') for r in x]
    # which_pot = None
    if which_pot=='FULL': 
        axis.plot(x,pot_full,ls='-.',color=col_pot)
        # axis.plot(x,pot_full,ls='-.',label='Depletion potential',color=col_pot)
    if which_pot=='AO': 
        axis.plot(x,pot_AO,ls='-.',label='AO potential',color=col_pot)


    # Get the histogram
    bins = bins                                                   # change n of bins here
    d_bin = (max(dist_colloids) - min(dist_colloids))/bins
    pdr, dr, bar = plt.hist(dist_colloids,bins,density=True)


    # Compute now the g(r)
    # g(r) prop to p(r) / (4 pi r^2 dr)
    gdr = [pdr[i]/(dr[i]**2 *4*np.pi*d_bin) for i in range(len(pdr))]
    # gdr /= np.mean(gdr[len(gdr)//4*3:]) # make the last points to be around 1
    # gdr /= np.mean(gdr[5:7]) # make the middle points to be around 1
    gdr /= np.mean(gdr[-3:])

    # Compute and plot the effective potential
    Ueff = -np.log(gdr)
    label=r'$\eta_{part}$'+f' = {packing_fract_part:.2}'
    # label=r'$\eta_{coll}$'+f' = {packing_fract_coll:.2}'
    axis.plot(dr[:-1],Ueff,'o-',markersize=3,label=label,color=col_dat)

        
    ddata = max(Ueff) - min(Ueff)
    ddr = max(dr[:-1]) - min(dr[:-1])
    ylim = [min(Ueff)-ddata/20,max(Ueff)+ddata/20]
    xlim = [min(dr[:-1])-ddr/20,max(dr[:-1])+ddr/20]
    return xlim,ylim


# Initialize figure
figsize = (360/96,270/96)
dpi = 200
fig, ax = plt.subplots(figsize=figsize,dpi=dpi)
fig_hist, ax_hist = plt.subplots(figsize=figsize,dpi=dpi)
# fig = plt.figure(figsize=figsize)
# ax=plt.gca()
# fig_hist = plt.figure(figsize=figsize)
# ax_hist=plt.gca()
xlims =[]
ylims =[]


# ######################
# r_coll = 1.8        ##
# n_part = 400        ##
# steps = 6e6         ###################################################################
# seeds=[120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139]
# xlim,ylim=plot_potential_hist(r_coll,n_part,steps,seeds,bins=25,axis=ax,figure=fig,color='#ad004A',which_pot='FULL',corr=1.3)
# xlims.append(xlim[0])
# xlims.append(xlim[1])
# ylims.append(ylim[0])
# ylims.append(ylim[1])

# ######################
# r_coll = 1.8        ##
# n_part = 600        ##
# steps = 6e6         ###################################################################
# seeds=[120,123,124,125,126,127,128,129,130,132,133,134,135,136,137,138,139]
# xlim,ylim=plot_potential_hist(r_coll,n_part,steps,seeds,bins=25,axis=ax,figure=fig,color='#ad6300',which_pot='FULL',corr=1.3)
# xlims.append(xlim[0])
# xlims.append(xlim[1])
# ylims.append(ylim[0])
# ylims.append(ylim[1])

######################
r_coll = 1.8        ##
n_part = 800        ##
steps = 6e6         ###################################################################
seeds=[120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139]
# seeds+=[220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239] 
xlim,ylim=plot_potential_hist(r_coll,n_part,steps,seeds,bins=25,axis=ax,figure=fig,color='#004aad',which_pot='FULL',corr=1.15)
xlims.append(xlim[0])
xlims.append(xlim[1])
ylims.append(ylim[0])
ylims.append(ylim[1])

# ######################
# r_coll = 1.2        ##
# n_part = 2200       ##
# steps = 6e6         ###################################################################
# seeds=[320,321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339]
# xlim,ylim=plot_potential_hist(r_coll,n_part,steps,seeds,bins=25,axis=ax,figure=fig,color='#004aad',which_pot='AO')
# xlims.append(xlim[0])
# xlims.append(xlim[1])
# ylims.append(ylim[0])
# ylims.append(ylim[1])



# ##################################################################################################
# #######################################        HIST        #######################################
# ##################################################################################################
# ax_hist.grid(True,which="both")
# ax_hist.set_ylabel(f"$p(r)$")
# ax_hist.set_xlabel(f"$r/2R$")
# fig_hist.tight_layout()
# fig_hist.savefig('analysis/Histogram-'+str(n_part)+'.png')


##################################################################################################
#######################################        Ueff        #######################################
##################################################################################################
# ddata = max(Ueff) - min(Ueff)
# ddr = max(dr[:-1]) - min(dr[:-1])
# ax.set_ylim([min(Ueff)-ddata/20,max(Ueff)+ddata/20])
# ax.set_xlim([min(dr[:-1])-ddr/20,max(dr[:-1])+ddr/20])
# ttl = f"Effective potential with theoretical predictions\nfor the values $\eta_p$={packing_fract_part:.2}, $\eta_c$={packing_fract_coll:.2}, q={q}"
# ax.set_title(ttl)
ylim = [min(ylims),max(ylims)]
xlim = [min(xlims),max(xlims)]
print()
ax.set_ylim(ylim)
ax.set_xlim(xlim)
ax.grid(True,which="both")
ax.legend()
ax.set_ylabel(f"$V'(r)$")
ax.set_xlabel(f"$r/2R$")
fig.tight_layout()
# fig.savefig('analysis/Effective_pot-'+str(n_part)+'.png')
fig.savefig('analysis/Effective_pot-'+'all'+'.png')
# fig.savefig('analysis/Effective_pot-'+'1.2-1.8'+'.png')
