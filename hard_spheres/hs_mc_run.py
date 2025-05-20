import os
from f2py_jit import jit
from config_io_module import read_cnf_atoms, write_cnf_atoms
import numpy as np

# # Compile the relevant modules for hard spheres
# f90 = jit(['maths_module.f90', 'mc_hs_module.f90', 'mc_hs_driver.f90'],
#            #'config_io_module.f90', 'maths_module.f90', 'initialize_module.f90', 'initialize.f90'],
#           flags='-O3 -ffast-math')
# Compile the relevant modules for hard spheres
f90_init = jit(['maths_module.f90', 
           'config_io_module.f90',
           'initialize_module.f90', 
           'initialize_driver.f90', 
        #    'mc_hs_module.f90', 
        #    'mc_hs_driver.f90'
           ],
          flags='-O3 -ffast-math')
f90_mc_hs = jit(['maths_module.f90', 
        #    'config_io_module.f90',
        #    'initialize_module.f90', 
        #    'initialize_driver.f90', 
           'mc_hs_module.f90', 
           'mc_hs_driver.f90',
           'gdr_module.f90'
           ],
          flags='-fbounds-check')
        #   flags='-O3 -ffast-math')


def initialize(density=0.1, nc = 3):
    f90_init.initialize_driver.nc = nc
    f90_init.initialize_driver.density = density
    f90_init.initialize_driver.initialize()

def run_nvt_hs(input_file='./cnf.inp', nsteps=1, dr_max=0.1):
    n, box, r = read_cnf_atoms(input_file)
    f90_mc_hs.mc_module.r = r.transpose()
    f90_mc_hs.mc_driver.box = box
    f90_mc_hs.mc_driver.dr_max = dr_max
    for step in range(nsteps):
        f90_mc_hs.mc_driver.run(10000)

def stamp_pos(input_file='./cnf.inp', output_file='./cnf.out'):
    n, box, r = read_cnf_atoms(input_file)
    write_cnf_atoms(output_file, n, box, f90_mc_hs.mc_module.r.transpose())
    print("Output file: " + "\"" + output_file + "\"")

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

def compute_gdr(filename='./cnf.out'):
    print("Computing g(r): ")
    n, box, r = read_cnf_atoms(filename)
    rho = n/box**3
    print(r.shape)
    # dr = 0.5
    dr = 0.04
    n_dim = int(box / dr)
    hist = np.zeros(n_dim, dtype=np.int64)
    bins = np.zeros(n_dim, dtype=np.float64)
    box = box * np.ones(r.shape[1])
    f90_mc_hs.gdr_module.gr_self(positions=r.transpose(),box=box,dr=dr,hist=hist,bins=bins)
    hist = [hist[i] / (4*np.pi*rho*bins[i]**2) for i in range(n_dim)] # normalize the hist
    # return dr, hist, bins

    import matplotlib.pyplot as plt
    
    # plt.bar(bins, hist, width=dr, align='center', edgecolor='black')
    # plt.title("Histogram of g(r)")
    plt.plot(bins,hist)
    plt.title("Plot of g(r)")
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.savefig("gdr.png")

    print("Histogram saved in gdr.png")



if __name__ == '__main__':
    import argh
    argh.dispatch_command(initialize)
    argh.dispatch_command(run_nvt_hs)
    argh.dispatch_command(stamp_pos)
    # argh.dispatch_command(visualize_3dmol(read_cnf_atoms("cnf.out"), cell=[2, 2, 2]))
    argh.dispatch_command(compute_gdr(filename="../0-1-lab/hard_spheres/cnf.out"))
    # argh.dispatch_command(compute_gdr)

