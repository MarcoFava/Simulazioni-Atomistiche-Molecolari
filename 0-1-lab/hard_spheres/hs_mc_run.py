import os
from f2py_jit import jit
from config_io_module import read_cnf_atoms, write_cnf_atoms

# Compile the relevant modules for hard spheres
f90 = jit(['maths_module.f90', 'mc_hs_module.f90', 'mc_hs_driver.f90'],
          flags='-O3 -ffast-math')

def run_nvt_hs(input_file='./cnf.inp', nsteps=1, dr_max=0.1):
    n, box, r = read_cnf_atoms(input_file)
    f90.mc_module.r = r.transpose()
    f90.mc_module.n = n
    f90.mc_driver.box = box
    f90.mc_driver.dr_max = dr_max
    for step in range(nsteps):
        f90.mc_driver.run(10000)

# Include these lines only if you want to execute the code from the command line
def stamp_pos(input_file='./cnf.inp', output_file='./cnf.out'):
    n, box, r = read_cnf_atoms(input_file)
    write_cnf_atoms(output_file, n, box, f90.mc_module.r.transpose())
    print("Output file: " + "\"" + output_file + "\"")


if __name__ == '__main__':
    import argh
    # argh.dispatch_command(initialize)
    argh.dispatch_command(run_nvt_hs)
    # argh.dispatch_command(stamp_pos)

