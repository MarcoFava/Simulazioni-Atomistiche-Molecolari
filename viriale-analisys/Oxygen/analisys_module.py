import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

########################################################
def read_data(filename="data.txt"):
    data = []
    with open(filename,"r") as f:
        for line in f:
            if "Temperature" in line or "liquid" in line or "infinite" in line or "undefined" in line:
                continue
            line = line.split()
            data.append(line[:-1])
    data = np.array(data,dtype=np.float64)
    return data

########################################################
def write_data(filename,data,first_row=None):
    '''
    Data must be stored in an numpy array column-shaped
    Describe the type of data in "first_row"
    '''
    with open(filename,"w") as f:
        if first_row != None:
            f.write(first_row + "\n")
        for k in range(np.shape(data)[0]):
            for i in range(np.shape(data)[1]):
                f.write(f'{data[k,i]:13.4e}')
            f.write('\n')

########################################################
def lin_regress(X,Y,grade=None):

    if grade==None:
        def f(x,m,q):
            return x*m + q
        fit = curve_fit(f,X,Y,p0=[0,1])
        fit = fit[0]
        m = fit[0]
        q = fit[1]

    else:
        fit = np.polyfit(X,Y,grade)
        m = fit[-2]
        q = fit[-1]

    return m, q

########################################################
def plot_EOS(X,Y,ttl="EOS",lgd="",grade=3):

    m, q = lin_regress(X,Y,grade=grade)
    x = np.linspace(min(X),max(X))

    lgd_dat = "dati " + lgd
    lgd_reg = "regressione " + lgd + r'B$_2$=' + '{:.2e}'.format(m) + r' [m$^3$/mol]'

    plt.figure()
    plt.plot(X,Y,'.-',label=lgd_dat)
    plt.plot(x,m*x+q,'-',label=lgd_reg)
    plt.xlabel(r'$\rho$ [mol/m$^3$]')
    plt.ylabel(r'P / (RT$\rho)$')
    # plt.ylabel(r'$\frac{P}{RT\rho}$')
    plt.title(ttl)
    plt.legend()
    plt.grid()
    plt.show()
    
    return

########################################################
def compute_B2(temperature=300,element='.',grade=None):
    filename = 'data_T' + str(temperature) + '.txt'
    filepath = element + '/data/'
    filename = filepath + filename

    data = read_data(filename=filename)
    T = data[:,0]
    P = data[:,1]
    d = data[:,2]
    red_P = P/8.314/T/d*1e6

    B2, q = lin_regress(d,red_P,grade=grade)

    return B2

########################################################
def virial_coefficient(potential, kT, *args, **kwargs):
    """
    Compute the second virial coefficient for an arbitrary pair potential `u(r)`
    at a temperature equal to `kT` over k_B

    The additional arguments `*args` and `**kwargs` are passed to the
    potential function.
    """
    def integrand(x):
        return - 2*np.pi * x**2 * (np.exp(- potential(x, *args, **kwargs) / kT) - 1)
    return integrate.quad(integrand, 0., np.inf)[0]

def lj(r, epsilon=1.0, sigma=1.0):
    """Lennard-Jones potential"""
    return 4*epsilon * ((sigma/r)**12 - (sigma/r)**6)

@np.vectorize
def sw(r, epsilon=1.0, sigma=1.0, R=1.0):
    """Square well potential"""
    if r <= sigma:
        return 1e10 # np.inf
    elif sigma < r < sigma * R:
        return - epsilon
    else:
        return 0

def virial_coefficient_sw(kT, epsilon=1.0, sigma=1.0, R=1.2):
    """Analytical expression of B_2 for the square well potential"""
    return 2 * np.pi / 3 * sigma**3 * (1 - (np.exp(epsilon / kT) - 1) * (R**3 - 1))

########################################################
def logx_grid(x1, x2, n):
    """Create a list of n numbers exponentially spaced from x1>0 to x2."""
    xx = (x2 / x1)**(1 / n)
    return [x1] + [x1 * xx**(i + 1) for i in range(1, n)]

########################################################
def virial_coefficient_vdw(T, a, b):
    return b - a/T

def virial_coefficient_lj(T, epsilon, sigma):
    return np.array([virial_coefficient(lj, _T, epsilon=epsilon, sigma=sigma) for _T in T])
