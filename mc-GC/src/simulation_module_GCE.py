import numpy as np

def run_simulation(function,parallel=True,clear_first=False, wait=False,**kwargs):
    import pantarei as rei

    scheduler = rei.scheduler.ThreadScheduler(backend='process')
    task = rei.Task(function,clear_first=clear_first)
    if parallel: job = rei.Thread(task, scheduler=scheduler)
    else: job = task

    data = job(kwargs=kwargs)
    if wait: return 0
    scheduler.wait()
    return data


def plot_hist(data):
    import matplotlib.pyplot as plt
    Ns = data['N']

    plt.hist(Ns)


def weight(eta,beta_0,beta_1=0,mu_0=0,mu_1=0,energy=None):
    if energy is None: energy = np.zeros_like(eta)
    return np.array([np.exp( -(beta_1-beta_0)*energy[j] + (beta_1*mu_1-beta_0*mu_0)*j + beta_0*eta[j] )  for j in range(len(eta))])


def read_eta(eta_namefile='save_past/eta_save.txt'):

    eta = []
    with open(eta_namefile,'r') as f:
        for line in f:
            val = line.split()[1]
            eta.append(float(val))

    return np.array(eta)

    
def save_eta(val,temp_eta_file='save_past/eta_save.txt',bins_divider=20):
    import matplotlib.pyplot as plt

    bins = int((max(val)-min(val))/bins_divider)
    p = plt.hist(val,bins)


    p_x = p[1].copy()
    p_y = p[0].copy()
    for k in range(len(p_x)-1):
        p_x[k] += (p_x[k+1]-p_x[k])/2


    prob = np.zeros(1000)
    for n in range(1000):
        for k in range(len(p_y)-1):
            m = (p_y[k+1]-p_y[k])/(p_x[k+1]-p_x[k])
            if n >= p_x[k] and n < p_x[k+1]:
                prob[n] = (n-p_x[k])*m + p_y[k]
        if prob[n] < min(p_y):
            prob[n] = min(p_y)


    x = np.linspace(0,1000,1000)
    plt.plot(x,prob,'-')
    dx = (max(val) - min(val))/20
    plt.xlim([min(val) - dx, max(val) + dx])

    prob /= prob.sum()
    eta = np.log(prob)

    with open(temp_eta_file,'w') as f:
        for k in range(eta.shape[0]):
            f.write(f'{k:6}   {eta[k]}\n')