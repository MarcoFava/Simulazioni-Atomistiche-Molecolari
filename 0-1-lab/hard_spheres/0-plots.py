import numpy as np
import matplotlib.pyplot as plt

# P_units = 'kT/sigma**3'
# Phi_units = 'pi/6*rho*sigma**3'

Rho = [0.01] + [0.1*(i+1) for i in range(13)]
L = [22.104191,10.259856,8.143253,7.113787,6.463305,6.000000,5.646216,
         5.363421,5.129928,4.932425,4.762203,4.613286,4.481405,4.363418]
Phi = [np.pi/6*Rho[i] for i in range(np.size(Rho))]
P = [0.010191,0.124086,0.308506,0.591740,1.004763,1.603520,
     2.561307,3.955812,6.127113,9.197833,10.148582,
     14.640719,23.563833,47.597309]

plt.plot(Phi,P,'*-',label="P")

# factor Beta*P/Rho
factor = [P[i]/Rho[i] for i in range(np.size(P))]
plt.plot(Phi,factor,'*-',label=r"$\frac{\beta P}{\rho}$")


# theoretical line
x = np.linspace(Phi[0], Phi[-1])
g = 1
P_theo = x + 2/3*np.pi*x**2*1

#    plt.plot(x,P_theo)

plt.xlabel(r'$\Phi$ [$\sigma^{-3}$]')
plt.ylabel(r'P [kT/$\sigma^3$]')
plt.legend()
plt.grid()
plt.show()
