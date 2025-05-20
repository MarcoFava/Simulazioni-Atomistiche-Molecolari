import numpy as np
import matplotlib.pyplot as plt

# P_units = 'kT/sigma**3'
# Phi_units = 'pi/6*rho*sigma**3'

P = []
Rho = []
with open("hard_spheres/EOS.txt", "r") as f:
    data = []
    for line in f:
            line = line.split()
            if line:
                line = [float(i) for i in line]
                Rho.append(line[0])
                P.append(line[1])

plt.plot(Rho,P,'-*')
# plt.xlabel(r'$\Phi$ [$\sigma^{-3}$]')
# plt.ylabel(r'P [kT/$\sigma^3$]')
# plt.legend()
plt.grid()
plt.show()
