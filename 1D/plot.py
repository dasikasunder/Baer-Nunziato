import numpy as np
import matplotlib.pyplot as plt 

W = np.loadtxt('sol.csv', skiprows=1, delimiter=',')

# Density of solid phase
plt.figure(1)
plt.plot(W[:,0], W[:,1], linestyle='-', linewidth =0.3, color='black', marker = 'o', mew=1.0, ms=1.0)
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho_s$")
plt.savefig("RHO_s.pdf")

# Density of gas phase
plt.figure(2)
plt.plot(W[:,0], W[:,4], linestyle='-', linewidth =0.3, color='black', marker = 'o', mew=1.0, ms=1.0)
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho_g$")
plt.savefig("RHO_g.pdf")


"""
plt.figure(1)
plt.plot(W[:,0], W[:,3], linestyle='-', linewidth =0.3, marker = 'o', mew=1.0, ms=1.0, label='solid-phase')
plt.plot(W[:,0], W[:,6], linestyle='-', linewidth =0.3, marker = 'o', mew=1.0, ms=1.0, label='gas-phase')
plt.xlabel(r"$x$")
plt.ylabel(r"$p$")
plt.legend()
plt.savefig("P.pdf")

plt.figure(2)
plt.plot(W[:,0], W[:,2], linestyle='-', linewidth =0.3, marker = 'o', mew=1.0, ms=1.0, label='solid-phase')
plt.plot(W[:,0], W[:,5], linestyle='-', linewidth =0.3, marker = 'o', mew=1.0, ms=1.0, label='gas-phase')
plt.xlabel(r"$x$")
plt.ylabel(r"$u$")
plt.legend()
plt.savefig("U.pdf")
"""
