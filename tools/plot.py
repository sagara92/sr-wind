import numpy as np
import matplotlib.pyplot as plt
c = 2.99e10
r, u, d, h = np.loadtxt('./data/out.dat').T
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(r, d / 1e-10, label=r'$\rho / (10^{-10} \mathrm{g/cm^3})$')
ax1.plot(r, h / c**2 - 1, label=r'$\mu / c^2$')
ax1.plot(r, u, label=r'$\Gamma \beta$')
ax1.plot(r[r>1e10], 1e19/r[r>1e10]**2, c='k', ls='--', lw=1.0, label=r'$\propto r^{-2}$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e8, 2e11)
ax1.set_xlabel(r'$r \ [\mathrm{cm}]$')
ax1.legend()
plt.show()
