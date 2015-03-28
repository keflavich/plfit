import numpy as np
import plfit
import pylab as plt

xmins = np.logspace(-1,1)
nel = 2000

fig1 = plt.figure(1)
fig1.clf()
ax1 = fig1.add_subplot(2,1,1)
ax2 = fig1.add_subplot(2,1,2)

for alpha in (1.5,2.5,3.5):
    results = []
    for xmin in xmins:
        data = plfit.plexp_inv(np.random.rand(nel), xmin, alpha)
        result = plfit.plfit(data, quiet=True, silent=True)
        results.append(result)

    fitted_xmins = [r._xmin for r in results]
    fitted_alphas = [r._alpha for r in results]
    fitted_alpha_errors = [r._alphaerr for r in results]

    ax1.loglog(xmins, fitted_xmins, 's',
               label='$\\alpha={0}$'.format(alpha),
               alpha=0.5)
    ax1.loglog(xmins, xmins, 'k--', alpha=0.5, zorder=-1)

    ax2.errorbar(xmins, fitted_alphas, yerr=fitted_alpha_errors, marker='s',
               label='$\\alpha={0}$'.format(alpha),
               alpha=0.5)
    ax2.set_xscale('log')

ax1.legend(loc='best')
ax2.set_xlabel("$X_{min}$ (input)")
ax1.set_ylabel("$X_{min}$ (measured)")
ax2.set_ylabel("$\hat{\\alpha}$ (measured)")

