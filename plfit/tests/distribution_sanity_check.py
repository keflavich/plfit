from __future__ import print_function
import numpy as np
import plfit
import pylab as plt

nel = 2000
alpha = 2.5
xmin = 1.0

data_prob = np.random.rand(nel)
data = plfit.plexp_inv(data_prob, xmin, alpha)
inds = np.argsort(data)
data_prob = data_prob[inds]
data = data[inds]

xdata = np.logspace(-2,2,10000)
plt.figure(1).clf()
plt.loglog(xdata, 1-plfit.plexp_cdf(xdata, xmin=xmin, alpha=alpha), 'k', linewidth=5, alpha=0.2, zorder=-1)
plt.loglog(xdata, 1-plfit.plexp_cdf(xdata, xmin=xmin, alpha=alpha, pl_only=True), 'g', linewidth=3, alpha=1, zorder=0)
plt.loglog(xdata, 1-plfit.plexp_cdf(xdata, xmin=xmin, alpha=alpha, exp_only=True), 'b', linewidth=3, alpha=1, zorder=0)
plt.plot(xmin, 1-plfit.plexp_cdf(xmin, xmin=xmin, alpha=alpha, exp_only=True), 'kx', markersize=20, linewidth=3, alpha=1, zorder=1)
plt.ylim(1e-2,1)

plt.figure(2).clf()
plt.loglog(data, 1-plfit.plexp_cdf(data, xmin=xmin, alpha=alpha), 'k', linewidth=10, alpha=0.1, zorder=-1)
plt.loglog(data, 1-plfit.plexp_cdf(data, xmin=xmin, alpha=alpha, pl_only=True), 'g', linewidth=3, alpha=0.5, zorder=0)
plt.loglog(data, 1-plfit.plexp_cdf(data, xmin=xmin, alpha=alpha, exp_only=True), 'b', linewidth=3, alpha=0.5, zorder=0)
plt.plot(xmin, 1-plfit.plexp_cdf(xmin, xmin=xmin, alpha=alpha, exp_only=True), 'kx', markersize=20, linewidth=3, alpha=1, zorder=1)
plt.loglog(data, 1-data_prob, 'r.', zorder=2, alpha=0.1)
plt.ylim((1-data_prob).min(),1)

result_fo = plfit.plfit(data, quiet=False, silent=False, usecy=False, usefortran=True )
plt.plot(result_fo._xmin, 1-plfit.plexp_cdf(result_fo._xmin, xmin=xmin,
                                            alpha=alpha, exp_only=True), 'ko',
         markerfacecolor='none', markeredgecolor='k', markersize=20,
         linewidth=3, alpha=1, zorder=1)
