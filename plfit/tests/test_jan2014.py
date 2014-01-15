import scipy.io
import plfit
import time

m = scipy.io.loadmat('AUD_Ret_1000.mat')
A = m['A'].squeeze()
Pnpy = plfit.plfit(A)
Pfor = plfit.plfit(A)
Pcy = plfit.plfit(A)
Py = plfit.plfit_py(A)

t0 = time.time()
r_for = Pfor.plfit(usefortran=True, verbose=True, quiet=False, discrete=False, nosmall=False)
t1 = time.time()
r_cy  = Pcy.plfit(usecy=True, verbose=True, quiet=False, discrete=False, nosmall=False)
t2 = time.time()
r_ppy = Py.plfit(nosmall=False)
t3 = time.time()
r_npy = Pnpy.plfit(usefortran=False, usecy=False, verbose=True, quiet=False, discrete=False, nosmall=False)
t4 = time.time()

print r_npy
print r_ppy
print r_cy
print r_for
print

import powerlaw
if 'results' not in locals():
    t5 = time.time()
    results = powerlaw.Fit(A[A>0])
    t6 = time.time()
print t6-t5,"powerlaw: ",results.power_law.alpha, results.power_law.xmin
print t4-t3,"npy: ",r_npy
print t3-t2,"ppy: ",r_ppy
print t2-t1,"for: ",r_cy
print t1-t0,"cy:  ",r_for
print


print results.power_law.alpha, results.power_law.xmin
print "npy: ",[(x1-x2) for x1,x2 in zip(r_npy,(results.power_law.xmin, results.power_law.alpha, ))]
print "ppy: ",[(x1-x2) for x1,x2 in zip(r_ppy,(results.power_law.xmin, results.power_law.alpha, ))]
print "for: ",[(x1-x2) for x1,x2 in zip(r_for,(results.power_law.xmin, results.power_law.alpha, ))]
print "cy:  ",[(x1-x2) for x1,x2 in zip(r_cy ,(results.power_law.xmin, results.power_law.alpha, ))]

from pylab import *

mpl.rc_file('/Users/adam/.matplotlib/ggplotrc')

figure(1)
clf()
plot(results.sigmas[np.isfinite(results.sigmas)], label='powerlaw')
plot(Py._sigma, label='py')
plot(Pfor._sigma, label='for')
plot(Pnpy._sigma, label='npy')
plot(Pcy._sigma, label='cy')
gca().set_ylim(0,5)
legend(loc='best')

figure(2)
clf()
plot(np.array(results.Ds)[np.isfinite(results.alphas)], linewidth=2, label='powerlaw')
plot(Py._xmin_kstest, linewidth=3, alpha=0.6, linestyle='--', label='py')
plot(Pnpy._xmin_kstest, linewidth=3, alpha=0.6, linestyle=':', label='npy')
plot(Pcy._xmin_kstest, linewidth=3, alpha=0.6, linestyle=':', label='cy')
plot(Pfor._xmin_kstest, linewidth=3, alpha=0.6, linestyle=':', label='for')
legend(loc='best')

print [(x,locals()[x]._xmin_kstest[-1]) for x in "Py,Pnpy,Pcy,Pfor".split(',')]

