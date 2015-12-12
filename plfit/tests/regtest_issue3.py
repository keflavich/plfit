import os

if not os.path.exists('tst.csv'):
    import requests
    result = requests.get('https://gist.githubusercontent.com/vfilimonov/1072e402e922712ad980/raw/27cc61d65590b382ec39120a1d25d4bd3abcfb4d/tst.csv')
    with open('tst.csv','w') as f:
        f.write(result.content)

import numpy as np
import plfit

y = np.genfromtxt('tst.csv', delimiter=',')
tst_fit_py = plfit.plfit(y, usecy=False, usefortran=False, discrete=False)
tst_fit_fo = plfit.plfit(y, usecy=False, usefortran=True, discrete=False)
tst_fit_cy = plfit.plfit(y, usecy=True, usefortran=False, discrete=False)

print("py: ",tst_fit_py._xmins.shape, tst_fit_py._xmin_kstest.shape)
print("cy: ",tst_fit_cy._xmins.shape, tst_fit_cy._xmin_kstest.shape)
print("fo: ",tst_fit_fo._xmins.shape, tst_fit_fo._xmin_kstest.shape)

def func(xmin):
    ff = plfit.plfit(y, xmin=xmin, quiet=True, silent=True)
    return ff._ks

tst_KS_plfit = [func(xmin) for xmin in tst_fit_py._xmins]

import pylab as pl
pl.plot(tst_fit_py._xmins, tst_KS_plfit, 'g-')
pl.plot(tst_fit_py._xmins, tst_fit_py._xmin_kstest, 'r-', alpha=0.5)
pl.plot(tst_fit_cy._xmins, tst_fit_cy._xmin_kstest, 'b--', linewidth=2, alpha=0.5)
pl.plot(tst_fit_fo._xmins, tst_fit_fo._xmin_kstest, 'k:', linewidth=2, alpha=0.5)
