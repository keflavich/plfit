import numpy as np
import plfit
import pylab as plt
import itertools

nel = 2000
alpha = 2.5
xmin = 1.0
data = plfit.plexp_inv(np.random.rand(nel), xmin, alpha)

result_py = plfit.plfit(data, quiet=False, silent=False, usecy=False, usefortran=False)
result_cy = plfit.plfit(data, quiet=False, silent=False, usecy=True,  usefortran=False)
result_fo = plfit.plfit(data, quiet=False, silent=False, usecy=False, usefortran=True )

for aa,bb in itertools.combinations((result_py, result_cy, result_fo), 2):
    np.testing.assert_almost_equal(aa._alpha, bb._alpha)
    np.testing.assert_array_almost_equal(aa._xmin_kstest, bb._xmin_kstest)
