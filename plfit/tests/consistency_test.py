from __future__ import print_function
import numpy as np
import plfit
import pylab as plt
import itertools

for ii in range(10):
    nel = 2000
    alpha = 2.5
    xmin = 1.0
    data = plfit.plexp_inv(np.random.rand(nel), xmin, alpha)

    result_py = plfit.plfit(data, quiet=False, silent=False, usecy=False, usefortran=False)
    result_cy = plfit.plfit(data, quiet=False, silent=False, usecy=True,  usefortran=False)
    result_fo = plfit.plfit(data, quiet=False, silent=False, usecy=False, usefortran=True )
    result_py.name = 'python'
    result_cy.name = 'cython'
    result_fo.name = 'fortran'

    for aa,bb in itertools.combinations((result_py, result_cy, result_fo), 2):
        np.testing.assert_almost_equal(aa._alpha, bb._alpha, 5)

        assert aa._ngtx == bb._ngtx
        # should be the same value and exact
        assert aa._xmin == bb._xmin

        maxdiff_xmin_kstest = np.max(np.abs((aa._xmin_kstest[:-1] - bb._xmin_kstest[:-1])))
        maxdiff_alpha_values = np.max(np.abs((aa._alpha_values[:-1] - bb._alpha_values[:-1])))
        print("comparing {0} to {1}".format(aa.name, bb.name))
        print("maxdiff xmin: ", maxdiff_xmin_kstest)
        print("maxdiff alpha: ", maxdiff_alpha_values)

        np.testing.assert_array_almost_equal(aa._xmin_kstest[:-1], bb._xmin_kstest[:-1], 5)

        np.testing.assert_array_almost_equal(aa._alpha_values[:-1], bb._alpha_values[:-1], 5)

