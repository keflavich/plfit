Power-law Distribution Fitting
==============================
This is a python implementation of a power-law distribution fitter.  The code
here was originally hosted on `agpy
<http://code.google.com/p/agpy/source/browse/wiki/PowerLaw.wiki>`_ but was
moved and re-packaged to make setup.py cleaner.  

`API Documentation <http://plfit.readthedocs.org/>`_

Installation
------------

I've attempted to make the setup.py file work nicely, but it includes some hacks, so if you run into trouble,
please report it on `github <github.com/keflavich/plfit>`_::

    git clone git@github.com:keflavich/plfit.git
    cd plfit
    python setup.py install

*If* ``python setup.py install`` doesn't work, you can try the following:

To install the cython function, run:
``python setup.py build_ext --inplace``

To install the fortran function::

    cd plfit/plfit/ 
    f2py -c fplfit.f -m fplfit --fcompiler=gfortran

Description
-----------

Aaron Clauset et al. address the issue of fitting power-laws to distributions
on `this website <http://www.santafe.edu/~aaronc/powerlaws/>`_ and in their paper
`Power-law distributions in empirical
data <http://code.google.com/p/agpy/source/browse/wiki/PowerLaw.wiki>`_.  I have
created a python implementation of their code because I didn't have matlab or R
and wanted to do some power-law fitting. 

Power-laws are very commonly used in astronomy and are typically used to
describe the initial mass function (IMF), the core mass function (CMF), and
often luminosity distributions.  Most distributions in astronomy tend to be
apparent power-laws because the source counts are too few or too narrow to
distinguish powerlaws from log-normal and other distributions.  But, to this
end, I've included the testing mechanism to test for consistency with a power
law as described in the above paper.

The python internal documentation is complete.  A brief description of relevant functions is included here for convenience:

plfit is implemented as a class.  This means that you import plfit, and declare an instance of the plfit class::

    import plfit
    X = rand(1000)
    myplfit = plfit.plfit(X)

The results of the fit are printed to the screen (if desired) and are stored as part of the object.

``alpha_`` and ``kstest_`` are functions used internally to determine the ks-statistic and alpha values as a function of xmin.

There are 3 predefined plotting functions:
  * ``alphavsks`` plots alpha on the y-axis vs. the ks statistic value on the
    x-axis with the 'best-fit' alpha value plotted with error bars.   These
    plots are a useful way to determine if other values of xmin are similarly
    good fits.
  * ``plotcdf`` plots the cumulative distribution function along with the
    best-fit power law
  * ``plotpdf`` plots a histogram of the PDF with the best fit power law.  It
    defaults to log binning (i.e. a linear power-law fit) but can do dN/dS and
    linear binning as well.

Other useful functions:
 * ``test_pl`` uses the fitted power-law as the starting point for a monte-carlo
   test of whether the powerlaw is an acceptable fit.  It returns a "p-value" that
   should be >0.1 if a power-law fit is to be considered (though a high p-value
   does not ensure that the distribution function is a power law!).

 * ``plexp_inv`` creates a cutoff power-law distribution with an exponential
   tail-off.  It is useful for tests.
 * ``pl_inv`` creates a pure cutoff power-law distribution
 * ``test_fitter`` uses the previous two functions to test the fitter's ability
   to return the correct xmin and alpha values for large numbers of iterations


The powerlaw fitter is very effective at returning the correct value of alpha
but not as good at returning the correct value of xmin.

There are 3 implementations of the code internals.  fplfit.f is a fortran
function, cplfit.pyx is a cython function, and plfit.py is the wrapper and
includes a python-only implementation that requires numpy.  FORTRAN is fastest,
follow closely by cython.  Python is ~3x slower.  

As of November 21, 2011, there is a pure python (i.e., no numpy) implementation
at <https://github.com/keflavich/plfit/blob/master/plfit/plfit_py.py> - you can just
put this file in your local working directory and import it, since it contains
no requirements beyond pure python.  It's slower and hobbled, but it works, and perhaps
will run fast with `pypy <http://pypy.org/>`_.


For usage *examples*, see
 * `<https://github.com/keflavich/plfit/blob/master/plfit/tests/clauset2009_tests.py>`_
 * `<https://github.com/keflavich/plfit/blob/master/plfit/tests/plfit_tests.py>`_
 * `<https://github.com/keflavich/plfit/blob/master/plfit/tests/speedcompare_plfit.py>`_

A very simple example::

    import plfit
    from numpy.random import rand,seed

    # generate a power law using the "inverse" power-law generator code
    X=plfit.plexp_inv(rand(1000),1,2.5)

    # use the numpy version to fit (usefortran=False is only needed if you installed the fortran version)
    myplfit=plfit.plfit(X,usefortran=False)
    # output should look something like this:
    # PYTHON plfit executed in 0.201362 seconds
    # xmin: 0.621393 n(>xmin): 263 alpha: 2.39465 +/- 0.0859979   Log-Likelihood: -238.959   ks: 0.0278864 p(ks): 0.986695

    # generate some plots
    from pylab import *
    figure(1)
    myplfit.plotpdf()

    figure(2)
    myplfit.plotcdf()


*If you use this code, please cite Clauset et al 2009 and consider posting a comment below.*  

Direction citations to the source are welcome!  The python translation has been cited in the following works (and perhaps others?):

 * http://adsabs.harvard.edu/abs/2011ApJ...735...51M
 * http://adsabs.harvard.edu/abs/2011ApJ...736..149G
 * http://www.rsc.org/suppdata/CC/c0/c0cc00366b/c0cc00366b.pdf


.. image:: https://d2weczhvl823v0.cloudfront.net/keflavich/plfit/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

