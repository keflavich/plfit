# intended to implement a power-law fitting routine as specified in.....
# http://www.santafe.edu/~aaronc/powerlaws/
#
# The MLE for the power-law alpha is very easy to derive given knowledge
# of the lowest value at which a power law holds, but that point is 
# difficult to derive and must be acquired iteratively.

"""
plfit.py - a python power-law fitter based on code by Aaron Clauset
http://www.santafe.edu/~aaronc/powerlaws/
http://arxiv.org/abs/0706.1062 "Power-law distributions in empirical data" 
Requires pylab (matplotlib), which requires numpy

example use:
from plfit import plfit

MyPL = plfit(mydata)
MyPL.plotpdf(log=True)


"""

import time
import random
import math

class plfit:
    """
    A Python implementation of the Matlab code http://www.santafe.edu/~aaronc/powerlaws/plfit.m
    from http://www.santafe.edu/~aaronc/powerlaws/

    See A. Clauset, C.R. Shalizi, and M.E.J. Newman, "Power-law distributions
    in empirical data" SIAM Review, 51, 661-703 (2009). (arXiv:0706.1062)
    http://arxiv.org/abs/0706.1062

    The output "alpha" is defined such that p(x) ~ (x/xmin)^-alpha
    """

    def __init__(self,x,**kwargs):
        """
        Initializes and fits the power law.  Can pass "quiet" to turn off 
        output (except for warnings; "silent" turns off warnings)
        """
        neg = [i<0 for i in x]
        if any(neg) > 0:
            print "Removed %i negative points" % (sum(neg))
            x = [i for i in x if i > 0]
        self.data = x
        self.plfit(**kwargs)


    def alpha_(self,x):
        def alpha(xmin,x=x):
            """
            given a sorted data set and a minimum, returns power law MLE fit
            data is passed as a keyword parameter so that it can be vectorized
            """
            x = [i for i in x if i>=xmin]
            n = sum(x)
            divsum = sum([math.log(i/xmin) for i in x])
            if divsum == 0: return float('inf')
            a = float(n) / divsum
            return a
        return alpha

    def kstest_(self,x):
        def kstest(xmin,x=x):
            """
            given a sorted data set and a minimum, returns power law MLE ks-test w/data
            data is passed as a keyword parameter so that it can be vectorized

            The returned value is the "D" parameter in the ks test...
            """
            x = [i for i in x if i>=xmin]
            n = float(len(x))
            if n == 0: return float('inf')
            divsum = sum([math.log(i/xmin) for i in x])
            if divsum == 0: return float('inf')
            a = float(n) / divsum
            cx = [float(i)/float(n) for i in xrange(n)]
            cf = [1-(xmin/i)**a for i in x]
            ks = max([abs(a-b) for a,b in zip(cf,cx)])
            return ks
        return kstest
    

    def plfit(self,nosmall=True,finite=False,quiet=False,silent=False,usefortran=False,usecy=False,
            xmin=None, verbose=False):
        """
        A Python implementation of the Matlab code http://www.santafe.edu/~aaronc/powerlaws/plfit.m
        from http://www.santafe.edu/~aaronc/powerlaws/

        See A. Clauset, C.R. Shalizi, and M.E.J. Newman, "Power-law distributions
        in empirical data" SIAM Review, 51, 661-703 (2009). (arXiv:0706.1062)
        http://arxiv.org/abs/0706.1062

        nosmall is on by default; it rejects low s/n points
        can specify xmin to skip xmin estimation

        There are 3 implementations of xmin estimation.  The fortran version is fastest, the C (cython)
        version is ~10% slower, and the python version is ~3x slower than the fortran version.
        Also, the cython code suffers ~2% numerical error relative to the fortran and python for unknown
        reasons.
        """
        x = self.data
        z = sorted(x)
        t = time.time()
        xmins = sorted(set(z))
        argxmins = [z.index(i) for i in xmins]
        self._nunique = len(xmins)
        if xmin is None:
            av  = map(self.alpha_(z),xmins)
            dat = map(self.kstest_(z),xmins)
            sigma = [(a-1)/math.sqrt(len(z)-i+1) for a,i in zip(av,argxmins)]
            if nosmall:
                # test to make sure the number of data points is high enough
                # to provide a reasonable s/n on the computed alpha
                goodvals = [s<0.1 for s in sigma]
                if False in goodvals: 
                    nmax = goodvals.index(False)
                    dat = dat[:nmax]
                    xmins = xmins[:nmax]
                    av = av[:nmax]
                else:
                    print "Not enough data left after flagging - using all positive data."
            if not quiet: print "PYTHON plfit executed in %f seconds" % (time.time()-t)
            self._av = av
            self._xmin_kstest = dat
            self._sigma = sigma
            xmin  = xmins[dat.index(min(dat))] 
        z     = [i for i in z if i >= xmin]
        n     = len(z)
        alpha = 1 + n / sum( [math.log(a/xmin) for a in z] ) 
        if finite:
            alpha = alpha*(n-1.)/n+1./n
        if n == 1 and not silent:
            print "Failure: only 1 point kept.  Probably not a power-law distribution."
            self._alpha = 0
            self._alphaerr = 0
            self._likelihood = 0
            self._ks = 0
            self._ks_prob = 0
            self._xmin = xmin
            return xmin,0
        if n < 50 and not finite and not silent:
            print '(PLFIT) Warning: finite-size bias may be present. n=%i' % n
        # ks = max(abs( numpy.arange(n)/float(n) - (1-(xmin/z)**(alpha-1)) ))
        ks = max( [abs( i/float(n) - (1-(xmin/b)**(alpha-1))) for i,b in zip(xrange(n),z)] )
        # Parallels Eqn 3.5 in Clauset et al 2009, but zeta(alpha, xmin) = (alpha-1)/xmin.  Really is Eqn B3 in paper.
        #L = n*log((alpha-1)/xmin) - alpha*sum(log(z/xmin))
        sl = sum([math.log(a/xmin) for a in z])
        L = (n*math.log((alpha-1)/xmin) - alpha*sl)
        #requires another map... Larr = arange(len(unique(x))) * log((av-1)/unique(x)) - av*sum
        self._likelihood = L
        self._xmin = xmin
        self._xmins = xmins
        self._alpha= alpha
        self._alphaerr = (alpha-1)/math.sqrt(n)
        self._ks = ks  # this ks statistic may not have the same value as min(dat) because of unique()
        #if scipyOK: self._ks_prob = scipy.stats.kstwobign.sf(ks*numpy.sqrt(n))
        self._ngtx = n
        if math.isnan(L) or math.isnan(xmin) or math.isnan(alpha):
            raise ValueError("plfit failed; returned a nan")

        if not quiet:
            if verbose: print "The lowest value included in the power-law fit, ",
            print "xmin: %g" % xmin,
            if verbose: print "\nThe number of values above xmin, ",
            print "n(>xmin): %i" % n,
            if verbose: print "\nThe derived power-law alpha (p(x)~x^-alpha) with MLE-derived error, ",
            print "alpha: %g +/- %g  " % (alpha,self._alphaerr), 
            if verbose: print "\nThe log of the Likelihood (the maximized parameter), ",
            print "Log-Likelihood: %g  " % L,
            if verbose: print "\nThe KS-test statistic between the best-fit power-law and the data, ",
            print "ks: %g" % (ks)

        return xmin,alpha


def plexp(x,xm=1,a=2.5):
    """
    CDF(x) for the piecewise distribution exponential x<xmin, powerlaw x>=xmin
    This is the CDF version of the distributions drawn in fig 3.4a of Clauset et al.
    """

    C = 1/(-xm/(1 - a) - xm/a + math.exp(a)*xm/a)
    Ppl = lambda(X): 1+C*(xm/(1-a)*(X/xm)**(1-a))
    Pexp = lambda(X): C*xm/a*math.exp(a)-C*(xm/a)*math.exp(-a*(X/xm-1))
    d=Ppl(x)
    d[x<xm]=Pexp(x)
    return d

def plexp_inv(P,xm,a):
    """
    Inverse CDF for a piecewise PDF as defined in eqn. 3.10
    of Clauset et al.  
    """

    C = 1/(-xm/(1 - a) - xm/a + math.exp(a)*xm/a)
    Pxm = 1+C*(xm/(1-a))
    pp = P
    x = xm*(pp-1)*(1-a)/(C*xm)**(1/(1-a)) if pp >= Pxm else (math.log( ((C*xm/a)*math.exp(a)-pp)/(C*xm/a)) - a) * (-xm/a)
    #x[P>=Pxm] = xm*( (P[P>=Pxm]-1) * (1-a)/(C*xm) )**(1/(1-a)) # powerlaw
    #x[P<Pxm] = (math.log( (C*xm/a*math.exp(a)-P[P<Pxm])/(C*xm/a) ) - a) * (-xm/a) # exp

    return x

def pl_inv(P,xm,a):
    """ 
    Inverse CDF for a pure power-law
    """
    
    x = (1-P)**(1/(1-a)) * xm
    return x

def test_fitter(xmin=1.0, alpha=2.5, niter=500, npts=1000, invcdf=plexp_inv,
        quiet=True, silent=True):
    """
    Tests the power-law fitter 

    Example (fig 3.4b in Clauset et al.):
    xminin=[0.25,0.5,0.75,1,1.5,2,5,10,50,100]
    xmarr,af,ksv,nxarr = plfit.test_fitter(xmin=xminin,niter=1,npts=50000)
    loglog(xminin,xmarr.squeeze(),'x')

    Example 2:
    xminin=[0.25,0.5,0.75,1,1.5,2,5,10,50,100]
    xmarr,af,ksv,nxarr = plfit.test_fitter(xmin=xminin,niter=10,npts=1000)
    loglog(xminin,xmarr.mean(axis=0),'x')

    Example 3:
    xmarr,af,ksv,nxarr = plfit.test_fitter(xmin=1.0,niter=1000,npts=1000)
    hist(xmarr.squeeze());
    # Test results:
    # mean(xmarr) = 0.70, median(xmarr)=0.65 std(xmarr)=0.20
    # mean(af) = 2.51 median(af) = 2.49  std(af)=0.14
    # biased distribution; far from correct value of xmin but close to correct alpha
    
    Example 4:
    xmarr,af,ksv,nxarr = plfit.test_fitter(xmin=1.0,niter=1000,npts=1000,invcdf=pl_inv)
    print("mean(xmarr): %0.2f median(xmarr): %0.2f std(xmarr): %0.2f" % (mean(xmarr),median(xmarr),std(xmarr)))
    print("mean(af): %0.2f median(af): %0.2f std(af): %0.2f" % (mean(af),median(af),std(af)))
    # mean(xmarr): 1.19 median(xmarr): 1.03 std(xmarr): 0.35
    # mean(af): 2.51 median(af): 2.50 std(af): 0.07

    """
    sz = niter
    xmarr,alphaf_v,ksv,nxarr = ([0]*sz,)*4
    for i in xrange(niter):
        randarr = [random.random() for k in xrange(npts)]
        fakedata = [invcdf(r,xmin,alpha) for r in randarr]
        TEST = plfit(fakedata,quiet=quiet,silent=silent,nosmall=True)
        alphaf_v[i] = TEST._alpha
        ksv[i] = TEST._ks
        nxarr[i] = TEST._ngtx
        xmarr[i] = TEST._xmin

    return xmarr,alphaf_v,ksv,nxarr




