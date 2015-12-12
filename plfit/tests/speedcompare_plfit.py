#import cplfit
import plfit
import time
from numpy.random import rand,seed
from numpy import unique,sort,array,asarray,log,sum,min,max,argmin,argmax,arange
import sys
import powerlaw
from agpy import readcol

try:
    ne = int(sys.argv[1])
    seed(1)
    X=plfit.plexp_inv(rand(ne),1,2.5)
    X[:100] = X[100:200]
except ValueError:
    X = readcol(sys.argv[1])

if len(sys.argv)>2:
    discrete = bool(sys.argv[2])
else:
    discrete=None

print("Cython")
t1=time.time(); p3=plfit.plfit(X,discrete=discrete,usefortran=False,usecy=True); print(time.time()-t1)
print("Fortran")
t1=time.time(); p1=plfit.plfit(X,discrete=discrete,usefortran=True); print(time.time()-t1)
print("Numpy")
t1=time.time(); p3=plfit.plfit(X,discrete=discrete,usefortran=False); print(time.time()-t1)

print("Jeff Alcott's Powerlaw")
t5=time.time(); p5=powerlaw.Fit(X,discrete=discrete); print(time.time()-t5)



print("Pure Python")
t4=time.time(); p4=plfit.plfit_py(X.tolist()); print(time.time()-t4)

