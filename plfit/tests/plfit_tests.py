import numpy as np
import plfit
from plfit import cplfit
from plfit import fplfit
import time
import pylab as plt
#from numpy.random import rand
#from numpy import unique,sort,array,asarray,log,sum,min,max,argmin,argmax,arange
import sys
""" 
Test code for fixed, varying xmin.  Alpha is set to 2.5, and argv[1] random power-law
distributions are then fit.  

Timing tests implemented in speedcompare_plfit.py
X=plfit.plexp_inv(rand(ne),1,2.5)
t1=time.time(); p3=plfit.plfit(X,usefortran=False,usecy=True); print time.time()-t1
t1=time.time(); p1=plfit.plfit(X); print time.time()-t1
t1=time.time(); p3=plfit.plfit(X,usefortran=False); print time.time()-t1
"""

if len(sys.argv) > 1:
    ntests = int(sys.argv[1])
    if len(sys.argv) > 2:
        nel = int(sys.argv[2])
        if len(sys.argv) > 3:
            xmin = float(sys.argv[3])
        else:
            xmin = 0.5
    else: 
        nel = 1000
else:
    nel = 1000
    xmin = 0.5
    ntests = 1000

a = np.zeros(ntests)
for i in range(ntests):
    X=plfit.plexp_inv(np.random.rand(nel),xmin,2.5)
    p=plfit.plfit(X,xmin=xmin,quiet=True,silent=True)
    a[i]=p._alpha

h,b = plt.hist(a,bins=30)[:2]
bx = (b[1:]+b[:-1])/2.0

from agpy import gaussfitter
p,m,pe,chi2 = gaussfitter.onedgaussfit(bx,h,params=[0,ntests/10.0,2.5,0.05],fixed=[1,0,0,0])

fig1 = plt.figure(1)
fig1.clf()
plt.plot(bx,m)

print("XMIN fixed: Alpha = 2.5 (real), %0.3f +/- %0.3f (measured)" % (p[2],p[3]))


a=plt.zeros(ntests)
xm=plt.zeros(ntests)
for i in range(ntests):
    data = plfit.plexp_inv(np.random.rand(nel),xmin,2.5)
    p = plfit.plfit(data,quiet=True,silent=True)
    a[i] = p._alpha
    xm[i] = p._xmin

fig2 = plt.figure(2)
fig2.clf()
h1,b1 = plt.hist(a,bins=30)[:2]
plt.xlabel('alpha')
bx1 = (b1[1:]+b1[:-1])/2.0

p1,m1,pe1,chi21 = gaussfitter.onedgaussfit(bx1,h1,params=[0,ntests/10.0,2.5,0.05],fixed=[1,0,0,0])
plt.plot(bx1,m1)
print("XMIN varies: Alpha = 2.5 (real), %0.3f +/- %0.3f (measured)" % (p1[2],p1[3]))

fig3 = plt.figure(3)
fig3.clf()
h2,b2 = plt.hist(xm,bins=30)[:2]
plt.xlabel('xmin')
bx2 = (b2[1:]+b2[:-1])/2.0

p2,m2,pe2,chi2 = gaussfitter.onedgaussfit(bx2,h2,params=[0,ntests/10.0,xmin,0.2],fixed=[1,0,0,0])
plt.plot(bx2,m2)
print("XMIN varies: XMIN = %0.3f (real), %0.3f +/- %0.3f (measured)" % (xmin,p2[2],p2[3]))
