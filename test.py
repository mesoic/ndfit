#!/usr/bin/python

import math
import random
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Import ndfit
import ndfit as ndfC

## # Define some fit functions
def fitfunc(dat,p): return p[2]-(p[0]**2)/(p[0]**2 +(dat[0]-p[1])**2)
def errfunc(dat,p): return fitfunc(dat,p)-dat[1]

# Generate random data
x = list(np.linspace(0,10,1000));
y = [fitfunc([i],[2,5,6])+float(random.randint(-5,5))/50 for i in x]

# Try to fit it
print "2D test"
data,params,step = zip(x,y),[4.,5.,4.],[0.01,0.01,0.01]
print "NDFIT ----------> C VERSION"
print "---------------------------"
a = datetime.now()
ndfC.convergence(0.04)
ndfC.maxdepth(1000)
ndfC.throttle_factor(60.0)
ndf = ndfC.run(fitfunc,errfunc,data,params,step,mode="short",throttle=True)
b = datetime.now()
print b-a
print ndf.getresult()
# And plot it
yfit = ndf.buildcurve(zip(x))
plt.plot(x,y)
plt.plot(x,yfit,'r')
plt.show()






