#!/usr/bin/python

import math
import random
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Import ndfit
import ndfit as ndfC
import multiprocessing as mp

## # Define some fit functions
def fitfunc(dat,p): return p[2]-(p[0]**2)/(p[0]**2 +(dat[0]-p[1])**2)
def errfunc(dat,p): return fitfunc(dat,p)-dat[1]

def worker(fitfunc,errfunc,data,params,step,i):
    a = datetime.now()
    ndfC.convergence(0.04)
    ndfC.maxdepth(1000)
    ndfC.throttle_factor(60.0)
    ndf = ndfC.run(fitfunc,errfunc,data,params,step,mode="short",throttle=True)
    b = datetime.now()
    #print b-a
    return (ndf,i)

result_list = []
def log_result(result):
    result_list.append(result)

def apply_async_with_callback(data,params,step):
    pool = mp.Pool()
    for i in range(16):
        pool.apply_async(worker,args=(fitfunc,errfunc,data,params,step,i),callback=log_result)
    pool.close()
    pool.join()
    print(result_list)

if __name__ == "__main__":
    # Generate random data
    x = list(np.linspace(0,10,1000));
    y = [fitfunc([i],[2,5,6])+float(random.randint(-5,5))/50 for i in x]

    # Try to fit it
    print "2D test"
    data,params,step = zip(x,y),[4.,5.,4.],[0.01,0.01,0.01]
    print "NDFIT ----------> C VERSION"
    print "---------------------------"
    apply_async_with_callback(data,params,step)




# And plot it
#yfit = ndf.buildcurve(zip(x))
#plt.plot(x,y)
#plt.plot(x,yfit,'r')
#plt.show()






