#!/usr/bin/python

import math
import random
import numpy as np
import matplotlib.pyplot as plt

# Import ndfit
import ndfit as ndf

## Define some fit functions. 
def fitfunc(dat,p,c):

    return c[0]*p[2]+(c[1]*p[0]**2)/(p[0]**2 +(dat[0]-p[1])**2)

def errfunc(dat,p,c): 
    return fitfunc(dat,p,c)-dat[1]

if __name__ == "__main__":
    
    params = [2.0,5.0,6.0]  # Parameters of the fit function
    consts = [1.3,1.5]      # Constants passed to the fit function 

    # Generate randomized data. note how everything passed to the fit function is 
    # passes as a list. 
    x = list(np.linspace(0,10,150));
    y = [fitfunc([i],params , consts)+float(random.randint(-5,5))/50 for i in x]

    # Try to fit it. 
    print("2D test")
    guess = [2.10, 5.80, 5.00]    # <--- Guess of the parameters
    step  = [0.01, 0.01, 0.01] # <--- Step size you want ndfit to take for params

    # ndfit ALWAYS expects a list of tuples (x,y) so zip your lists. 
    data   = list(zip(x,y))

    # Set the ndfit parameters
    ndf.convergence(0.04)     # <--- minimum convergence
    ndf.maxdepth(1000)        # <--- maximum recursion depth (number of iterations) 
    ndf.throttle_factor(60.0) # <--- Dynamic fitting ... higher = faster convergences
    
    # RUN ndfit .... this returns an NDFIT object. 
    NDF = ndf.run(fitfunc, errfunc, data, guess, consts, step, mode="full",throttle=True)

    #print(dir(NDF))         # <--- show the list of things that you have in the NDF object
    print("-------------- NDFIT RESULT IS: -----------------")
    print(NDF.getresult())    # <--- Print the result of the fit (entropy , [resulting parameters])
    curve = NDF.buildcurve(x) # <--- Build the final curve from the original x data 
    print("----------- Compare with %s --------------"%(params))
    
    #plt.figure(1)
    #plt.plot(NDF.getentropy())
    #plt.xlabel("Resursion Depth")
    #plt.ylabel("Fit Entropy")

    plt.figure(2) 
    plt.plot(x,y)     # <--- plot the original data
    plt.plot(x,curve, "k", lw=2) # <--- plot the result of the fit
    plt.show()
    plt.xlabel("x Data")
    plt.ylabel("y Data")

    





