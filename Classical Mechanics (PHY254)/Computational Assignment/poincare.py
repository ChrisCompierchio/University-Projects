from scipy import interpolate
from scipy import optimize
from numpy import *

import math
from pylab import *


def poincaresection(x,xdot,phi,t,omega,phimod2pi):

    xinterp = interpolate.UnivariateSpline(t, x, s=0, k=3)
    xdotinterp = interpolate.UnivariateSpline(t, xdot, s=0, k=3)
    
    tmax=t[-1]
    
    numzeros=int(floor((tmax*omega)/(2*math.pi)))
        
    phiroot=0
    i=1
    
    xpoincare=zeros(numzeros)
    xdotpoincare=zeros(numzeros)
    
    while True:
        phiroot=2*math.pi*i/omega+phimod2pi/omega
        
        if phiroot>tmax:
            break
            
        xpoincare[i-1]=xinterp(phiroot)
        xdotpoincare[i-1]=xdotinterp(phiroot)
        i=i+1

    cutoff = int(numzeros/2)
    return [xpoincare[cutoff:],xdotpoincare[cutoff:]]
    
