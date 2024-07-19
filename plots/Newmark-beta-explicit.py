###################################################
#
# 1-DOF spring-mass-damper model with Newmark-beta method
#
###################################################


import matplotlib.pyplot as plt
import numpy as np
from pylab import *


## Input parameters
#
m  = 1.0;                        # total mass, m=ms+mf
k  = 4.0*np.pi**2;               # stiffness
c  = 0.0;                        # damping coefficient

d0 = 1.0;                        # initial displacement
v0 = 0.0;                        # initial velocity
a0 = (-c*v0 - k*d0)/m;           # initial acceleration


## Intermediate (auxiliary) parameters
#
w  = np.sqrt(k/m);               # (circular) natural frequency
xi = c/(2.0*np.sqrt(k*m));       # damping ratio
wd = w*np.sqrt(1.0-xi*xi);       # damped natural frequency
T  = 2.0*np.pi/w;                # time period


B1 = d0
B2 = (v0+xi*w*d0)/wd

B = np.sqrt(B1*B1+B2*B2)
phi = np.arctan(B2/B1)


dt = T/100.0

tt = np.arange(0.0, 5*T+dt, dt)

Nt = size(tt)
print(Nt)

## parameters 
gamma = 0.5
beta  = 0.25

##
# solution phase starts from here

dispExct = np.zeros((Nt,1), dtype=float)
dispNumOpenFOAM  = np.zeros((Nt,1), dtype=float)
dispNumCorrect  = np.zeros((Nt,1), dtype=float)


# initialise variables used in the solution
dispPrev = d0
veloPrev = v0
accePrev = a0

disp = d0
velo = v0
acce = a0


# store the solutions at t=0
dispNumOpenFOAM[0] = d0
dispExct[0] = np.exp(-xi*w*tt[0])*B*np.cos(wd*tt[0]-phi)


## Current implementation in OpenFOAM
for ii in range(1,Nt):
    
    dispExct[ii] = np.exp(-xi*w*tt[ii])*B*np.cos(wd*tt[ii]-phi)

    force = 0.0

    #Explicit scheme
    acce = (force - c*veloPrev - k*dispPrev)/m
    velo = veloPrev + (1.0-gamma)*dt*accePrev + gamma*dt*acce
    disp = dispPrev + dt*veloPrev + 0.5*dt*dt*( (1.0-2*beta)*accePrev + 2.0*beta*acce )


    dispNumOpenFOAM[ii] = disp

    # store solution variables
    dispPrev = disp
    veloPrev = velo
    accePrev = acce


# initialise variables used in the solution
dispPrev = d0
veloPrev = v0
accePrev = a0

# store the solutions at t=0
dispNumCorrect[0] = d0

## Correct implementation of explicit version of Newmark scheme
# For the explicit scheme, beta=0.0
#
for ii in range(1,Nt):
    
    dispExct[ii] = np.exp(-xi*w*tt[ii])*B*np.cos(wd*tt[ii]-phi)

    force = 0.0

    #Explicit scheme
    dispPred = dispPrev + dt*veloPrev + dt*dt*0.5*accePrev;
    veloPred = veloPrev + dt*(1.0-gamma)*accePrev;

    acce = (force - c*veloPred - k*dispPred)/(m+c*dt*gamma);

    disp = dispPred;
    velo = veloPred + (dt*gamma)*acce;

    dispNumCorrect[ii] = disp

    # store solution variables
    dispPrev = disp
    veloPrev = velo
    accePrev = acce


plt.plot(tt, dispExct,        '-',  color='b', label="Exact",    linewidth=1)
plt.plot(tt, dispNumOpenFOAM, '-',  color='k', label="OpenFOAM", linewidth=1)
plt.plot(tt, dispNumCorrect,  '--', color='r', label="Correct",  linewidth=3)

plt.axis([0, 5, -3, 3])
plt.grid('on')

plt.legend(loc='upper left', markerscale=1.0, ncol=1, handlelength = 3.33, numpoints=1, fontsize=12)
#plt.tight_layout()

plt.show()
outfile = 'graph.png'
plt.savefig(outfile, dpi=500)







