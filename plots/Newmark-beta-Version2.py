###########################################################
#
# 1-DOF spring-mass-damper model with Newmark-beta method
#
# Displacement is found first, and then velocity
# and acceleration are calculated.
#
# Author: Dr Chennakesava Kadapa
# Date: 29-07-2024
#
###########################################################

import matplotlib.pyplot as plt
import numpy as np
from pylab import *


## Input parameters
#
m  = 1.0;                        # mass
k  = 4.0*np.pi**2;               # stiffness
c  = 0.0;                        # damping

d0 = 0.0;                        # initial displacement
v0 = 5.0;                        # initial velocity
a0 = (-c*v0 - k*d0)/m;           # initial acceleration

## Intermediate (auxiliary) parameters
#
w  = np.sqrt(k/m);               # (circular) natural frequency
xi = c/(2.0*np.sqrt(k*m));       # damping ratio
wd = w*np.sqrt(1.0-xi*xi);       # damped natural frequency
T  = 2.0*np.pi/w;                # time period


B1 = d0;
B2 = (v0+xi*w*d0)/wd;

B = np.sqrt(B1*B1+B2*B2);
phi = np.arctan(B2/B1);


dt = T/40.0

tt = np.arange(0.0, 10*T+dt, dt)

Nt = size(tt)
print(Nt)

## parameters for second-order accuracy and unconditonal stability
gamma = 0.5
beta  = 0.25

# beta = 0 for the explicit scheme but this version does not work due to division by beta.

# solution phase starts from here
dispExct = np.zeros((Nt,1), dtype=float)
dispNum  = np.zeros((Nt,1), dtype=float)
veloNum  = np.zeros((Nt,1), dtype=float)


# initialise variables used in the solution
dispPrev = d0
veloPrev = v0
accePrev = a0

disp = d0
velo = v0
acce = a0

# store the solutions at t=0
dispNum[0] = d0;
veloNum[0] = v0;
dispExct[0] = np.exp(-xi*w*tt[0])*B*np.cos(wd*tt[0]-phi);

Keff = m + (gamma*dt)* c + (beta*dt*dt)* k

for ii in range(1,Nt):
    dispExct[ii] = np.exp(-xi*w*tt[ii])*B*np.cos(wd*tt[ii]-phi)

    force = 0.0
    Feff = force + m*( dispPrev + dt*veloPrev + (0.5-beta)*dt*dt*accePrev )
    Feff = Feff  + c*(gamma*dt*dispPrev - (beta-gamma)*dt*dt*veloPrev - (beta-0.5*gamma)*dt*dt*dt*accePrev)

    disp = Feff/Keff
    acce =   (1/beta/dt**2)*(disp-dispPrev) -    (1/beta/dt)*veloPrev - (1/2/beta-1)*accePrev
    velo = veloPrev + (1.0-gamma)*dt*accePrev + gamma*dt*acce

    dispNum[ii] = disp;
    veloNum[ii] = velo;

    # store solution variables
    dispPrev = disp;
    veloPrev = velo;
    accePrev = acce;
    
plot(tt, dispExct,'k')
plot(tt, dispNum, 'b')
plt.show()







