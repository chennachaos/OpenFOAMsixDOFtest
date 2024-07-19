import matplotlib.pyplot as plt
import numpy as np
from pylab import *

params = {'axes.labelsize': 14,
     'legend.fontsize': 14,
     'xtick.labelsize': 14,
     'ytick.labelsize': 14}
     
rcParams.update(params)     


def  processinputfile(filename):
    with open(filename) as f:
        inpdata = f.readlines()[3:]

    #print(inpdata)
    #result1=[item.split("\n")[0] for item in inpdata]

    nrows = shape(inpdata)[0]
    timearray = np.zeros(nrows)
    dataarray = np.zeros( (nrows,3) )

    ind=0
    #for item in inpdata:
    for i in range(nrows):
        item = inpdata[i]
        #print(item)
        result1 = item.split('\n')[0]
        #print(result1)
        result2 = result1.split('\t')
        #print(result2)
        timearray[ind] = float(result2[0])
        #
        result3=result2[1].replace('(','')
        result4=result3.replace(')','')
        dataarray[ind,:] = result4.split(' ')
        ind = ind+1
        #
        dispY = dataarray[:,1]
    return timearray, dispY

def  exact_solution(timearray):
    m  = 1.0
    k  = 4.0*np.pi**2
    c  = 0.0

    d0 = 0.0
    v0 = 5.0
    a0 = (-c*v0 - k*d0)/m

    w  = np.sqrt(k/m)
    xi = c/(2.0*np.sqrt(k*m))
    wd = w*np.sqrt(1.0-xi*xi)

    B1 = d0
    B2 = (v0+xi*w*d0)/wd

    B = np.sqrt(B1*B1+B2*B2)
    phi = np.arctan(B2/B1)
    print(phi)

    Nt = shape(timearray)[0]
    disp_exact = np.zeros(Nt)

    for ii in range(Nt):
    #for t in timearray:
        t = timearray[ii]
        #print(t)
        disp_exact[ii] = np.exp(-xi*w*t)*B*np.cos(wd*t-phi)

    return disp_exact

#filename='../postProcessing/sixDoFRigidBodyState/0/sixDoFRigidBodyState-old.dat'
#file_dt0p002 = 'sixDoFRigidBodyState-CrankNicolson-dt0p002.dat'
#file_dt0p001 = 'sixDoFRigidBodyState-CrankNicolson-dt0p001.dat'

#file_dt0p002 = 'sixDoFRigidBodyState-Newmark-dt0p002.dat'
#file_dt0p001 = 'sixDoFRigidBodyState-Newmark-dt0p001.dat'

file_dt0p002 = 'sixDoFRigidBodyState-Symplectic-dt0p002.dat'
file_dt0p001 = 'sixDoFRigidBodyState-Symplectic-dt0p001.dat'


timearray_dt0p002, dispY_dt0p002 = processinputfile(file_dt0p002)
timearray_dt0p001, dispY_dt0p001 = processinputfile(file_dt0p001)

disp_exact = exact_solution(timearray_dt0p001)

#plt.figure(2)
plt.plot(timearray_dt0p001, disp_exact,    '-', color='r', label="Analytical",        linewidth=2.0, markersize=7.0)
plt.plot(timearray_dt0p002, dispY_dt0p002, '-', color='k', label=r'$\Delta t=0.002$', linewidth=1.0, markersize=7.0)
plt.plot(timearray_dt0p001, dispY_dt0p001, '-', color='b', label=r'$\Delta t=0.001$', linewidth=1.0, markersize=7.0)

plt.xlabel('Time',fontsize=14)
plt.ylabel('Displacement ',fontsize=14)

plt.axis([0, 10, -1.5, 1.5])
plt.grid('on')

plt.legend(loc='upper left', markerscale=1.0, ncol=3, handlelength = 2.0, numpoints=1, fontsize=12)
plt.tight_layout()

plt.show()
outfile = 'graph.png'
plt.savefig(outfile, dpi=500)


