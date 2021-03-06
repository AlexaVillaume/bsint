
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import bsint

# Integration Time for the 3-body system
t0 = 0
t1 = 1e6
eccentricity = 0.03
num_bodies = 1

def create_initial_values(x_init):
    Y0 = np.zeros((4*num_bodies))
    Y0[0] = x_init # x-position
    Y0[1] = 0.0    # x-velocity
    Y0[2] = 0.0    # y-position
    Y0[3] = 0.0    # y-velocity

    return Y0

def compute_energy(Y0):
    """
    The energy of the system is described by the kinetic energy term and the
    potential for a non-axismmetric potential in cartesian coordinates.
    """

    energy = 0.5*(Y0[1]**2. + Y0[3]**2.) + np.log(np.sqrt(Y0[0]**2 + Y0[2]**2)) + \
            np.log(1 - (eccentricity*Y0[0])/(np.sqrt(Y0[0]**2 + Y0[2]**2)))

    return energy

def compute_vy(x_init):
    """
    This is derived from setting the energy of the system to be a certain value.

    It'd be better to leave this general and have it compute from the inital values
    vector and pass energy as an argument.
    """

    Y0 = create_initial_values(x_init)
    Y0[3] = (2*(-1 - 0.5*np.log(Y0[0]**2) - np.log(1 - (eccentricity))))**(1/2.)

    return Y0

def derivs(t, Y0):
    x_acc = - ((Y0[0] - eccentricity*(Y0[0]**2 + Y0[2]**2)**(1/2))/((Y0[0]**2 + Y0[2]**2) - \
                eccentricity*Y0[0]*(Y0[0]**2 + Y0[2]**2)**(1/2.)))
    y_acc = - ((Y0[2])/((Y0[0]**2 + Y0[2]**2) - \
                eccentricity*Y0[0]*(Y0[0]**2 + Y0[2]**2)**(1/2.)))

    return np.array([Y0[1], x_acc, Y0[3], y_acc])

def create_sos(x_initials):
    """
    Create a surface of section plot for each initial x-position.

    bsint will return a nD array of:
    x-position, x-velocity, y-position, y-velocity
    """

    initial_conds = compute_vy(x_initials)
    out, tout = bsint.bsintegrate(derivs, initial_conds, t0, t1, tacc=1e-14, mxstep=20000)
    pick = []
    for i in range(len(out[:,2])):
        if out[:,2][i] > 1.e-3 and out[:,2][i-1] < 1.e-3:
            pick.append(i)
    #pick = np.where(np.abs(out[:,2]) < 1.e-3)
    print out[pick,0]
    plt.plot(out[pick,0].T, out[pick,1].T, ls='none', marker='.', markersize=2)

x_initials = np.linspace(0.01, 0.6, 20)
for value in x_initials:
    create_sos(value)

plt.ylabel('x-velocity', fontsize=16)
plt.xlabel('x-position', fontsize=16)
plt.show()

