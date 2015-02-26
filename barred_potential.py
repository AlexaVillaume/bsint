"""
For problem 2 on problem set 3 of Winter 2015 Galaxy Dynamics class.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import bsint

num_bodies = 1      # See pythagorean_three_body.py
t0 = 0
t1 = 1e6            # endpoints for the time array
ecce = np.sqrt(0.01)
R_c = 0.5           # don't know what this represents
V_0 = 1    # don't know what this represents

def create_initial_values(x_init):
    Y0 = np.zeros((4*num_bodies))
    Y0[0] = x_init # x-position
    Y0[1] = 0.0    # x-velocity
    Y0[2] = 0.0    # y-position
    Y0[3] = 0.0    # y-velocity

    return Y0

def galaxy_potential(init_conds):

    r_sqr = init_conds[0]**2 + init_conds[2]**2
    gal_pot = 0.5*V_0*np.log(R_c**2 + r_sqr)*(1 + (ecce**2)*((init_conds[0]**2 - init_conds[2]**2)/r_sqr))

    return gal_pot

def bar_potential(time, init_conds):
    """
    The bar has a seperate influence on the test particle, distince from
    that of the potential of the galaxy.

    See the solution set for winter term 2015 Galaxy Dynamics problem set 3
    for derivation of this term.
    """

    bar_strength = 1.25
    other_thing = (1e-2)**2
    r = np.sqrt(init_conds[0]**2 + init_conds[2]**2)
    # not quite sure if this will work as a seperate function with the time
    # dependency
    bar_pot = other_thing*(((init_conds[0]**2 - \
              init_conds[1]**2)/r)*np.cos(bar_strength*time) - \
              ((2*init_conds[0]*init_conds[2])/r)*(np.sin(bar_strength*time)))

    return bar_pot

def compute_vy(x_init):
    """
    Set the energy of the system to a certain value and compute
    the y velocity given that energy.
    """

    init_conds = create_initial_values(x_init)
    bar_pot = bar_potential(t0, init_conds)
    gal_pot = galaxy_potential(init_conds)
    energy = 5

    init_conds[3] = (2*(energy - bar_pot - gal_pot))**0.5 - init_conds[1]

    return init_conds

def derivs(t, Y0):
    """
    See the solution set for winter term 2015 Galaxy Dynamics problem set 3
    for derivation of this term.
    """

    r_sqr = init_conds[0]**2 + init_conds[2]**2
    r_b = init_conds[0]**2 - init_conds[2]**2
    term1 = R_c**2 + r_sqr

    #bar_pot = bar_potential()

    x_acc = -(V_0**2 * ()/term1 + 0.5*V_0**2 *()*np.log(term1))

    #y_acc =

    #return np.array([Y0[1], x_acc, Y0[3], y_acc])

    return 0

def integrate(x_initial):
    initial_conds = compute_vy(x_initial)
    #out, tout = bsint.bsintegrate(derivs, initial_conds, t0, t1, tacc=1e-14, mxstep=20000)

if __name__ == '__main__':
    x_initial = 1
    integrate(x_initial)


