"""
For problem 2 on problem set 3 of Winter 2015 Galaxy Dynamics class.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import bsint

num_bodies = 1      # See pythagorean_three_body.py
t0 = 0
t1 = 10e2            # endpoints for the time array
ecce = np.sqrt(0.01)
R_c = 0.5           # don't know what this represents
V0_sqr = 1    # don't know what this represents
Omega_b = 1.25

def create_initial_values(x_init):
    Y0 = np.zeros((4*num_bodies))
    Y0[0] = x_init # x-position
    Y0[1] = 0.0    # x-velocity
    Y0[2] = 0.0    # y-position
    Y0[3] = 0.0    # y-velocity

    return Y0

def co_rotation():
    return np.sqrt(V0_sqr - (R_c**2)*(Omega_b**2))/Omega_b

def inner_linblad():
    return -(-np.sqrt(2)*np.sqrt(V0_sqr) + 2*np.sqrt(V0_sqr)/2*Omega_b

def outer_liblad():
    return (np.sqrt(2)*np.sqrt(V0_sqr) + 2*np.sqrt(V0_sqr)/2*Omega_b

def compute_vy(x_init):
    """
    Use the rotation curve to get initial y-velocity
    """

    init_conds = create_initial_values(x_init)
    init_conds[3] = init_conds[0]*Omega_b

    return init_conds

def derivs(t, init_conds):
    """
    See the solution set for winter term 2015 Galaxy Dynamics problem set 3
    for derivation of this term.
    """

    r_sqr = init_conds[0]**2 + init_conds[2]**2
    r_b = init_conds[0]**2 - init_conds[2]**2
    term1 = R_c**2 + r_sqr


    x_acc = -(V0_sqr*init_conds[0]*(ecce**2 * (r_b/r_sqr) + 1)/term1) + \
            (0.5*V0_sqr)*(((2*(ecce**2)*init_conds[0])/r_sqr) - \
            (2*(ecce**2)*init_conds[0]*r_b)/r_sqr**2)*np.log(term1)

    y_acc = -(V0_sqr*init_conds[2]*(ecce**2 * (r_b/r_sqr) + 1)/term1) + \
            (0.5*V0_sqr)*(((2*(ecce**2)*init_conds[2])/r_sqr) - \
            (2*(ecce**2)*init_conds[2]*r_b)/r_sqr**2)*np.log(term1)


    return np.array([init_conds[1], x_acc, init_conds[3], y_acc])

def integrate(x_initial, ax, label):
    init_conds = compute_vy(x_initial)
    out, tout = bsint.bsintegrate(derivs, init_conds, t0, t1, tacc=1e-14, mxstep=20000)
    ax.plot(out[:,0], out[:,2], color='#424242', alpha=0.5)
    ax.annotate(label, xy=(1,0), xycoords='axes fraction', xytext=(0.95, 0.95),
                ha='right', va='top', fontsize=12)
    ax.set_xlabel('x-position', fontsize=16)
    ax.set_ylabel('y-position', fontsize=16)

if __name__ == '__main__':
    fig = plt.figure(figsize=(22,8))
    x_initials = [cortation(), inner_linblad(), outer_linblad()]
    labels = ['Co-rotation', 'Inner Linblad', 'Outer Linblad']
    for i, (x_init, label) in enumerate(zip(x_initials, labels)):
        ax = plt.subplot(1,3,i)
        integrate(x_init, ax, label)
    plt.tight_layout()
    plt.show()

