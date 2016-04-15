import sys
import numpy as np
import matplotlib.pyplot as plt
import bsint

t0 = 0
t1 = 10e5
energy = 1
num_bodies = 1
alpha = 1e-6

def create_initial_values(x_init):
    Y0 = np.zeros((4*num_bodies))
    Y0[0] = x_init # x-position
    Y0[1] = 0.0    # x-velocity
    Y0[2] = 0.0    # y-position
    Y0[3] = 0.0    # y-velocity

    return Y0

def compute_vy(b, Y0):
    '''
    Compute the v-component of the velocity,
    assume x-velocity is initally zero,
    from the expression for totall energy
    '''
    vy = (2*energy)**0.5

    Y0[3] = vy
    return Y0

def derivs(t, Y0):
    x_acc = -(alpha*Y0[0]*6)/(Y0[0]**2 + Y0[2]**2)**(4.)
    y_acc = -(alpha*Y0[2]*6)/(Y0[0]**2 + Y0[2]**2)**(4.)

    return np.array([Y0[1], x_acc, Y0[3], y_acc])

if __name__=='__main__':

    x_inits = np.linspace(1e-3, 1e-1, num=5) # b_0
    for b in x_inits:
        initial_conds = create_initial_values(b)
        initial_conds = compute_vy(b, initial_conds)
        out, tout = bsint.bsintegrate(derivs, initial_conds,
                t0, t1)
        plt.plot(out[:,0], out[:,2], label='b={}'.format(b), lw=2)

    plt.xlabel("x position", fontsize=20)
    plt.ylabel("y position", fontsize=20)
    plt.show()
    #plt.legend(frameon=False, loc='upper left')
    #plt.savefig('scatter_orbits.pdf')


