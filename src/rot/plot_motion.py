import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from equations_rot import Rot
from params import rot_params_no_fric, rot_params


def sysode(state, t):
    theta, dtheta = state
    u = 0
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


# Rot
if __name__ == '__main__':
    t = np.linspace(0, 20, 10000)
    rot = Rot(**rot_params_no_fric, alpha0=0.0)
    # rot = Rot(**rot_params)
    theta0 = 200
    init_state = [theta0, 0]

    sol = odeint(sysode, init_state, t)
    theta, dtheta = sol[:, 0], sol[:, 1]
    ddtheta = rot.ddtheta(sol.T)

    Es = list(map(rot.E, sol))
    K1s = list(map(rot.K1, sol))
    K2s = list(map(rot.K2, sol))
    Ps = list(map(rot.P, sol))
    Alphas = list(map(rot.alpha, sol))
    r0, L0, R, I1, I2, m, h, g = rot.params
    Ts = rot.T(sol.T) + I2 * rot.J(sol.T) / R * rot.ddtheta(sol.T)

    plt.title("TSA tension T(t)")
    plt.plot(t, Ts)
    plt.grid()
    plt.show()

    plt.title("Alphas = x/R")
    plt.plot(t, Alphas)
    plt.grid()
    plt.show()

    plt.title("Motion (rotational)")
    plt.plot(t, theta, label=r'$\theta$')
    plt.plot(t, dtheta, label=r'$\dot{\theta}$')
    plt.plot(t, ddtheta, label=r'$\ddot{\theta}$')
    plt.legend()
    plt.grid()
    plt.show()

    plt.title("Energy (rotational)")
    plt.plot(t, Es, label=r'Full energy')
    plt.plot(t, K1s, label=r'Motor motion energy')
    plt.plot(t, K2s, label=r'Load motion energy')
    plt.plot(t, Ps, label=r'Load potential energy')
    plt.legend()
    plt.grid()
    plt.show()
