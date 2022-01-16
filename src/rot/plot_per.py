import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_rot import Rot
from equations_per import period_from_data, period_from_integral
from params import rot_params_no_fric

plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


if __name__ == '__main__':
    t = np.linspace(0, 20, 10000)
    rot = Rot(**rot_params_no_fric)
    # T1 = T_integral()

    pers = []
    pers_non = []
    thetas0 = np.linspace(25, 200, 100)
    for theta0 in thetas0:
        init_state = [theta0, 0]
        theta = odeint(sys_ode, init_state, t, )[:, 0]
        pers.append(period_from_data(theta, t))
        pers_non.append(period_from_integral(theta0, rot))

    plt.title(r'Period ($T$) of initial state ($\theta_0$)')
    plt.plot(thetas0, pers, linewidth=2, label=r'$T$ from simulations')
    plt.plot(thetas0, pers_non, linewidth=2, linestyle='-.', label=r'$T$ nonlinear numerical')
    # plt.hlines(T1, min(thetas0), max(thetas0), color='red', linestyles='--', linewidth=2.0,
    #            label=r'$T$ rotear formula')
    plt.legend()
    plt.xlabel(r'$\theta_0$, rad')
    plt.ylabel(r'$T$, sec')
    plt.show()
