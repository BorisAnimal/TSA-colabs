import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_small_lin, period_from_data, period_from_integral
from params import lin_params_no_fric, rot_params_no_fric

plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = lin.ddtheta(state, u)
    return dtheta, ddtheta


def experiment(lin):
    t = np.linspace(0, 80, 10000)
    T1 = period_small_lin(lin)
    pers = []
    pers_non = []
    thetas0 = np.linspace(1, 215, 100)
    for theta0 in thetas0:
        init_state = [theta0, 0]
        theta = odeint(sys_ode, init_state, t, )[:, 0]
        pers.append(period_from_data(theta, t))
        pers_non.append(period_from_integral(theta0, lin))

    plt.title(r'Period ($T$) of initial state ($\theta_0$)')
    plt.plot(thetas0, pers, linewidth=2, label=r'$T$ from simulations')
    plt.plot(thetas0, pers_non, linewidth=2, linestyle='-.', label=r'$T$ nonlinear numerical')
    plt.hlines(T1, min(thetas0), max(thetas0), color='red', linestyles='--', linewidth=2.0,
               label=r'$T$ linear formula')
    plt.legend()
    plt.xlabel(r'$\theta_0$, rad')
    plt.ylabel(r'$T$, sec')
    plt.grid()
    plt.show()


if __name__ == '__main__':
    for i in np.linspace(1.0, 100, 10):
        lin = Lin(**lin_params_no_fric)
        # lin = Lin(**rot_params_no_fric)
        lin.params[3] *= i  # I1
        # lin.params[1] *= i  # I1
        experiment(lin)
