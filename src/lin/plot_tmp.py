import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_small_lin, period_from_data, period_from_integral
from params import lin_params_no_fric

# plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = lin.ddtheta(state, u)
    return dtheta, ddtheta


def experiment(lin):
    t = np.linspace(0, 5, 10000)
    T1 = period_small_lin(lin)
    pers = []
    pers_non = []
    thetas0 = np.linspace(1, 215, 100)
    for theta0 in thetas0:
        init_state = [theta0, 0]
        sol = odeint(sys_ode, init_state, t, )
        theta, dtheta = sol[:, 0], sol[:, 1]
        pers.append(period_from_data(theta, t))
        pers_non.append(period_from_integral(theta0, lin))

        Es = list(map(lin.E, sol))
        K1s = list(map(lin.K1, sol))
        K2s = list(map(lin.K2, sol))
        Ps = list(map(lin.P, sol))

        # plt.title("Motion lin")
        # plt.plot(t, theta, label=r'$\theta$')
        # plt.plot(t, dtheta, label=r'$\dot{\theta}$')
        # plt.legend()
        # plt.grid()
        # plt.show()
        #
        # plt.title("Energy")
        # plt.plot(t, Es, label=r'Full energy')
        # plt.plot(t, K1s, label=r'Motor motion energy')
        # plt.plot(t, K2s, label=r'Load motion energy')
        # plt.plot(t, Ps, label=r'Load potential energy')
        # plt.legend()
        # plt.grid()
        # plt.show()

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
    for i in np.linspace(2.0, 2.0, 1):
        lin = Lin(**lin_params_no_fric)
        lin.params[3] *= i  # I1
        experiment(lin)
