import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_from_data, period_range_from_data, period_from_integral
from params import lin_params_no_fric
import os
import pickle as pkl

from signal_processing import preprocess

plt.rcParams['figure.figsize'] = [5, 4]


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


if __name__ == '__main__':
    # dirpath = "../data_lin_2,5/"
    # lin_params_no_fric['m'] = 2.65
    dirpath = "../data_lin_1,25/"
    lin_params_no_fric['m'] = 2.150
    lin_params_no_fric['I1'] = 22.48e-7
    theta0s = []
    pers = []
    # Scatter dots
    for file in os.listdir(dirpath):
        filepath = os.path.join(dirpath, file)

        with open(filepath, "rb") as f:
        theta = np.array(data['theta'])
        x = np.array(data['x'])
        dtheta = np.array(data['dtheta'])
        E = np.array(data['E'])
        t = np.array(data['t'])
        T = np.array(data['T'])
        T0 = np.percentile(T, 3)
        print("T0: ", T0)

        t_max = max(t)
        t_start = 0.25*t_max  # sec
        Ftheta = preprocess(theta, t, show_plots=False)
        Fdtheta = preprocess(dtheta, t, show_plots=False)
        Fx = preprocess(x, t, show_plots=False)
        FE = preprocess(E, t, show_plots=False)
        t = np.linspace(t_start, 0.9*t_max, len(t))
        theta = Ftheta(t)
        dtheta = Fdtheta(t)
        x = Fx(t)

        per = period_range_from_data(theta, t, distance=100, debug=True)
        print("Period from data:", per)
        theta0 = (abs(np.percentile(theta, 1)) + abs(np.percentile(theta, 99))) / 2
        theta0s.append(theta0)
        pers.append(per)

    for theta0, (MIN, MEAN, MAX) in zip(theta0s, pers):
        plt.vlines(theta0, MIN, MAX, 'b')  # , label='experimental data'
        plt.scatter(theta0, MEAN, c='b')  # , label='experimental data'

    t = np.linspace(0, 20, 10000)
    rot = Lin(**lin_params_no_fric)
    # T1 = T_integral()
    pers = []
    pers_non = []
    thetas0 = np.linspace(min(theta0s) * 0.9, max(theta0s) * 1.3, 100)
    for theta0 in thetas0:
        init_state = [theta0, 0]
        theta = odeint(sys_ode, init_state, t, )[:, 0]
        pers.append(period_from_data(theta, t))
        pers_non.append(period_from_integral(theta0, rot))

    plt.title(r'Period ($T$) of initial state ($\theta_0$)')
    plt.plot(thetas0, pers, linewidth=3, label=r'$T$ from simulations')
    plt.plot(thetas0, pers_non, linewidth=3, linestyle='-.', label=r'$T$ with stiffness energy')
    # plt.hlines(T1, min(thetas0), max(thetas0), color='red', linestyles='--', linewidth=2.0,
    #            label=r'$T$ rotear formula')
    plt.legend()
    plt.xlabel(r'$\theta_0$, rad')
    plt.ylabel(r'$T$, sec')
    # plt.ylim(0, 12)
    plt.show()
