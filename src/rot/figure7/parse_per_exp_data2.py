import os
import pickle as pkl

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit

from equations_per import period_range_from_data, period_from_integral
from equations_rot import Rot
from params import rot_params_no_fric
from signal_processing import preprocess

plt.rcParams['figure.figsize'] = [5, 4]


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


if __name__ == '__main__':
    dirpath = "../../data2,5/"
    rot_params_no_fric['m'] = 2.6
    # dirpath = "../../data1,25/"
    # rot_params_no_fric['m'] = 1.3
    theta_amps = []
    alpha_amps = []
    pers = []
    # Scatter dots
    for num, file in enumerate(os.listdir(dirpath)):
        print(file)
        filepath = os.path.join(dirpath, file)
        with open(filepath, "rb") as f:
            data = pkl.load(f)
        theta = np.array(data['theta'])
        alpha = np.array(data['alpha'])
        dtheta = np.array(data['dtheta'])
        E = np.array(data['E'])
        t = np.array(data['t'])
        T = np.array(data['T'])
        # T0 = T[1]
        alpha0 = alpha[1]
        # T0 = np.percentile(T, 3)
        # print("T0: ", T0)

        t_start = 5  # sec
        t_end = max(t)
        Ftheta = preprocess(theta, t, show_plots=False)
        Fdtheta = preprocess(dtheta, t, show_plots=False)
        Falpha = preprocess(alpha, t, show_plots=False)
        FE = preprocess(E, t, show_plots=False)
        t = np.linspace(t_start, t_end, 1000)
        theta = Ftheta(t)
        dtheta = Fdtheta(t)
        alpha = Falpha(t)

        per = period_range_from_data(theta, t)
        print("Period from data:", per)
        theta_amp = (abs(np.percentile(theta, 1)) + abs(np.percentile(theta, 99))) / 2
        theta_amps.append(theta_amp)

        rot = Rot(U_coef=0.0, **rot_params_no_fric)
        alpha_amp = rot.alpha([theta_amp, 0.0])
        alpha_amps.append(alpha_amp)
        pers.append(per)
        (MIN, MEAN, MAX) = per
        label1 = "Experiment period range" if num == 0 else None
        plt.vlines(alpha_amp, MIN, MAX, 'b', zorder=9)  # , label='experimental data'
        plt.scatter(alpha_amp, MEAN, c='b', zorder=10, label=label1)  # , label='experimental data'

        # rot = Rot(T0=0, alpha0=min(alpha), U_coef=0.0, **rot_params_no_fric)
        # init_state = [theta[0], dtheta[0]]
        # theta_sim, dtheta_sim = odeint(sys_ode, init_state, t, ).T
        # (MIN, MEAN, MAX) = period_range_from_data(theta_sim, t)
        # print("Model per range:", (MIN, MEAN, MAX))
        # if MIN > 1.0:
        #     label2 = "Simulation period from experiment init state" if num == 0 else None
        #     plt.vlines(theta_amp, MIN, MAX, 'r')  # , label='experimental data'
        #     plt.scatter(theta_amp, MEAN, c='r', label=label2)  # , label='experimental data'

    # Fit model on data
    idx = np.argsort(alpha_amps)
    xs = np.array(alpha_amps)[idx]
    ys = np.array([i[1] for i in pers])[idx]
    def func1(x, a, b, c):
        return a + b * (1 - np.cos(c * x))

    def func2(x, a, b, c):
        return a*x**2 + b*x + c

    func = func2

    popt, _ = curve_fit(func, xs, ys)
    xs2 = np.arange(0.7, 2.2, 0.05)
    plt.plot(xs2, func(xs2, *popt), c='gray', linewidth=3, label="Fitted model")

    t = np.linspace(0, 20, 10000)
    rot = Rot(U_coef=0.0, **rot_params_no_fric)
    # rot = Rot(**rot_params_no_fric)
    # T1 = T_integral()
    pers = []
    pers_non = []
    thetas0 = np.linspace(min(theta_amps) * 0.5, min(max(theta_amps) * 1.1, 270), 100)
    alphas0 = []
    for theta_amp in thetas0:
        init_state = [theta_amp, 0]
        theta, dtheta = odeint(sys_ode, init_state, t, ).T
        alphas0.append(rot.alpha(init_state))
        # pers.append(period_from_data(theta, t))
        pers_non.append(period_from_integral(theta_amp, rot))

    plt.title(r'Period ($T$) of initial state ($\alpha_0$)')
    # plt.plot(thetas0, pers, linewidth=3, label=r'$T$ from simulations')
    # plt.plot(thetas0, pers_non, linewidth=3, linestyle='-.', label=r'$T$ integral without stiffness energy')
    plt.plot(alphas0, pers_non, linewidth=3, c='r', label=r'$T$, sec by integral equation')
    # plt.hlines(T1, min(thetas0), max(thetas0), color='red', linestyles='--', linewidth=2.0,
    #            label=r'$T$ rotear formula')
    plt.legend()
    plt.xlabel(r'$\alpha_0$, rad')
    plt.ylabel(r'$T$, sec')
    plt.grid()
    # plt.ylim(top=6)
    plt.show()
