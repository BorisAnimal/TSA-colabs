import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks

from equations_lin import Lin
from equations_rot import Rot
from equations_per import period_from_data, period_range_from_data, period_from_integral
from params import rot_params_no_fric
import os
import pickle as pkl

from signal_processing import preprocess

plt.rcParams['figure.figsize'] = [5, 4]
plt.rcParams['axes.facecolor'] = 'white'


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


if __name__ == '__main__':
    filepath = "../data/111data_ebc.pkl"
    rot_params_no_fric['m'] = 2.6
    with open(filepath, "rb") as f:
        data = pkl.load(f)
    theta = np.array(data['theta'])
    alpha = np.array(data['alpha'])
    dtheta = np.array(data['dtheta'])
    E = np.array(data['E'])
    t = np.array(data['t'])
    T = np.array(data['T'])
    alpha0 = alpha[1]

    t_start = 5  # sec
    t_end = max(t)
    Ftheta = preprocess(theta, t, show_plots=False)
    Fdtheta = preprocess(dtheta, t, show_plots=False)
    Falpha = preprocess(alpha, t, show_plots=False)
    FE = preprocess(E, t, show_plots=False)
    t = np.linspace(t_start, t_end, 10000)
    theta = Ftheta(t)
    dtheta = Fdtheta(t)
    alpha = Falpha(t)

    ## Step 2
    t_start = 5  # sec
    tau = np.linspace(t_start, min(90, max(t)), 1000)
    alpha = Falpha(tau)
    peaks, _ = find_peaks(-abs(alpha), height=np.percentile(-abs(alpha), 60), distance=20)
    tpeaks = np.take(tau, peaks)

    pers = (tpeaks[1:] - tpeaks[:-1]) * 2
    pers = [i for i in pers if i > 0.3]
    print(pers)
    theta0s = []
    for i in range(len(peaks) - 1):
        thetas = theta
    # bar(tpeaks[:len(pers)], pers)
    # grid()
    # title("Full period for data")
    # xlabel("time")
    # ylabel("period")
    # show()
    #
    # scatter(tpeaks, Falpha(tpeaks), c='r')
    # plot(tau, alpha)
    # grid()
    # title("Alpha of process")
    # xlabel("time")
    # ylabel("alpha, rad")
    # show()

    # t = np.linspace(0, 20, 10000)
    # rot = Rot(T0=10, alpha0=alpha0, U_coef=0.0, **rot_params_no_fric)
    # # rot = Rot(T0=0.0, **rot_params_no_fric)
    # # T1 = T_integral()
    # pers = []
    # pers_non = []
    # thetas0 = np.linspace(min(theta_amps) * 0.1, min(max(theta_amps) * 1.3, 270), 100)
    # for theta_amp in thetas0:
    #     init_state = [theta_amp, 0]
    #     theta, dtheta = odeint(sys_ode, init_state, t, ).T
    #     pers.append(period_from_data(theta, t))
    #     pers_non.append(period_from_integral(theta_amp, rot))
    #
    # plt.title(r'Period ($T$) of initial state ($\theta_0$)')
    # plt.plot(thetas0, pers, linewidth=3, label=r'$T$ from simulations')
    # plt.plot(thetas0, pers_non, linewidth=3, linestyle='-.', label=r'$T$ integral without stiffness energy')
    # # plt.hlines(T1, min(thetas0), max(thetas0), color='red', linestyles='--', linewidth=2.0,
    # #            label=r'$T$ rotear formula')
    # plt.legend()
    # plt.xlabel(r'$\theta_0$, rad')
    # plt.ylabel(r'$T$, sec')
    # plt.grid()
    # # plt.ylim(0, 12)
    # plt.show()