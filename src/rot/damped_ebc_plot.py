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
    filepath = "../data/112data_ebc.pkl"
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
    t_end = min(max(t), 55)
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

    bar(tpeaks[:len(pers)], pers)
    grid()
    title("Full period for data")
    xlabel("time")
    ylabel("period")
    show()

    scatter(tpeaks, Falpha(tpeaks), c='r')
    plot(tau, alpha)
    grid()
    title("Alpha of process")
    xlabel("time")
    ylabel("alpha, rad")
    show()
