import pickle as pkl

import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from params import rot_params_no_fric
from signal_processing import preprocess

plt.rcParams['figure.figsize'] = [5, 4]
plt.rcParams['axes.facecolor'] = 'white'

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
    tau = np.linspace(t_start, min(50, max(t)), 1000)
    alpha = Falpha(tau)

    # select only positive/negative part of rotations
    theta = Ftheta(tau)
    ids = np.array(theta) <= 0
    alpha = np.array(alpha)[ids]
    tau = tau[ids]
    # remove time steps in skipped ranges
    for i, t in enumerate(tau[1:]):
        dt = t - tau[i]
        if dt > 0.1:
            tau[i + 1:] = tau[i + 1:] - dt

    peaks, _ = find_peaks(-abs(alpha), height=np.percentile(-abs(alpha), 60), distance=10)
    tpeaks = np.take(tau, peaks)
    apeaks = np.array([max(alpha[i:j]) for (i, j) in zip(peaks[:-1], peaks[1:])])

    pers = (tpeaks[1:] - tpeaks[:-1]) * 2
    pers = [i for i in pers if i > 0.3]
    print("Periods:", pers)
    theta0s = []


    def func(x, a, b, c):
        return a + b * (1 - np.cos(c * x))


    popt, pcov = curve_fit(func, apeaks, pers)

    grid()
    # scatter(tpeaks[:len(pers)], pers)
    # plot([apeaks[0], apeaks[-1]], [pers[0], pers[-1]], c='gray')
    plot(alpha, func(alpha, *popt), c='gray', label="Fitted sinusoid")
    scatter(apeaks, pers, c='r', label='Period sample')
    xlim(left=apeaks[-1] * 0.9)
    title("Full period for data")
    xlabel(r"$\alpha$, rad")
    ylabel(r"$T$, sec")
    legend()
    show()

    grid()
    scatter(tpeaks, alpha[peaks], c='r', label='End of wave')
    plot(tau, alpha, label='Joint\'s angle')
    title("Oscillator angle under descending EPC")
    xlabel(r"$t$, sec")
    ylabel(r"$\alpha$, rad")
    legend()
    show()
