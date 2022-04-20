import pickle as pkl

import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from scipy.signal import find_peaks

from params import rot_params_no_fric
from signal_processing import preprocess

plt.rcParams['figure.figsize'] = [5, 4]
plt.rcParams['axes.facecolor'] = 'white'

if __name__ == '__main__':
    filepath = "../../data/data_ebc(15sec).pkl"
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
    t_end = min(max(t), 22.5)
    Ftheta = preprocess(theta, t, show_plots=False)
    Fdtheta = preprocess(dtheta, t, show_plots=False)
    Falpha = preprocess(alpha, t, show_plots=False)
    FE = preprocess(E, t, show_plots=False)
    t = np.linspace(t_start, t_end, 10000)
    theta = Ftheta(t)
    dtheta = Fdtheta(t)
    alpha = Falpha(t)

    ## Step 2
    tau = np.linspace(t_start, t_end, 1000)
    theta = Ftheta(tau)

    peaks, _ = find_peaks(-abs(theta), height=np.percentile(-abs(theta), 60), distance=30)
    peaks = peaks[:-2]
    tpeaks = np.take(tau, peaks)
    ypeaks = np.take(theta, peaks)

    hpers = (tpeaks[1:] - tpeaks[:-1])
    print("HALF Periods:", hpers)
    theta0s = []

    # Main plot
    grid()
    xticks(np.linspace(5, 22.5, 8))
    yticks(np.linspace(-200, 200, 5))
    xlim(5,22.5)
    ylim(-250, 250)
    # vlines(15, -200, 200)
    vlines(tpeaks, -190, 190, alpha=0.3, linewidth=2)
    plot(tau, theta, linewidth=3)
    scatter(tpeaks, ypeaks, c='r', zorder=10)
    ylabel(r"$\theta$, rad")
    xlabel(r"$t$, sec")
    show()
