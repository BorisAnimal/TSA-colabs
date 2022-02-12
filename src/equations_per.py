import numpy as np
from scipy import integrate
from scipy.signal import find_peaks

from equations_lin import Lin


def period_from_data(theta, t):
    """
    Checks when sign of time series flips.
    :returns mean period - single float
    """
    return period_range_from_data(theta, t)[1]


def period_range_from_data(theta, t, distance=None, debug=False):
    """
    Checks when sign of time series flips.
    :returns tuple of periods (min,mean,max) - tuple of floats
    """
    y = -abs(theta)
    hmin = np.percentile(y, 95)
    hmax = np.max(y)
    peak_ids = find_peaks(y, distance=distance, height=(hmin, hmax))[0]
    cross_times = np.take(t, peak_ids)
    cross_y = np.take(y, peak_ids)
    if debug:
        import matplotlib.pyplot as plt
        plt.plot(t, y)
        plt.scatter(cross_times, cross_y)
        plt.show()
    periods = 2 * np.diff(cross_times)
    # print("periods:", periods)
    if len(periods) > 0:
        mean = periods.mean()
        MIN = periods.min()
        MAX = periods.max()
        return (MIN, mean, MAX)
    else:
        print("Full oscillation is not recognised for data1,25")
        return (0.0, 0.0, 0.0)


def __period_range_from_data(theta, t):
    """
    Checks when sign of time series flips.
    :returns tuple of periods (min,mean,max) - tuple of floats
    """
    cross_times = []
    Is = []
    for i, v in enumerate(theta[1:] * theta[:-1] < 0):
        if v:
            cross_times.append(t[i])
            Is.append(i)
    import matplotlib.pyplot as plt
    plt.plot(t, theta)
    plt.scatter(np.take(t, Is), np.take(theta, Is))
    plt.show()
    periods = 2 * np.diff(cross_times)
    # print("periods:", 2 * half_periods)
    if len(periods) > 0:
        mean = periods.mean()
        MIN = np.take(periods, periods < mean).mean()
        MAX = np.take(periods, periods > mean).mean()
        return (MIN, mean, MAX)
    else:
        print("Full oscillation is not recognised for data1,25")
        return (0.0, 0.0, 0.0)


def period_from_integral(theta0, lin: Lin):
    def f(theta):
        D = lin.D([theta, 0])
        E0 = lin.P([theta0, 0])
        P = lin.P([theta, 0])
        return np.sqrt(D / (E0 - P))

    quad = integrate.quad(f, 0, theta0)[0]
    return 2 ** 1.5 * quad


def period_small_lin(lin):
    r0, L0, R, I1, I2, m, h, g = lin.params
    return 2 * np.pi / r0 * np.sqrt(I1 * L0 / (m * g))
