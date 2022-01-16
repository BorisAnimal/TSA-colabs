import numpy as np
from scipy import integrate

from equations_lin import Lin


def period_from_data(theta, t):
    """
    Checks when sign of time series flips.
    """
    z = []
    for i, v in enumerate(theta[1:] * theta[:-1] < 0):
        if v:
            z.append(t[i])
    a = []
    for i in range(len(z) - 1):
        a.append(z[i + 1] - z[i])
    if len(a) > 0:
        T_aprox = 2 * np.mean(a)
    else:
        print("Full oscillation is not recognised for data")
        T_aprox = 0.0
    return T_aprox


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
