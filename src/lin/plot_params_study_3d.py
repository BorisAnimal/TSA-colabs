import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_small_lin, period_from_data, period_from_integral
from params import lin_params_no_fric, lin_params
from utils import harmonic_trajectory_builder


# plt.rcParams['figure.figsize'] = [12, 4]


def f(theta0, I):
    params = lin_params_no_fric.copy()
    params['I1'] = I
    lin = Lin(**params)
    return period_from_integral(theta0, lin)  # self oscillations per by integral


if __name__ == '__main__':
    theta0s = np.linspace(1, 200, 30)
    Is = np.linspace(1e-7, 1e-3, 60)
    X, Y, Z = [], [], []
    for theta0 in theta0s:
        for I in Is:
            X.append(theta0)
            Y.append(I)
            Z.append(f(theta0, I))

    ax = plt.axes(projection='3d')
    ax.scatter(X, Y, Z, c=Z, cmap='viridis', linewidth=0.5)

    ax.set_xlabel(r"$\theta_0$")
    ax.set_ylabel(r"$Inertia$")
    ax.set_zlabel(r"$T_{integral}$")
    plt.show()
