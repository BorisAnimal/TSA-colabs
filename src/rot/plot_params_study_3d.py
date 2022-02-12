import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.integrate import odeint

from equations_rot import Rot
from equations_per import period_from_data, period_from_integral
from params import rot_params_no_fric, rot_params
from utils import harmonic_trajectory_builder


plt.rcParams['figure.figsize'] = [5, 5]


def f(theta0, I):
    params = rot_params_no_fric.copy()
    params['I1'] = I
    rot = Rot(**params)
    return period_from_integral(theta0, rot)  # self oscillations per by integral


if __name__ == '__main__':
    theta0s = np.linspace(80, 300, 20)
    Is = np.linspace(1e-6, 1e-3, 20)
    X, Y, Z = [], [], []
    for theta0 in tqdm(theta0s):
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
