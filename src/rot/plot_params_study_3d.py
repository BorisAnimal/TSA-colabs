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
    rot = Rot(**rot_params_no_fric)
    theta0s = np.linspace(80, 270, 20)
    Is = np.linspace(1e-5, 1e-3, 20)
    X, Y, Z = [], [], []
    for theta0 in tqdm(theta0s):
        for I in Is:
            X.append(rot.alpha((theta0, 0.0)))
            Y.append(I)
            Z.append(f(theta0, I))

    ax = plt.axes(projection='3d')
    transform = lambda a: np.array(a).reshape((20, 20))
    x = transform(X)
    y = transform(Y)
    z = transform(Z)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(x, y, z, alpha=0.7, cmap='viridis', edgecolor='none')

    ax.set_xlabel(r"$\alpha_{amp}$")
    ax.set_ylabel(r"$Inertia$")
    ax.set_zlabel(r"$T_{integral}$, $sec$")

    # ax.set_title('Surface plot')
    plt.show()
