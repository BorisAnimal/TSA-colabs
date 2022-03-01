import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_rot import Rot
from equations_per import period_from_data, period_from_integral
from params import rot_params_no_fric, rot_params
from utils import harmonic_trajectory_builder

plt.rcParams['figure.figsize'] = [5, 2]


def f(theta0, I):
    params = rot_params_no_fric.copy()
    params['I1'] = I
    rot = Rot(**params)
    return period_from_integral(theta0, rot)  # self oscillations per by integral


def params_plot(major_param: tuple, minor_param: tuple, y_axis_title, title="",
                rote_styles=("-", ":", "-.", "--", ".")):
    label1, x1 = major_param
    label2, x2 = minor_param
    rote_styles = iter(rote_styles)

    for i in x2:
        y = []
        for j in x1:
            y.append(f(**{label1: j, label2: i}))
        plt.plot(x1, y, next(rote_styles), linewidth=3, label=f"{label2} = {i}")
    plt.title(title)
    plt.xlabel(label1)
    plt.ylabel(y_axis_title)
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    theta0s = np.linspace(50, 300, 30)
    Is = np.linspace(1e-7, 1e-4, 60)

    params_plot(
        ("theta_amp", theta0s),
        ("I", [1e-7, 1e-6, 1e-5, 1e-4, ]),
        r"$T$ sec",
        "Self oscillations period (Rot) dependence on params",
    )

    params_plot(
        ("I", Is),
        ("theta_amp", [50, 100, 200, 300]),
        r"$T$ sec",
        "Self oscillations period (Rot) dependence on params",
    )
