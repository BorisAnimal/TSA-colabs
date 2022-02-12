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


def params_plot(major_param: tuple, minor_param: tuple, y_axis_title, title="",
                line_styles=("-", ":", "-.", "--", ".")):
    label1, x1 = major_param
    label2, x2 = minor_param
    line_styles = iter(line_styles)

    for i in x2:
        y = []
        for j in x1:
            y.append(f(**{label1: j, label2: i}))
        plt.plot(x1, y, next(line_styles), label=f"{label2} = {i}")
    plt.title(title)
    plt.xlabel(label1)
    plt.ylabel(y_axis_title)
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    theta0s = np.linspace(100, 200, 60)
    Is = np.linspace(1e-7, 1e-4, 60)

    params_plot(
        ("theta0", theta0s),
        ("I", [1e-7, 1e-6,3.848e-6, 1e-5 ]),
        r"$T$ sec",
        "Self oscillations period dependence on params",
    )

    params_plot(
        ("I", Is),
        ("theta0", [50, 100, 150, 200]),
        r"$T$ sec",
        "Self oscillations period dependence on params",
    )
