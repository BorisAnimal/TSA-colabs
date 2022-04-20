import matplotlib.pyplot as plt
import numpy as np

from equations_per import period_from_integral
from equations_rot import Rot
from params import rot_params_no_fric

plt.rcParams['figure.figsize'] = [5, 4.5]


def params_plot(*xs_ys_labels, x_axis_title="", y_axis_title="", title="",
                line_styles=("-", ":", "-.", "--", ".")):
    line_styles = iter(line_styles)
    for x, y, label in xs_ys_labels:
        plt.plot(x, y, next(line_styles), linewidth=3, label=label)
    plt.title(title)
    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    thetas1 = np.linspace(50, 270, 30)
    Is1 = np.linspace(1e-6, 1e-4, 60)

    thetas2 = [100, 175, 270]
    Is2 = [1e-5, 5e-5, 1e-4]

    xs_ys_labels = []
    for theta in thetas2:
        ys = []
        for I in Is1:
            params = rot_params_no_fric.copy()
            params['I1'] = I
            rot = Rot(**params)
            ys.append(period_from_integral(theta, rot))
        xs_ys_labels.append((Is1, ys, r'$\alpha$ = ' + f'{rot.alpha([theta, 0.0]):.2f}'))
    params_plot(
        *xs_ys_labels,
        x_axis_title=r"$I$, $kgm^2$",
        y_axis_title=r"$T$, $sec$",
        title="Self oscillations period (Rot) dependence on params",
    )

    ## 2nd plot
    xs_ys_labels = []
    for I in Is2:
        params = rot_params_no_fric.copy()
        params['I1'] = I
        rot = Rot(**params)
        ys = []
        for theta in thetas1:
            ys.append(period_from_integral(theta, rot))
        xs_ys_labels.append((
            rot.alpha([thetas1, np.zeros_like(thetas1)]),
            ys,
            r'$I$ = ' + f'{I}'
        ))
    params_plot(
        *xs_ys_labels,
        x_axis_title=r"$\alpha_{amp}$, $rad$",
        y_axis_title=r"$T$, $sec$",
        title="Self oscillations period (Rot) dependence on params",
    )
