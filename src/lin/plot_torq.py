import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_small_lin, period_from_data, period_from_integral
from params import lin_params_no_fric, lin_params

plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t, T, amp=5e-3):
    theta, dtheta = state
    if T == 0:
        u = 0
    else:
        u = amp * np.sin(2 * np.pi * t / T)
    ddtheta = lin.ddtheta(state, u)
    return dtheta, ddtheta


if __name__ == '__main__':
    t = np.linspace(0, 10, 10000)
    lin = Lin(**lin_params_no_fric)

    theta0 = 100  # initial TSA angle
    T1 = period_small_lin(lin)  # small self oscillations per
    tmp = odeint(sys_ode, [theta0, 0], t, args=(0,))[:, 0]
    T2 = period_from_data(tmp, t)  # self oscillations per by frictionless simulation
    T3 = period_from_integral(theta0, lin)  # self oscillations per by integral

    lin = Lin(**lin_params)
    # for T in [T1, T2, T3]:
    for T in [T3]:
        sol = odeint(sys_ode, [0, 0], t, args=(T,))
        thetas, dthetas = sol.T

        Es = list(map(lin.E, sol))
        K1s = list(map(lin.K1, sol))
        K2s = list(map(lin.K2, sol))
        Ps = list(map(lin.P, sol))
        E_ideal = lin.E([theta0, 0])

        plt.title(f'Period of sin control = {T} sec')
        plt.plot(t, thetas, linewidth=2, label=r'$\theta$')
        plt.plot(t, dthetas, linewidth=2, label=r'$\dot{\theta}$')
        plt.legend()
        plt.grid()
        plt.show()

        plt.title("Phase plot")
        plt.plot(thetas, dthetas)
        plt.grid()
        plt.show()

        plt.title("Energy")
        plt.hlines(E_ideal, min(t), max(t), label=r'Expected energy level')
        plt.plot(t, Es, label=r'Full energy')
        plt.plot(t, K1s, label=r'Motor motion energy')
        plt.plot(t, K2s, label=r'Load motion energy')
        plt.plot(t, Ps, label=r'Load potential energy')
        plt.legend()
        plt.show()
