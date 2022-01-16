import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_rot import Rot
from equations_per import period_from_data, period_from_integral
from params import rot_params_no_fric, rot_params

plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t, T, amp=5e-3):
    theta, dtheta = state
    if T == 0:
        u = 0
    else:
        u = amp * np.sin(2 * np.pi * t / T)
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


if __name__ == '__main__':
    t = np.linspace(0, 100, 10000)
    rot = Rot(**rot_params_no_fric)

    theta0 = 100  # initial TSA angle
    # T1 = T_small_rot(rot)  # small self oscillations per
    tmp = odeint(sys_ode, [theta0, 0], t, args=(0,))[:, 0]
    T2 = period_from_data(tmp, t)  # self oscillations per by frictionless simulation
    T3 = period_from_integral(theta0, rot)  # self oscillations per by integral
    print(T2, T3)
    rot = Rot(**rot_params)
    # for T in [T1, T2, T3]:
    for T in [T3]:
        sol = odeint(sys_ode, [0, 0], t, args=(T,))
        thetas, dthetas = sol.T

        Es = list(map(rot.E, sol))
        K1s = list(map(rot.K1, sol))
        K2s = list(map(rot.K2, sol))
        Ps = list(map(rot.P, sol))
        E_ideal = rot.E([theta0, 0])

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
