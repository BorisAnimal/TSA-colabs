import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_small_lin, period_from_data, period_from_integral
from params import lin_params_no_fric, lin_params
from utils import harmonic_trajectory_builder

plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t, lin, control):
    theta, dtheta = state
    u = control(state, t)
    ddtheta = lin.ddtheta(state, u)
    return dtheta, ddtheta


def passive_control(state, t):
    return 0


def pd_tracker_builder(trajectory, u_max=3e-3, kp=0.8, kd=0.1):
    def control(state, t):
        theta, dtheta = state
        theta_d, dtheta_d, ddtheta_d = trajectory(t)
        u = kp * (theta_d - theta) + kd * (dtheta_d - dtheta)
        if abs(u) > u_max:
            u = u_max * np.sign(u)
        return u

    return control


if __name__ == '__main__':
    t = np.linspace(0, 10, 10000)
    lin = Lin(**lin_params_no_fric)

    theta0 = 100  # initial TSA angle
    T1 = period_small_lin(lin)  # small self oscillations per
    tmp = odeint(sys_ode, [theta0, 0], t, args=(lin, passive_control))[:, 0]
    T2 = period_from_data(tmp, t)  # self oscillations per by frictionless simulation
    T3 = period_from_integral(theta0, lin)  # self oscillations per by integral

    lin = Lin(**lin_params)
    # for T in [T1, T2, T3]:
    for T in [T3]:
        # T = T/2
        # T = T*2**0.5
        pd_control = pd_tracker_builder(
            harmonic_trajectory_builder(A0=0, A=100, freq=1 / T),
            u_max=6e-3,
            kp=0.8,
            kd=0.2,
        )
        sol = odeint(sys_ode, [0, 0], t, args=(lin, pd_control))
        thetas, dthetas = sol.T

        Es = list(map(lin.E, sol))
        K1s = list(map(lin.K1, sol))
        K2s = list(map(lin.K2, sol))
        Ps = list(map(lin.P, sol))
        E_ideal = lin.E([theta0, 0])
        us = [pd_control(state, tt) for state, tt in zip(sol, t)]

        thetas = np.array(thetas)
        dthetas = np.array(dthetas)
        us = np.array(us)

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
        plt.grid()
        plt.show()

        plt.title("Torque of time")
        plt.plot(t, us)
        plt.grid()
        plt.show()

        plt.title("Torque of phase")
        plt.plot(thetas, us)
        # plt.plot(abs(thetas), abs(us), )
        plt.grid()
        plt.show()
