import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from tqdm import tqdm

from equations_rot import Rot
from equations_per import period_from_data, period_from_integral
from params import rot_params_no_fric, rot_params
from utils import harmonic_trajectory_builder

plt.rcParams['figure.figsize'] = [12, 4]


def sys_ode(state, t, rot, control):
    theta, dtheta = state
    u = control(state, t)
    ddtheta = rot.ddtheta(state, u)
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
    t = np.linspace(0, 100, 10000)
    dt = t[1] - t[0]
    rot = Rot(**rot_params_no_fric)

    theta0 = 100  # initial TSA angle
    # T1 = T_small_rot(rot)  # small self oscillations per
    tmp = odeint(sys_ode, [theta0, 0], t, args=(rot, passive_control))[:, 0]
    T2 = period_from_data(tmp, t)  # self oscillations per by frictionless simulation
    T3 = period_from_integral(theta0, rot)  # self oscillations per by integral
    TT = T3

    rot = Rot(**rot_params)
    spend_energy = []
    max_amp = []
    Ts = TT * np.linspace(0.9, 1.9, 50)
    for T in tqdm(Ts):
        pd_control = pd_tracker_builder(
            harmonic_trajectory_builder(A0=0, A=100, freq=1 / T),
            u_max=5e-3,
            kp=4,
            kd=1,
        )
        sol = odeint(sys_ode, [0, 0], t, args=(rot, pd_control))
        us = [pd_control(state, tt) for state, tt in zip(sol, t)]
        us = np.array(us)
        max_amp.append(max(abs(sol[:, 0])))
        spend_energy.append(sum(abs(us * dt)))
    T_opt = Ts[np.argmin(spend_energy)]
    print(f"T_self {TT} | T_opt {T_opt} (ratio {T_opt / TT})")

    plt.title(f'Torque integral spend for oscillations')
    plt.plot(Ts, spend_energy, linewidth=2)
    plt.vlines(TT, min(spend_energy), max(spend_energy), "r", linewidth=2, label=r"$T_{self}$")
    plt.vlines(T_opt, min(spend_energy), max(spend_energy), "g", label=r"$T_{opt}$")
    plt.xlabel(r"$T_{traj}$")
    plt.ylabel(r"Cumulative torque, $E$")
    plt.grid()
    plt.legend()
    plt.show()

    T_opt = Ts[np.argmax(max_amp)]
    print(f"T_self {TT} | T_opt {T_opt} (ratio {T_opt / TT})")
    plt.title(f'Max amp during oscillations')
    plt.plot(Ts, max_amp, linewidth=2)
    plt.vlines(TT, min(max_amp), max(max_amp), "r", linewidth=2, label=r"$T_{self}$")
    plt.vlines(T_opt, min(max_amp), max(max_amp), "g", label=r"$T_{opt}$")
    plt.xlabel(r"$T_{traj}$")
    plt.ylabel(r"Max amp, rad")
    plt.grid()
    plt.legend()
    plt.show()
