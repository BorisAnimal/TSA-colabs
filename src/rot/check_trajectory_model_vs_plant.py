import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from equations_lin import Lin
from equations_rot import Rot
from equations_per import period_from_data, period_range_from_data, period_from_integral
from params import rot_params_no_fric
import os
import pickle as pkl

from signal_processing import preprocess

plt.rcParams['figure.figsize'] = [5.67, 4]


def sys_ode(state, t):
    theta, dtheta = state
    u = 0.0
    ddtheta = rot.ddtheta(state, u)
    return dtheta, ddtheta


def theta_inverse(E, theta_init, plant: Lin):
    theta_init *= 0.5
    P = plant.P((theta_init, 0.0))
    assert P < E, str((P, E))
    c = 1.0
    while P < E:
        c += 0.001
        theta = c * theta_init
        P = plant.P((theta, 0.0))
    print(P, E, c)
    return theta


if __name__ == '__main__':
    dirpath = "../data2,5/"
    rot_params_no_fric['m'] = 2.65
    # dirpath = "../data1,25/"
    # rot_params_no_fric['m'] = 1.3
    theta0s = []
    pers = []
    # Scatter dots
    for file in os.listdir(dirpath):
        filepath = os.path.join(dirpath, file)
        with open(filepath, "rb") as f:
            data = pkl.load(f)
        theta = np.array(data['theta'])
        alpha = np.array(data['alpha'])
        dtheta = np.array(data['dtheta'])
        E = np.array(data['E'])
        t = np.array(data['t'])
        T = np.array(data['T'])

        t_start = 5  # sec
        Ftheta = preprocess(theta, t, show_plots=False)
        Fdtheta = preprocess(dtheta, t, show_plots=False)
        Falpha = preprocess(alpha, t, show_plots=False)
        FE = preprocess(E, t, show_plots=False)
        FT = preprocess(T, t, show_plots=False)
        t = np.linspace(t_start, max(t), 10000)
        theta = Ftheta(t)
        dtheta = Fdtheta(t)
        alpha = Falpha(t)
        T = FT(t)

        per = period_range_from_data(theta, t)
        print("Period from data:", per)

        plant = Rot(alpha0=alpha[0], **rot_params_no_fric)
        E0 = plant.E((Ftheta(t_start), Fdtheta(t_start)))
        print("E0s comparison:", E0, FE(t_start))
        theta0 = theta_inverse(E0, max(theta), plant)

        print("Period from params by integral:", period_from_integral(theta0, plant))
        theta0s.append(theta0)
        pers.append(per)

        # t = np.linspace(0, 20, len(theta))
        rot = Rot(alpha0=min(alpha), **rot_params_no_fric)
        init_state = [theta[0], dtheta[0]]
        theta_sim, dtheta_sim = odeint(sys_ode, init_state, t, ).T

        plt.plot(t, theta, linewidth=4, label='experiment')
        plt.plot(t, theta_sim, linewidth=3, label='simulation')
        plt.grid()
        plt.legend()
        plt.title('theta')
        # plt.ylabel('dtheta')
        # plt.xlabel('theta')
        plt.show()

    # plt.plot(t,dtheta, linewidth=4, label='experiment')
    # plt.plot(t,dtheta_sim, linewidth=3, label='simulation')
    # plt.grid()
    # plt.legend()
    # plt.title('dtheta')
    # # plt.ylabel('dtheta')
    # # plt.xlabel('theta')
    # plt.show()
    #
    # plt.plot(t,alpha, linewidth=4, label='experiment')
    # plt.plot(t,rot.alpha([theta_sim, dtheta_sim]), linewidth=3, label='simulation')
    # plt.grid()
    # plt.legend()
    # plt.title(file + ' alpha')
    # # plt.ylabel('dtheta')
    # # plt.xlabel('theta')
    # plt.show()
    #
    # plt.plot(t, T, linewidth=4, label='experiment')
    # plt.plot(t, rot.T([theta_sim, dtheta_sim]), linewidth=3, label='simulation')
    # plt.grid()
    # plt.legend()
    # plt.title('Tension force T')
    # # plt.ylabel('dtheta')
    # # plt.xlabel('theta')
    # plt.show()
