import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from equations_lin import Lin
from equations_per import period_from_data, period_from_integral
from params import lin_params_no_fric


def sysode(state, t):
    theta, dtheta = state
    u = 0
    ddtheta = lin.ddtheta(state, u)
    return dtheta, ddtheta


# Lin
if __name__ == '__main__':
    t = np.linspace(0, 3, 10000)
    lin = Lin(**lin_params_no_fric)
    theta0 = 150
    init_state = [theta0, 0]

    sol = odeint(sysode, init_state, t)
    theta, dtheta = sol[:, 0], sol[:, 1]

    print("Period from simulation:", period_from_data(theta, t))
    print("Period from integral:", period_from_integral(theta0, lin))

    Es = list(map(lin.E, sol))
    K1s = list(map(lin.K1, sol))
    K2s = list(map(lin.K2, sol))
    Ps = list(map(lin.P, sol))

    plt.title("Motion lin")
    plt.plot(t, theta, label=r'$\theta$')
    plt.plot(t, dtheta, label=r'$\dot{\theta}$')
    plt.legend()
    plt.grid()
    plt.show()

    plt.title("Energy")
    plt.plot(t, Es, label=r'Full energy')
    plt.plot(t, K1s, label=r'Motor motion energy')
    plt.plot(t, K2s, label=r'Load motion energy')
    plt.plot(t, Ps, label=r'Load potential energy')
    plt.legend()
    plt.grid()
    plt.show()


