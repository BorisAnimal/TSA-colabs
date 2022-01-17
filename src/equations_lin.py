from numpy import *
import numpy as np


class Lin:
    def __init__(self, r0, L0, R, I1, I2, m, h, g=9.81, b_theta=0, tau_c=0, b_x=0, F_c=0, alpha0=0.0, mu=0.0):
        self.params = [r0, L0, R, I1, I2, m, h, g]
        self.fric_params = [b_theta, tau_c, b_x, F_c]
        self.mu = mu
        self.alpha0 = alpha0

    def x(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        return L0 - sqrt(L0 ** 2 - (r0 * theta) ** 2)

    def J(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        return theta * r0 ** 2 / sqrt(L0 ** 2 - (theta * r0) ** 2)

    def dJ(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        tmp = L0 ** 2 - (r0 * theta) ** 2
        # if tmp == 0:
        #     return 0
        return dtheta * r0 ** 2 * (1 / np.sqrt(tmp) + (r0 * theta) ** 2 / tmp ** 1.5)

    def T(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        b_theta, tau_c, b_x, F_c = self.fric_params
        dJ = self.dJ(state)
        dx = self.J(state) * dtheta
        return m * dJ * dtheta + m * g + b_x * dx + F_c * tanh(dx * 100)

    def ddtheta(self, state, u=0):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        b_theta, tau_c, b_x, F_c = self.fric_params
        J = self.J(state)
        T = self.T(state)
        return 1 / (I1 + m * J ** 2) * (u - J * T - b_theta * dtheta - tau_c * tanh(dtheta * 100))

    def K(self, state):
        return self.K1(state) + self.K2(state)

    def K1(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        return I1 * dtheta ** 2 / 2

    def K2(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        J = self.J(state)
        return m * J ** 2 * dtheta ** 2 / 2

    def D(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        J = self.J(state)
        return I1 + m * J ** 2

    def P(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        x = self.x(state)
        return m * g * x

    def E(self, state):
        return self.K(state) + self.P(state)
