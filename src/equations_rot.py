from numpy import *

from equations_lin import Lin


class Rot(Lin):

    def alpha(self, state):
        r0, L0, R, I1, I2, m, h, g = self.params
        alpha = self.x(state) / R
        return alpha + self.alpha0

    def dalpha(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        dx = self.J(state) * dtheta
        dalpha = dx / R
        return dalpha

    def T(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        b_theta, tau_c, b_x, F_c = self.fric_params
        dJ = self.dJ(state)
        dx = self.J(state) * dtheta
        dalpha = self.dalpha(state)
        alpha = self.alpha(state)
        return I2 * dJ / R * dtheta + m * g * h / R * sin(alpha) + b_x * dalpha + F_c * tanh(dalpha * 100)

    def ddtheta(self, state, u=0):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        b_theta, tau_c, b_x, F_c = self.fric_params
        mu = self.mu
        J = self.J(state)
        T = self.T(state)
        return 1 / (I1 + (J + mu) * I2 * J / R) * (u - (mu + J) * T - b_theta * dtheta - tau_c * tanh(dtheta * 100))

    def K(self, state):
        return self.K1(state) + self.K2(state)

    def K2(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        J = self.J(state)
        return I2 * (J * h / R) ** 2 / 2 * dtheta ** 2

    def D(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        J = self.J(state)
        return I1 + I2 * (J * h / R) ** 2

    def P(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        alpha = self.alpha(state)
        return m * g * h * (1 - cos(alpha)) - m * g * h * (1 - cos(self.alpha0))
