from numpy import *

from equations_lin import Lin


class Rot(Lin):

    def __init__(self, T0=0.0, alpha0=0.0, U_coef=0.0, **kwargs):
        super().__init__(**kwargs)
        self.T0 = T0
        self.alpha0 = alpha0
        self.U_coef = U_coef

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
        dalpha = self.dalpha(state)
        alpha = self.alpha(state)
        return self.T0 + I2 * dJ / R * dtheta + m * g * h / R * sin(alpha) + b_x * dalpha + F_c * tanh(dalpha * 100)

    # def hooke(self, state):
    #     theta, dtheta = state
    #     r0, L0, R, I1, I2, m, h, g = self.params
    #     Kr = self.Kr
    #     Kl = self.Kl
    #
    #     T = self.T(state)
    #     dS = 2 * (r0 / (L0 ** 2 - theta ** 2 * r0 ** 2)) ** 2
    #     dS *= (2 * L0 ** 2 - r0 ** 2 * theta ** 2) / Kr * theta ** 3 + L0 ** 2 * theta / Kl
    #     return 0.5 * dS * T ** 2

    # def J2(self, state):
    #     theta, dtheta = state
    #     r0, L0, R, I1, I2, m, h, g = self.params
    #     Kr = self.Kr
    #     Kl = self.Kl
    #
    #     dS = 2 * (r0 / (L0 ** 2 - theta ** 2 * r0 ** 2)) ** 2
    #     dS *= (2 * L0 ** 2 - r0 ** 2 * theta ** 2) / Kr * theta ** 3 + L0 ** 2 * theta / Kl
    #     return dS

    def J2(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        Kr = self.Kr
        Kl = self.Kl

        dS = 2 * (r0 ** 2 * theta) / (L0 ** 2 - r0 ** 2 * theta ** 2) ** 2
        dS *= (L0 ** 2 / Kl + 2 * L0 ** 2 * theta ** 2 / Kr - r0 ** 2 * theta ** 4 / Kr)
        return dS

    def ddtheta1(self, state, u=0):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        b_theta, tau_c, b_x, F_c = self.fric_params
        mu = self.mu
        J = self.J(state)
        dJ = self.dJ(state)
        J2 = self.J2(state)

        tau = tau_c + b_theta * dtheta  # + TODO: + mu*T и соответственно пересчитать в символьном коде

        sqr = -sign(theta) * R * (I1 ** 2 * R ** 2 + 2 * I1 * I2 ** 2 * J * J2 * dJ * dtheta + 2 * I1 * I2 * J * R * (
                J + J2 * g * m) + I2 ** 2 * J ** 4 +
                                  2 * I2 ** 2 * J ** 2 * J2 * (u - tau)) ** 0.5
        return (sqr - I2 * J ** 2 * R - I2 * J * J2 * R * g * m - I1 * R ** 2 - I2 ** 2 * J * J2 * dJ * dtheta) / (
                I2 ** 2 * J ** 2 * J2)

    def ddtheta2(self, state, u=0):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        b_theta, tau_c, b_x, F_c = self.fric_params
        mu = self.mu
        J = self.J(state)
        T = self.T(state)
        return 1 / (I1 + (J + mu) * I2 * J / R) * (u - (mu + J) * T - b_theta * dtheta - tau_c * tanh(dtheta * 100))

    def ddtheta(self, state, u=0):
        theta, dtheta = state
        return self.ddtheta2(state, u)

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

    def P(self, state, alpha=None):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        if alpha is None:
            alpha = self.alpha(state)
        U = self.U(state)
        return U + m * g * h * (1 - cos(alpha)) - m * g * h * (1 - cos(self.alpha0))

    def U(self, state):
        theta, dtheta = state
        r0, L0, R, I1, I2, m, h, g = self.params
        Kr, Kl = self.Kr, self.Kl
        T0 = self.T((0, 0))
        return self.U_coef * 0.5 / (L0 ** 2 - theta ** 2 * r0 ** 2) * \
               (theta ** 4 * r0 ** 2 / Kr + L0 ** 2 / Kl) * (self.T(state) ** 2 - T0 ** 2)
