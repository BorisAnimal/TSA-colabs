from math import sqrt

from scipy import integrate


def T_rot(theta0, r0, L0, m, I, g=9.81):  # R,L,m,I,
    def f(theta):
        J = theta * r0 ** 2 / sqrt(L0 ** 2 - theta ** 2 * r0 ** 2)  # TODO: check energy members
        D = I + m * J ** 2
        E0 = m * g * (L0 - sqrt(L0 ** 2 - theta0 ** 2 * r0 ** 2))
        P = m * g * (L0 - sqrt(L0 ** 2 - theta ** 2 * r0 ** 2))
        return sqrt(D / (E0 - P))

    quad = integrate.quad(f, 0, theta0)[0]
    return 4 * quad / sqrt(2)

# Deprecated
# def T_lin(theta_amp, r0, L0, m, I, g=9.81):
#     def f(theta):
#         D = I + m * theta ** 2 * r0 ** 2 / (L0 ** 2 - theta ** 2 * r0 ** 2)
#         E0 = m * g * (L0 - sqrt(L0 ** 2 - theta_amp ** 2 * r0 ** 2))
#         P = m * g * (L0 - sqrt(L0 ** 2 - theta ** 2 * r0 ** 2))
#         return sqrt(D / (E0 - P))
#
#     quad = integrate.quad(f, 0, theta_amp)[0]
#     return 4 * quad / sqrt(2)
