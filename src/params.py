from copy import copy

# lin_params = {
#     'r0': 0.74 * 1e-3,
#     "L0": 160 * 1e-3,
#     "R": 50 * 1e-3,
#     "I1": 9.48e-7,
#     "I2": None,
#     "m": 2.0,
#     "h": 110 * 1e-3,
#     "g": 9.81,
#     "b_theta": 1.31 * 1e-6,
#     "tau_c": 1.91 * 1e-3,
#     "b_x": 9.46,
#     "F_c": 4.11,
# }

lin_params = {
    'r0': 0.74 * 1e-3,
    "L0": 160 * 1e-3,
    "R": None,
    "I1": 3.848e-6,
    "I2": None,
    "m": 2.150,
    "h": None,
    "g": 9.8,
    "b_theta": 1.31 * 1e-6,
    "tau_c": 1.91 * 1e-3,
    "b_x": 9.46,
    "F_c": 4.11,
}

lin_params_no_fric = copy(lin_params)
lin_params_no_fric['b_theta'] = 0.0
lin_params_no_fric['tau_c'] = 0.0
lin_params_no_fric['b_x'] = 0.0
lin_params_no_fric['F_c'] = 0.0

# rot_params = {
#     'r0': 0.74 * 1e-3,
#     "L0": 350 * 1e-3,
#     "R": 50 * 1e-3,
#     "I1": 0.12 * 1e-3,
#     "I2": 0.15,
#     # "m": 1.3,  # 2.5,
#     "m": 2.6,  # 2.5,
#     "h": 235 * 1e-3,
#     "g": 9.81,
#     "b_theta": 4e-2,  # 2.58 * 1e-6,
#     "tau_c": 0.75,  # * 1e-3,
#     "b_x": 9.46,  # ??
#     "F_c": 4.11,  # ??
#     "mu": 0.31,
#     "Kr": 11750 / (0.74 * 1e-3),
#     "Kl": 9980 / (350 * 1e-3),
# }

rot_params = {
    'r0': 0.74 * 1e-3,
    "L0": 310 * 1e-3,
    "R": 46 * 1e-3,
    # "R": 50 * 1e-3,
    "I1": 0.095 * 1e-3,
    "I2": 0.15,
    # "m": 1.3,  # 2.5,
    "m": 2.6,  # 2.5,
    "h": 235 * 1e-3,
    "g": 9.81,
    "b_theta": 4e-2,  # 2.58 * 1e-6,
    "tau_c": 0.75,  # * 1e-3,
    "b_x": 9.46,  # ??
    "F_c": 4.11,  # ??
    "mu": 0.31,
    "Kr": 11750,
    "Kl": 9980,
}

rot_params['Kr'] /= rot_params['r0']
rot_params['Kl'] /= rot_params['L0']

rot_params_no_fric = copy(rot_params)
rot_params_no_fric["b_theta"] = 0
rot_params_no_fric["tau_c"] = 0
rot_params_no_fric["b_x"] = 0
rot_params_no_fric["F_c"] = 0
rot_params_no_fric["mu"] = 0
# rot_params_no_fric["Kr"] = float('inf')
# rot_params_no_fric["Kl"] = float('inf')
