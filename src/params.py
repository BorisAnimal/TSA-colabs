from copy import copy

lin_params = {
    'r0': 0.74 * 1e-3,
    "L0": 160 * 1e-3,
    "R": 50 * 1e-3,
    "I1": 9.48e-7,
    "I2": None,
    "m": 2.0,
    "h": 110 * 1e-3,
    "g": 9.81,
    "b_theta": 1.31 * 1e-6,
    "tau_c": 1.91 * 1e-3,
    "b_x": 9.46,
    "F_c": 4.11,
}

lin_params_no_fric = {
    'r0': 0.74 * 1e-3,
    "L0": 160 * 1e-3,
    "R": 50 * 1e-3,
    "I1": 9.48e-7,
    "I2": None,
    "m": 2.0,
    "h": 110 * 1e-3,
    "g": 9.81,
    "b_theta": 0,
    "tau_c": 0,
    "b_x": 0,
    "F_c": 0,
}

rot_params = {
    'r0': 0.74 * 1e-3,
    "L0": 350 * 1e-3,
    "R": 50 * 1e-3,
    "I1": 0.11995 * 1e-3,
    "I2": 0.1477,
    "m": 1.3, #2.5,
    "h": 240 * 1e-3,
    "g": 9.81,
    "b_theta": 4e-2, #2.58 * 1e-6,
    "tau_c": 0.75, # * 1e-3,
    "b_x": 9.46,  # ??
    "F_c": 4.11,  # ??
    "mu": 0.31,
}

rot_params_no_fric = copy(rot_params)
rot_params_no_fric["b_theta"] = 0
rot_params_no_fric["tau_c"] = 0
rot_params_no_fric["b_x"] = 0
rot_params_no_fric["F_c"] = 0
