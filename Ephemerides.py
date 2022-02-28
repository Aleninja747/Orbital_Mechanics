import numpy as np


def sun(time_JD):
    T_ut1 = (time_JD - 51544.5) / 36525
    lambda_sun = 280.460 + 36000.771 * T_ut1  # º
    M_sun = 357.52772333 + 35999.05034 * T_ut1  # º
    lambda_ecliptic = lambda_sun + 1.914666471 * np.sin(np.deg2rad(M_sun)) + 0.019994643 * np.sin(np.deg2rad(2*M_sun))
    r_sun = 1.000140612 - 0.016708617 * np.cos(np.deg2rad(M_sun)) - 0.000139589 * np.cos(np.deg2rad(2*M_sun))
    epsilon = 23.439291 - 0.0130042 * T_ut1
    rx = r_sun * np.cos(np.deg2rad(lambda_ecliptic))
    ry = r_sun * np.sin(np.deg2rad(lambda_sun)) * np.cos(np.deg2rad(epsilon))
    rz = r_sun * np.sin(np.deg2rad(lambda_sun)) * np.sin(np.deg2rad(epsilon))
    r_sun_MOD = np.array([rx, ry, rz])
    return r_sun_MOD
