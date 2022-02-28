import numpy as np


def calc_period(a, mu):
    return 2 * np.pi * np.sqrt(a ** 3 / mu)


def polar_to_xy(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return x, y


def r_polar(theta, a, e):
    if e < 1:
        return (a * (1 - e ** 2)) / (1 + e * np.cos(theta))
    return (a * (e ** 2 - 1)) / (1 + e * np.cos(theta))


def vis_viva(r, a, mu, e):
    if e < 1:
        return np.sqrt(mu * (2 / r - 1 / a))
    return np.sqrt(mu * (2 / r + 1 / a))


def calc_specific_energy(mu, r, v, j2=False, j2_value=0.00108248, radius=6378.1363):
    if j2:
        return np.dot(v, v) / 2 - mu / np.linalg.norm(r) + mu / np.linalg.norm(r) * j2_value / 2 * \
               (radius / np.linalg.norm(r)) ** 2 * (3 * (r[2] / np.linalg.norm(r)) ** 2 - 1)
    return np.dot(v, v) / 2 - mu / np.linalg.norm(r)


def calc_spec_potential(mu, r):
    return -mu / np.linalg.norm(r)


def calc_spec_kinetic(v):
    return (np.linalg.norm(v) ** 2) / 2


def e_vector(r, v, mu=1):
    return (np.cross(v, np.cross(r, v))) / mu - r / np.linalg.norm(r)


def calc_true_anomaly(r, v, h, e, mu=1):
    true_anomaly = np.arccos((np.linalg.norm(h) ** 2 / mu) * (1 / (np.linalg.norm(r) * e)) - (1 / e))
    if np.dot(r, v) < 0 and e < 1:
        true_anomaly = 2 * np.pi - true_anomaly
    return true_anomaly


def calc_sma_rv(r, v, mu):
    e = np.linalg.norm(e_vector(r, v, mu))
    energy = calc_specific_energy(mu, r, v)
    if e < 1:
        return - mu / (2 * energy)
    return mu / (2 * energy)


def calc_sma(r, e, true_anomaly):
    if e < 1:
        return np.linalg.norm(r) * (1 + e * np.cos(true_anomaly)) / (1 - e ** 2)
    return np.linalg.norm(r) * (1 + e * np.cos(true_anomaly)) / ((e ** 2) - 1)


def calc_ra_rp(a, e):
    if e < 1:
        rp = a * (1 - e)
        ra = a * (1 + e)
    else:
        rp = a * (e - 1)
        ra = -a * (e + 1)
    return rp, ra


def calc_eccentric_anomaly(true_anomaly, e):
    if e < 1:
        return 2 * np.arctan2(np.sqrt(1 - e) * np.tan(true_anomaly / 2), np.sqrt(1 + e))
    return 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(true_anomaly / 2))


def calc_mean_anomaly(eccentric_anomaly, e):
    if e < 1:
        return eccentric_anomaly - e * np.sin(eccentric_anomaly)
    return - eccentric_anomaly + e * np.sinh(eccentric_anomaly)


def eccentric2true(eccentric_anomaly, e):
    if e < 1:
        return 2 * np.arctan2(np.sqrt(1 + e) * np.tan(eccentric_anomaly / 2), np.sqrt(1 - e))
    return 2 * np.arctan2(np.sqrt(e + 1) * np.tanh(eccentric_anomaly / 2), np.sqrt(e - 1))


def newton_true_anomaly(mean_anomaly, e, error=1e-10, max_iter=200):
    if mean_anomaly < np.pi:
        ecc_anomaly_0 = mean_anomaly + e / 2
    else:
        ecc_anomaly_0 = mean_anomaly - e / 2
    previous = ecc_anomaly_0
    current = previous - (previous - e * np.sin(previous) - mean_anomaly) / (1 - e * np.cos(previous))
    iteration = 0
    while abs(current - previous) > error and iteration < max_iter:
        iteration += 1
        previous = current
        current = previous - (previous - e * np.sin(previous) - mean_anomaly) / (1 - e * np.cos(previous))

    true_anomaly = eccentric2true(current, e)
    return true_anomaly


def newton_true_anomaly_h(mean_anomaly_h, e, error=1e-10, max_iter=1000):
    previous = mean_anomaly_h
    current = previous - ((-previous + e * np.sinh(previous) - mean_anomaly_h) / (-1 + e * np.cosh(previous)))
    iteration = 0
    while abs(current - previous) > error and iteration < max_iter:
        iteration += 1
        previous = current
        current = previous - (-previous + e * np.sinh(previous) - mean_anomaly_h) / (-1 + e * np.cosh(previous))
    true_anomaly = eccentric2true(current, e)
    return true_anomaly


def calendar2julian(year, month, day, hour, minute, second):
    return 367 * year - int((7 * (year + int((month + 9) / 12))) / 4) + int(275 * month / 9) + day + 1721013.5 + 1 / 24 \
           * (hour + (1 / 60) * (minute + second / 60))


def theta_gmst(T):
    theta = 67310.54841 + (876600 * 3600 + 8640184.812866) * T + (0.093104 * T ** 2) - (6.2e-6 * T ** 3)
    theta = theta - int(theta / 86400) * 86400
    theta = theta / 240
    return theta


def julian_cent(mjd):
    return (mjd - 51544.5) / 36525


def ERA(mjd):
    return 2 * np.pi * (0.7790572732640 + 1.00273781191135448 * (mjd + 2400000.5 - 2451545))



