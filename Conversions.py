import numpy as np
import tools as t
import math as m


def cart2coes(r, v, mu, degrees=False, units='DU', print_results=False):
    singularities = [False, False]
    # first is circular, second is equatorial
    e_vec = t.e_vector(r, v, mu)
    e = np.linalg.norm(e_vec)

    a = t.calc_sma_rv(r, v, mu)

    k = np.array([0, 0, 1])
    h = np.cross(r, v)
    n = np.cross(k, h)
    arg = error_check_angle(np.dot(k, h) / np.linalg.norm(h))
    inclination = np.arccos(arg)

    i = np.array([1, 0, 0])

    if virtually_0(e) and virtually_0(inclination):
        arg = error_check_angle(np.dot(i, r) / np.linalg.norm(r))
        true_lon = np.arccos(arg)
        if r[1] < 0:
            true_lon = 2 * np.pi - true_lon
        singularities[0] = True
        singularities[1] = True
        coes = [a, e, inclination, true_lon]

    elif virtually_0(e):
        arg = error_check_angle(np.dot(n, r) / (np.linalg.norm(n) * np.linalg.norm(r)))
        arg_lat = np.arccos(arg)
        if r[2] < 0:
            arg_lat = 2 * np.pi - arg_lat

        arg = error_check_angle(np.dot(n, i) / np.linalg.norm(n))
        raan = np.arccos(arg)
        if n[1] < 0:
            raan = 2 * np.pi - raan

        singularities[0] = True
        coes = [a, e, inclination, raan, arg_lat]

    elif virtually_0(inclination):
        arg = error_check_angle(np.dot(i, e_vec) / e)
        lon_periapsis = np.arccos(arg)
        if e_vec[1] < 0:
            lon_periapsis = 2 * np.pi - lon_periapsis

        arg = error_check_angle(np.dot(r, e_vec) / (np.linalg.norm(r) * e))
        true_anomaly = np.arccos(arg)
        if np.dot(r, v) < 0:
            true_anomaly = 2 * np.pi - true_anomaly

        singularities[1] = True
        coes = [a, e, inclination, lon_periapsis, true_anomaly]

    else:
        arg = error_check_angle(np.dot(n, i) / np.linalg.norm(n))
        raan = np.arccos(arg)
        if n[1] < 0:
            raan = 2 * np.pi - raan

        arg = error_check_angle(np.dot(e_vec, n) / (e * np.linalg.norm(n)))
        arg_periapsis = np.arccos(arg)
        if e_vec[2] < 0:
            arg_periapsis = 2 * np.pi - arg_periapsis

        arg = error_check_angle(np.dot(r, e_vec) / (np.linalg.norm(r) * e))
        true_anomaly = np.arccos(arg)
        if np.dot(r, v) < 0:
            true_anomaly = 2 * np.pi - true_anomaly

        coes = [a, e, inclination, raan, arg_periapsis, true_anomaly]

    if degrees:
        for j in range(2, len(coes)):
            coes[j] = np.rad2deg(coes[j])
    if print_results:
        print_coes(coes, singularities, units, degrees)
    return coes


def error_check_angle(arg):
    if abs(arg - 1) < 1e-12 and arg > 1:
        arg = 1
    elif abs(arg + 1) < 1e-12 and arg < -1:
        arg = -1

    return arg


def print_coes(coes, singularities, units, degrees=False):
    if degrees:
        angle = 'º'
    else:
        angle = 'rad'

    print('Semi-mayor axis:', coes[0], units)

    print('Eccentricity:', coes[1])
    print('Inclination:', coes[2], angle)
    if singularities[0] and singularities[1]:
        print('True Longitude:', coes[3], angle)
    elif singularities[0]:
        print('Right Ascension of the Ascending Node:', coes[3], angle)
        print('Argument of Latitude:', coes[4], angle)
    elif singularities[1]:
        print('Longitude of Periapsis:', coes[3], angle)
        print('True Anomaly:', coes[4], angle)
    else:
        print('Right Ascension of the Ascending Node:', coes[3], angle)
        print('Argument of Periapsis:', coes[4], angle)
        print('True Anomaly:', coes[5], angle)


def coes2cart(coes, mu, degree=False):
    if degree:
        for j in range(2, len(coes)):
            coes[j] = np.deg2rad(coes[j])
    if virtually_0(coes[1]) and virtually_0(coes[2]):
        p = coes[0]
        r_ijk = p * np.array([np.cos(coes[3]), np.sin(coes[3]), 0])
        v_ijk = np.sqrt(mu / p) * np.array([-np.sin(coes[3]), np.cos(coes[3]), 0])
    elif virtually_0(coes[1]):
        p = coes[0]
        r_pqw = p * np.array([[np.cos(coes[4])], [np.sin(coes[4])], [0]])
        v_pqw = np.sqrt(mu / p) * np.array([[-np.sin(coes[4])], [np.cos(coes[4])], [0]])

        q_pqw_ijk = np.matmul(euler_r3(-coes[3]), euler_r1(-coes[2]))
        for i in range(len(q_pqw_ijk)):
            for j in range(len(q_pqw_ijk[i])):
                if virtually_0(q_pqw_ijk[i][j]):
                    q_pqw_ijk[i][j] = 0
        r_ijk = np.matmul(q_pqw_ijk, r_pqw).T[0]
        v_ijk = np.matmul(q_pqw_ijk, v_pqw).T[0]
    elif virtually_0(coes[2]):
        if coes[1] < 1:
            p = coes[0] * (1 - coes[1] ** 2)
        else:
            p = coes[0] * (coes[1] ** 2 - 1)
        r_pqw = t.r_polar(coes[4], coes[0], coes[1]) * np.array([[np.cos(coes[4])], [np.sin(coes[4])], [0]])
        v_pqw = np.sqrt(mu / p) * np.array([[-np.sin(coes[4])], [coes[1] + np.cos(coes[4])], [0]])
        q_pqw_ijk = euler_r3(coes[3]).T
        for i in range(len(q_pqw_ijk)):
            for j in range(len(q_pqw_ijk[i])):
                if virtually_0(q_pqw_ijk[i][j]):
                    q_pqw_ijk[i][j] = 0
        r_ijk = np.matmul(q_pqw_ijk, r_pqw).T[0]
        v_ijk = np.matmul(q_pqw_ijk, v_pqw).T[0]
    else:
        if coes[1] < 1:
            p = coes[0] * (1 - coes[1] ** 2)
        else:
            p = coes[0] * (coes[1] ** 2 - 1)
        r_pqw = t.r_polar(coes[5], coes[0], coes[1]) * np.array([[np.cos(coes[5])], [np.sin(coes[5])], [0]])
        v_pqw = np.sqrt(mu / p) * np.array([[-np.sin(coes[5])], [coes[1] + np.cos(coes[5])], [0]])
        q_pqw_ijk = np.matmul(np.matmul(euler_r3(-coes[3]), euler_r1(-coes[2])), euler_r3(-coes[4]))
        for i in range(len(q_pqw_ijk)):
            for j in range(len(q_pqw_ijk[i])):
                if virtually_0(q_pqw_ijk[i][j]):
                    q_pqw_ijk[i][j] = 0
        r_ijk = np.matmul(q_pqw_ijk, r_pqw).T[0]
        v_ijk = np.matmul(q_pqw_ijk, v_pqw).T[0]

    return r_ijk, v_ijk


def virtually_0(arg, tol=1e-12):
    if abs(arg) < tol:
        return True
    return False


def euler_r3(angle):
    return np.array([np.array([m.cos(angle), m.sin(angle), 0]), np.array([-m.sin(angle), m.cos(angle), 0]),
                     np.array([0, 0, 1])])


def euler_r2(angle):
    return np.array([np.array([m.cos(angle), 0, -m.sin(angle)]), np.array([0, 1, 0]), np.array([m.sin(angle), 0,
                                                                                                  m.cos(angle)])])


def euler_r1(angle):
    return np.array([np.array([1, 0, 0]), np.array([0, m.cos(angle), m.sin(angle)]), np.array([0, -m.sin(angle),
                                                                                                 m.cos(angle)])])


def PN_mat(x, dx, y, dy, s):
    X = np.deg2rad(x / 3600) + np.deg2rad(dx / 3600)
    Y = np.deg2rad(y / 3600) + np.deg2rad(dy / 3600)
    a = 0.5 + 0.125 * (X ** 2 + Y ** 2)

    pnR1 = np.array([1 - a * X ** 2, -a * X * Y, X])
    pnR2 = np.array([-a * X * Y, 1 - a * Y ** 2, Y])
    pnR3 = np.array([-X, -Y, 1 - a * (X ** 2 + Y ** 2)])

    PN = np.array([pnR1, pnR2, pnR3])
    R3s = euler_r3(np.deg2rad(s / 3600))
    return np.matmul(PN, R3s)


def W_mat(sprime, xp, yp):
    print(sprime/3600, xp/3600, yp/3600)
    r3 = euler_r3(np.deg2rad(-sprime / 3600))
    r2 = euler_r3(np.deg2rad(xp / 3600))
    r1 = euler_r3(np.deg2rad(yp / 3600))

    return np.array(np.matmul(r3, np.matmul(r2, r1)))


def MOD_2_GCCRF(julian_century):
    zeta = np.deg2rad(
        (2306.2181 * julian_century + 0.30188 * julian_century ** 2 + 0.017998 * julian_century ** 3) / 3600)
    theta = np.deg2rad(
        (2004.3109 * julian_century - 0.42665 * julian_century ** 2 - 0.041833 * julian_century ** 3) / 3600)
    z = np.deg2rad(
        (306.2181 * julian_century + 1.09468 * julian_century ** 2 + 0.018203 * julian_century ** 3) / 3600)

    return np.matmul(euler_r3(zeta), np.matmul(euler_r2(-theta), euler_r3(z)))


def J2000_2_GCRF():
    delta_alpha0 = np.deg2rad(.0146/3600)
    xi0 = np.deg2rad(-0.16617/3600)
    eta0 = np.deg2rad(-0.0068192/3600)

    return np.matmul(euler_r3(-delta_alpha0), euler_r3(-xi0), euler_r3(eta0))
