import numpy as np

perts_dic = {
    'J2': False,
    'TB': (False, 'None', 0),
    'J3': False,
    'Drag': (False, 0, 0),
    'SRP': (False, 0, 0)
}


def J2_Pert(r, mu=398600.4415, r_planet=6300, j2_value=0.0010826267):
    norm_r = np.linalg.norm(r)
    z2 = r[2] ** 2
    r2 = norm_r ** 2
    tx = r[0] / norm_r * (5 * z2 / r2 - 1)
    ty = r[1] / norm_r * (5 * z2 / r2 - 1)
    tz = r[2] / norm_r * (5 * z2 / r2 - 3)
    return 1.5 * j2_value * mu * r_planet ** 2 / norm_r ** 4 * np.array([tx, ty, tz])


def J3_Perts(r, mu=398600.4415, r_planet=6300, j3_value=-2.5327e-6):
    a_j3 = np.zeros(3)
    r_norm = np.linalg.norm(r)
    a_j3[0] = -5 * j3_value * mu * r_planet ** 3 * r[0] / (2 * r_norm ** 7) * (3 * r[2] - 7 * r[2] ** 3 / r_norm ** 2)
    a_j3[1] = -5 * j3_value * mu * r_planet ** 3 * r[1] / (2 * r_norm ** 7) * (3 * r[2] - 7 * r[2] ** 3 / r_norm ** 2)
    a_j3[2] = -5 * j3_value * mu * r_planet ** 3 / (2 * r_norm ** 7) * (6 * r[2] ** 2 - 7 * r[2] ** 4 / r_norm ** 2 - 3 / 5 * r_norm ** 2)

    return a_j3


def third_body_pert(mu_tb, r_sat_tb, r_main_tb):
    a_tb = mu_tb * ((r_sat_tb / (np.linalg.norm(r_sat_tb) ** 3)) - (r_main_tb / (np.linalg.norm(r_main_tb) ** 3)))
    # print('Third Body: ', a_tb)
    return a_tb


def drag(r, v, omega_planet, rho, Cd, A2m):
    omega_planet = np.array([0, 0, omega_planet])
    v_rel = v - np.cross(omega_planet, r)  # in km/s
    v_rel *= 1000  # in m/s
    a_drag = - 1/2 * Cd*A2m * rho * np.linalg.norm(v_rel) * v_rel  # in m/s^2
    a_drag /= 1000  # in km/s^2
    # print('Drag: ', a_drag)
    return a_drag


def stand_atmo(r, radius=True):
    if not np.isscalar(r):
        r = np.linalg.norm(r)
    if radius:
        R_earth = 6300
        altitude = r - R_earth
    else:
        altitude = r
    if altitude > 1000:
        rho0 = 3.019e-15
        h0 = 1000
        H = 268
    elif altitude > 900:
        rho0 = 5.245e-15
        h0 = 900
        H = 181.05
    elif altitude > 800:
        rho0 = 1.170e-14
        h0 = 800
        H = 124.64
    elif altitude > 700:
        rho0 = 3.614e-14
        h0 = 700
        H = 88.667
    elif altitude > 600:
        rho0 = 1.454e-13
        h0 = 600
        H = 71.835
    elif altitude > 500:
        rho0 = 6.967e-13
        h0 = 500
        H = 63.822
    elif altitude > 450:
        rho0 = 1.585e-12
        h0 = 450
        H = 60.828
    elif altitude > 400:
        rho0 = 3.725e-12
        h0 = 400
        H = 58.515
    elif altitude > 350:
        rho0 = 9.518e-12
        h0 = 350
        H = 53.298
    elif altitude > 300:
        rho0 = 2.418e-11
        h0 = 300
        H = 53.628
    elif altitude > 250:
        rho0 = 7.248e-11
        h0 = 250
        H = 45.546
    elif altitude > 200:
        rho0 = 2.789e-10
        h0 = 200
        H = 37.105
    elif altitude > 180:
        rho0 = 5.464e-10
        h0 = 180
        H = 29.740
    elif altitude > 150:
        rho0 = 2.070e-9
        h0 = 150
        H = 22.523
    elif altitude > 140:
        rho0 = 3.845e-9
        h0 = 140
        H = 16.149
    elif altitude > 130:
        rho0 = 8.484e-9
        h0 = 130
        H = 12.636
    elif altitude > 120:
        rho0 = 2.438e-8
        h0 = 120
        H = 9.473
    elif altitude > 110:
        rho0 = 9.661e-8
        h0 = 110
        H = 7.263
    elif altitude > 100:
        rho0 = 5.297e-7
        h0 = 100
        H = 5.877
    elif altitude > 90:
        rho0 = 3.396e-6
        h0 = 90
        H = 5.382
    elif altitude > 80:
        rho0 = 1.905e-5
        h0 = 80
        H = 5.799
    elif altitude > 70:
        rho0 = 8.770e-5
        h0 = 70
        H = 6.549
    elif altitude > 60:
        rho0 = 3.206e-4
        h0 = 60
        H = 7.714
    elif altitude > 50:
        rho0 = 1.057e-3
        h0 = 50
        H = 8.382
    elif altitude > 40:
        rho0 = 3.972e-3
        h0 = 40
        H = 7.554
    elif altitude > 30:
        rho0 = 1.774e-2
        h0 = 30
        H = 6.682
    elif altitude > 25:
        rho0 = 3.899e-2
        h0 = 25
        H = 6.349
    else:
        rho0 = 1.225
        h0 = 0
        H = 7.249

    return rho0*np.exp(-(altitude-h0)/H)


def srp(r, r_sun, C_r, A2m):
    SF = 1367  # W/m^2
    c = 299792458  # m/s
    psrp = SF/c  # kg s^2/m
    r_ss = (r_sun - r) * 1000  # m
    r_norm = np.linalg.norm(r_ss)  # m
    a_srp = (-psrp * C_r * A2m) * r_ss / r_norm  # m/s^2
    # print('SRP: ', a_srp/1000)
    return a_srp/1000
