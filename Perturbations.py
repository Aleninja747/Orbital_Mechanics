import numpy as np

perts_dic = {
    'J2': False,
    'TB': (False, 'None', 0)
}


def J2_Pert(r, mu, r_planet=6378.1363, j2_value=0.00108248):
    norm_r = np.linalg.norm(r)
    z2 = r[2] ** 2
    r2 = norm_r ** 2
    tx = r[0] / norm_r * (5 * z2 / r2 - 1)
    ty = r[1] / norm_r * (5 * z2 / r2 - 1)
    tz = r[2] / norm_r * (5 * z2 / r2 - 3)
    return 1.5 * j2_value * mu * r_planet ** 2 / norm_r ** 4 * np.array([tx, ty, tz])


def third_body_pert(mu_tb, r_sat_tb, r_main_tb):
    return mu_tb * ((r_sat_tb / (np.linalg.norm(r_sat_tb)**3)) - (r_main_tb / (np.linalg.norm(r_main_tb)**3)))
