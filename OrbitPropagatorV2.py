import numpy as np
import Perturbations as perts
import Ephemerides as eph
import Conversions as conv
from scipy.integrate import ode

no_perts = perts.perts_dic
j2_value = 0.00108248
R = 6378.1363
Omega_earth = 7.2921e-5
AU = 149597870.7


def diff_eq(t, y, initial_time, mu, perts_dic):
    # unpack the state
    r = y[:3]
    v = y[3:]

    # Initialize the dY array
    dy = np.zeros(y.shape)

    # calculate acceleration vector
    a = -mu/(np.linalg.norm(r)**3) * r

    #  Include J2 perturbation
    if perts_dic['J2']:
        a += perts.J2_Pert(r, mu, R, j2_value)
    if perts_dic['J3']:
        a += 1
    if perts_dic['TB'][0]:
        r_tb = 0
        if perts_dic['TB'][1] == 'sun':
            r_tb = eph.sun(initial_time.julian())
            r_tb = AU * np.matmul(conv.MOD_2_GCCRF(initial_time.julian_cent()), r_tb)
        r_sat = r_tb-r
        a += perts.third_body_pert(perts_dic['TB'][2], r_sat, r_tb)
    if perts_dic['Drag'][0]:
        rho = perts.stand_atmo(r)
        a += perts.drag(r, v, Omega_earth, rho, perts_dic['Drag'][1], perts_dic['Drag'][2], perts_dic['Drag'][3])
    dy[:3] = v
    dy[3:] = a

    return dy


def propagate_orbit(y0, initial_time, tspan, dt, mu=398600.4415, perts_dic=no_perts):
    solver = ode(diff_eq).set_f_params(initial_time, mu, perts_dic)
    solver.set_integrator('dopri5', rtol=1e-15, atol=1e-20)
    solver.set_initial_value(y0, 0)

    y = [np.insert(y0, 0, 0)]
    while solver.successful() and solver.t < tspan:
        if abs(tspan-solver.t) < dt:
            dt = abs(tspan-solver.t)
        solver.integrate(solver.t + dt)
        initial_time.add_days(dt/86400)
        y.append(np.insert(solver.y, 0, solver.t))

    y = np.array(y)

    return y
