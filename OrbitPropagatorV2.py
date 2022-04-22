import numpy as np
import Perturbations as perts
import Ephemerides as eph
import Conversions as conv
from scipy.integrate import ode

no_perts = perts.perts_dic
j2_value = 0.00108248
R = 6300
Omega_earth = 7.2921e-5
AU = 149597870


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
        a += perts.J2_Pert(r, mu)
    if perts_dic['J3']:
        a += perts.J3_Perts(r, mu)
        print('J2+J3: ', a+(mu/(np.linalg.norm(r)**3) * r))
    if perts_dic['TB'][0]:
        r_tb = 0
        if perts_dic['TB'][1] == 'sun':
            r_tb = eph.sun(initial_time.julian())
            r_tb = AU * np.matmul(conv.MOD_2_GCCRF(initial_time.julian_cent()), r_tb)
            #print('third body (sun):', r_tb)
        r_sat = r_tb-r
        #print('r_sat: ', r_sat)
        #print("Sun", perts.third_body_pert(perts_dic['TB'][2], r_sat, r_tb))
        a += perts.third_body_pert(perts_dic['TB'][2], r_sat, r_tb)
    if perts_dic['Drag'][0]:
        rho = perts.stand_atmo(r)
        print("drag", perts.drag(r, v, Omega_earth, rho, perts_dic['Drag'][1], perts_dic['Drag'][2]))
        a += perts.drag(r, v, Omega_earth, rho, perts_dic['Drag'][1], perts_dic['Drag'][2])
    if perts_dic['SRP'][0]:
        r_sun = eph.sun(initial_time.julian())
        print('SRP', perts.srp(r, r_sun, perts_dic['SRP'][1], perts_dic['SRP'][2]))
        a += perts.srp(r, r_sun, perts_dic['SRP'][1], perts_dic['SRP'][2])
    dy[:3] = v
    dy[3:] = a

    return dy


def propagate_orbit(y0, initial_time, tspan, dt, mu=398600.4415, perts_dic=no_perts):
    solver = ode(diff_eq).set_f_params(initial_time, mu, perts_dic)
    solver.set_integrator('dopri5', rtol=1e-17, atol=1e-25)
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
