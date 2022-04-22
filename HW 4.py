import numpy as np
import matplotlib.pyplot as plt
import Perturbations as perts
import Conversions as conv
import OrbitPropagatorV2 as op
import tools as t
from Ephemerides import sun
from TimeClassV2 import Time

AU = 149597870.7

rs = np.linspace(6600, 50000, 1000)
moon_pert = []
sun_pert = []
mu_moon = 4.902e3
mu_sun = 1.327e11
r_moon = np.array([384.4e3, 0])
r_sun = np.array([149.6e6, 0])
for r in rs:
    r_sat_moon = r_moon - np.array([r, 0])
    moon_pert.append(np.linalg.norm(perts.third_body_pert(mu_moon, r_sat_moon, r_moon)))
    r_sat_sun = r_sun - np.array([r, 0])
    sun_pert.append(np.linalg.norm(perts.third_body_pert(mu_sun, r_sat_sun, r_sun)))

fig, ax = plt.subplots()
ax.plot(rs, moon_pert, label='Moon')
ax.plot(rs, sun_pert, label='Sun')
ax.set_title('Moon vs. Sun Perturbations at Different Distances')
ax.set_xlabel('Satellite Radius')
ax.set_ylabel('Magnitude of Perturbation [$km/s^2$]')
ax.legend()
plt.show()

# 4
test_time = Time('UTC', day=2451545.0)
r_sun_MOD = AU * sun(test_time.julian())
r_sun_gcrf = np.matmul(conv.MOD_2_GCCRF(test_time.julian_cent()), r_sun_MOD)
print(r_sun_gcrf)
time = Time('UTC', day=2457793.5)
r_sun_MOD = sun(time.julian())
r_sun_gcrf = np.matmul(conv.MOD_2_GCCRF(time.julian_cent()), r_sun_MOD)

print(r_sun_gcrf)

current_time = 0
rs = []
times = []
while current_time <= 365.25:
    r_sun_MOD = sun(time.julian())
    r_sun_gcrf = np.matmul(conv.J2000_2_GCRF().T, np.matmul(conv.MOD_2_GCCRF(time.julian_cent()), r_sun_MOD))
    rs.append(r_sun_gcrf)
    times.append(current_time)
    current_time += 0.25
    time.add_days(0.25)

rs = np.array(rs)
fig, ax = plt.subplots(3, 1)
ax[0].set_title('Sun position in J2000 Frame')
y_labels = ['$r_{GCRF_x}$ [AU]', '$r_{GCRF_y}$ [AU]', '$r_{GCRF_z}$ [AU]']
for i in range(len(ax)):
    ax[i].plot(times, rs[:, i])
    ax[i].set_ylabel(y_labels[i])
    ax[i].set_xlabel('Time [Days]')
plt.show()

# 5
initial_time = Time('UTC', day=2451545.0)
r_tb = sun(initial_time.julian())
r_tb = AU * np.matmul(conv.MOD_2_GCCRF(initial_time.julian_cent()), r_tb)
print('julian century: ', initial_time.julian_cent(), 'r_sun: ', r_tb)
tspan = 365.25*86400
dt = 0.5*86400
perts_dic = perts.perts_dic
perts_dic['TB'] = (True, 'sun', mu_sun)
r = np.array([39066, 221558, 268116])
v = np.array([-1.19828, 0.211289, 0])
print(t.calc_period(t.calc_sma_rv(r, v, 398600), 398600))
print(t.calc_specific_energy(398600, r, v))
y0 = np.concatenate([r, v], axis=None)
print('test: ', op.diff_eq(0, y0, initial_time, 398600, perts_dic=perts_dic))
ys = op.propagate_orbit(y0, initial_time, tspan, dt, 398600.4415, perts_dic=perts_dic)
ts = ys[:, 0]
rs = ys[:, 1:4]
vs = ys[:, 4:]

coes_vec = []
coes = conv.cart2coes(r, v, 398600, True, units='km', print_results=True)
hs_vec = []
for i in range(len(rs)):
    temp_coes = conv.cart2coes(rs[i], vs[i], 1, degrees=True)
    hs_vec.append(np.cross(rs[0], vs[0])-np.cross(rs[i], vs[i]))
    coes_delta = [0] * 6
    for j in range(len(temp_coes)):
        coes_delta[j] =temp_coes[j] #coes[j]
    coes_vec.append(coes_delta)

coes_vec = np.array(coes_vec)
fig, axs = plt.subplots(5, 1)
labels = ["$\Delta$a [DU]", "$\Delta e$", "$\Delta i$", "$\Delta \Omega$", "$\Delta \omega$", "$\Delta \\nu$"]
for i in range(len(axs)):
    axs[i].plot(ts, coes_vec[:, i])
    axs[i].set_ylabel(labels[i])
    axs[i].set_xlabel("Time [TU]")
axs[0].set_title("Change in Orbital Elements vs. Time")

hs_vec = np.array(hs_vec)
fig, ax = plt.subplots()
ax.plot(hs_vec[:, 0], hs_vec[:, 1])
plt.show()
