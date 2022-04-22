import numpy as np
import matplotlib.pyplot as plt
import Perturbations as perts
import Conversions as conv
import OrbitPropagatorV2 as op
import tools as t
from Ephemerides import sun
from TimeClassV2 import Time

AU = 149597870.7
mu_sun = 1.327e11

# 5
initial_time = Time('UTC', day=2451546.0)
r_tb = sun(initial_time.julian())
r_tb = AU * np.matmul(conv.MOD_2_GCCRF(initial_time.julian_cent()), r_tb)
print('julian century: ', initial_time.julian_cent(), 'r_sun: ', r_tb)
tspan = 365.25*86400
dt = 0.5*86400
perts_dic = perts.perts_dic
perts_dic['TB'] = (True, 'sun', mu_sun)
r = np.array([39066, 221558, 268116])
v = np.array([-1.19828, 0.211289, 0])
print(np.isscalar(r))
print(t.calc_period(t.calc_sma_rv(r, v, 398600), 398600))
print(t.calc_specific_energy(398600, r, v))
y0 = np.concatenate([r, v], axis=None)
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
        coes_delta[j] = coes[j] - temp_coes[j]
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