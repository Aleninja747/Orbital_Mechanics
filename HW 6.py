import numpy as np
import matplotlib.pyplot as plt
import Perturbations as perts
import Conversions as conv
import OrbitPropagatorV2 as op
from TimeClassV2 import Time

mu_sun = 1.327e11
Cd = 2
A2m = 0.1
C_r = 1.5
mu = 398600

# # 2)
# r0 = np.array([100, 5000, -6000])
# ap = perts.J2_Pert(r0) + perts.J3_Perts(r0)
# print('======2=======')
# print('a_p: ', ap)
# 3)

r = np.array([2781, 5318, -5629])
r_sun = np.array([2379260.000, 148079334.000, -1009936.000])
a_srp = perts.srp(r, r_sun, C_r, A2m)
print('======3=======')
print('a_srp: ', a_srp)
# # 4)
# alts = np.linspace(0, 1000, num=10000)
# rhos = []
# for alt in alts:
#     rhos.append(perts.stand_atmo(alt, radius=False))
#
# fig, ax = plt.subplots()
# ax.plot(alts, rhos)
# ax.set_yscale('log')
# ax.set_ylabel('Density [kg/m^3]')
# ax.set_xlabel('Height [km]')
# plt.show()

# 5)
coes = [6800, 0.005, 71, 300, 78, 0]
r0, v0 = conv.coes2cart(coes, degree=True)
y0 = np.concatenate([r0, v0], axis=None)
initial_time = Time('UTC', day=2458200.5)
time = initial_time.julian()
perts_dic = perts.perts_dic
perts_dic['J2'] = True
perts_dic['J3'] = True
perts_dic['TB'] = (True, 'sun', mu_sun)
perts_dic['Drag'] = (True, Cd, A2m)
perts_dic['SRP'] = (True, C_r, A2m)
tspan = 24*3600
dt = 60

print('initial value:', op.diff_eq(0, y0, initial_time, mu, perts_dic))
ys1 = op.propagate_orbit(y0, initial_time, tspan, dt)
ts1 = ys1[:, 0]
rs1 = ys1[:, 1:4]
vs1 = ys1[:, 4:]
ys2 = op.propagate_orbit(y0, initial_time, tspan, dt, perts_dic=perts_dic)
ts2 = ys2[:, 0]
rs2 = ys2[:, 1:4]
vs2 = ys2[:, 4:]

print('======5.1=======')
print('Keplerian:')
print(' r: ', rs1[-1])
print(' v: ', vs1[-1])
print('High Fidelity:')
print(' r: ', rs2[-1])
print(' v: ', vs2[-1])

r_diff = []
v_diff = []

for i in range(len(rs1)):
    r_diff.append(rs1[i, :]-rs2[i, :])
    v_diff.append(vs1[i, :]-vs2[i, :])

r_diff_1 = np.array(r_diff)
v_diff_1 = np.array(v_diff)

fig, axs = plt.subplots(3, 1)
for i in range(3):
    axs[i].plot(ts1, r_diff_1[:, i])

plt.show()
