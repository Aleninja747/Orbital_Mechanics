import numpy as np
import matplotlib.pyplot as plt
import Perturbations as perts
import OrbitPropagatorV2 as op
import Conversions as conv
from TimeClassV2 import Time


mu_sun = 1.327e11
initial_time = Time('UTC', day=2451545.0)
tspan = 20000
dt = 0.2
perts_dic = perts.perts_dic
perts_dic['TB'] = (True, 'sun', mu_sun)
coes = [6700, 0.07, 69, 0, -14, 0]
mu = 398600.4415
r, v = conv.coes2cart(coes, mu, degree=True)
y0 = np.concatenate([r, v], axis=None)
ys = op.propagate_orbit(y0, initial_time, tspan, dt, perts_dic=perts_dic)
ts = ys[:, 0]
rs = ys[:, 1:4]
vs = ys[:, 4:]

fig = plt.figure(figsize=(16, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(rs[:, 0], rs[:, 1], rs[:, 2])
plt.show()

hs_vec = []
for i in range(len(rs)):
    hs_vec.append(np.cross(rs[0], vs[0])-np.cross(rs[i], vs[i]))

hs_vec = np.array(hs_vec)
fig, ax = plt.subplots()
ax.plot(hs_vec[:, 0], hs_vec[:, 1])
plt.show()
