import Conversions as conv
import OrbitPropagatorV2 as op
import numpy as np
import matplotlib.pyplot as plt
import tools as t
import Perturbations as perts
from TimeClassV2 import Time

r = np.array([-2.2, -1.8, 2.2])
v = np.array([0.22, -0.528, -0.22])

# Problem 4
coes = conv.cart2coes(r, v, 1, True, print_results=True)

# Problem 5
y0 = np.concatenate([r, v], axis=None)
init_time = Time('UTC', day=2451545.0)
ys = op.propagate_orbit(y0, init_time, 2000, 0.2, 1)
ts = ys[:, 0]
rs = ys[:, 1:4]
vs = ys[:, 4:]
fig, axs = plt.subplots(3, 1)
labels = ["X", "Y", "Z"]
for i in range(len(axs)):
    axs[i].plot(ts, rs[:, i])
    axs[i].set_xlabel('Time [TU]')
    axs[i].set_ylabel(labels[i] + " Position [DU]")
axs[0].set_title("Position Components vs. Time")

fig = plt.figure(figsize=(16, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(rs[:, 0], rs[:, 1], rs[:, 2])
plt.show()

# Problem 6
e0 = t.calc_specific_energy(1, rs[0], vs[0])
es = []
for i in range(len(rs)):
    es.append(e0 - t.calc_specific_energy(1, rs[i], vs[i]))
fig, ax = plt.subplots()
ax.plot(ts, es)
ax.set_xlabel("Time [TU]")
ax.set_ylabel("$\Delta \epsilon$ $[DU^2/ TU^2]$")
ax.set_title("Energy Change vs. Time")
plt.show()

# Problem 7
coes_vec = []
for i in range(len(rs)):
    temp_coes = conv.cart2coes(rs[i], vs[i], 1, degrees=True)
    coes_delta = [0] * 6
    for j in range(len(temp_coes)):
        coes_delta[j] = coes[j] - temp_coes[j]
    coes_vec.append(coes_delta)

coes_vec = np.array(coes_vec)
fig, axs = plt.subplots(6, 1)
labels = ["$\Delta$a [DU]", "$\Delta e$", "$\Delta i$", "$\Delta \Omega$", "$\Delta \omega$", "$\Delta \\nu$"]
for i in range(len(axs)):
    axs[i].plot(ts, coes_vec[:, i])
    axs[i].set_ylabel(labels[i])
    axs[i].set_xlabel("Time [TU]")
axs[0].set_title("Change in Orbital Elements vs. Time")
plt.show()
