import numpy as np
import matplotlib.pyplot as plt
import Perturbations as perts
import Conversions as conv
import OrbitPropagatorV2 as op
import tools as t
from Ephemerides import sun
from TimeClassV2 import Time

# 4)
alts = np.linspace(0, 1000, num=10000)
rhos = []
for alt in alts:
    rhos.append(perts.stand_atmo(alt, radius=False))

fig, ax = plt.subplots()
ax.plot(alts, rhos)
ax.set_xscale('log')
plt.show()
