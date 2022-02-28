import numpy as np
import tools as t
from OrbitPropagator import OrbitPropagator as op
import planetary_data as pd

# Problem 2
r1 = np.array([4000, 4000, 4000])
v1 = np.array([5, 5, 5])
y = r1.tolist() + v1.tolist()
dt = 200  # seconds
cb = pd.earth
tspan = 100
orbit_propagator1 = op(r1, v1, tspan, dt)
dY = orbit_propagator1.dif_eq(0, y)
print('2')
print(dY)
print()

# 2.1
r = np.array([6667.998, -831.435, -1531.486])
v = np.array([-1.210185, -7.384715, -1.138699])
y = r.tolist() + v.tolist()
tspan = 3 * t.calc_period(t.calc_sma_rv(r, v, cb['mu']), cb['mu'])
orbit_propagator = op(r, v, tspan, dt)
dY = orbit_propagator.dif_eq(0, y)
print('2.1')
print(dY)

# Problem 3.1
orbit_propagator.propagate_orbit()
r = orbit_propagator.rs[len(orbit_propagator.rs)-1]
v = orbit_propagator.vs[len(orbit_propagator.vs)-1]
print(r, '\n', v, '\n')
