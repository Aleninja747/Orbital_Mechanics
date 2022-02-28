import numpy as np
import Groundtrack as gt
import tools as t
import Conversions as conv

# 4.1
coes = [6700, 0.07, 69, 0, -14, 0]
mu = 398600.4415
r, v = conv.coes2cart(coes, mu, degree=True)
y0 = np.concatenate([r, v], axis=None)
tspan = 3 * t.calc_period(coes[0], mu)
dt = 60
omega_earth = 2*np.pi/86164
latitudes, longitudes = gt.generate_groundtracks(y0, 0, omega_earth, tspan, dt)
gt.plot_groundtrack(latitudes, longitudes)
long_dif = abs(longitudes[-1]-longitudes[0])/3
period = np.deg2rad(long_dif)/omega_earth
print('4.1')
print('Estimated Period:', period, 's')
print('Period from sma:', tspan/3, 's\n')

# 4.2
coes = [50000, 0.01, 139, 0, -63, 0]
r, v = conv.coes2cart(coes, mu, degree=True)
y0 = np.concatenate([r, v], axis=None)
tspan = 2 * t.calc_period(coes[0], mu)
latitudes, longitudes = gt.generate_groundtracks(y0, 0, omega_earth, tspan, dt)
gt.plot_groundtrack(latitudes, longitudes)
long_dif = abs(longitudes[-1]-longitudes[240])
period = (np.deg2rad(long_dif)+np.pi*2)/omega_earth
print('4.2')
print('Estimated Period:', period, 's')
print('Period from sma:', tspan/2, 's\n')

# 4.3
a = gt.calc_sma_rgt(6, 5)
coes = [a, 0, 71, 194, 0]
r, v = conv.coes2cart(coes, mu, degree=True)
y0 = np.concatenate([r, v], axis=None)
tspan = 7 * t.calc_period(coes[0], mu)
latitudes, longitudes = gt.generate_groundtracks(y0, 0, omega_earth, tspan, dt)
gt.plot_groundtrack(latitudes, longitudes)
print('4.3')
print('The semi-mayor axis has a value of :', a, 'km')
