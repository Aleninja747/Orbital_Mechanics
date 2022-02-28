import OrbitPropagatorV2 as op
import tools as t
import numpy as np
import matplotlib.pyplot as plt
import Conversions as conv

# Problem 2
r1 = np.array([4600, 3250, 3250])
v1 = np.array([-5.537, 3.915, 3.9154])
mu = 398600.4415
y1 = np.concatenate((r1, v1), axis=None)
dy1 = op.diff_eq(0, y1, mu)
print('Problem 2')
print(dy1, '\n')

# 2.1
r = np.array([6667.998, -831.435, -1531.486])
v = np.array([-1.210185, -7.384715, -1.138699])
y = np.concatenate((r, v), axis=None)
dy = op.diff_eq(0, y, mu)
print('2.1')
print(dy, '\n')

# Problem 3
energy0 = t.calc_specific_energy(mu, r, v)
tspan = 3 * t.calc_period(t.calc_sma_rv(r, v, mu), mu)
y_result = op.propagate_orbit(y, tspan, 20, mu)
print('3.1')
print('Final Position: ', y_result[-1][1:4], 'km')
print('Final Velocity: ', y_result[-1][4:], 'km/s\n')
ts = []
rs = []
vs = []
acs = []
energies = []

for i in range(len(y_result)):
    ts.append(y_result[i][0]/3600)
    rs.append(np.linalg.norm(y_result[i][1:4]))
    vs.append(np.linalg.norm(y_result[i][4:]))
    energies.append(energy0 - t.calc_specific_energy(mu, rs[i], vs[i]))
    acs.append(-mu/rs[i]**2)

fig, axs = plt.subplots(3, 1)
axs[0].plot(ts, rs, color='#ff9d52')
axs[1].plot(ts, vs, color='#2397fc')
axs[2].plot(ts, acs, color='#40ff50')
axs[0].set_title("Magnitude of Radius, Velocity and Acceleration vs.Time")
axs[2].set_xlabel('Time [h]')
axs[0].set_ylabel('Radius [km]')
axs[1].set_ylabel('Velocity [km/s]')
axs[2].set_ylabel('Acceleration [km/s^2]')
plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(ts, energies, color="#ed2c09")
ax2.set_title('Specific Energy change vs. time')
ax2.set_xlabel('Time [h]')
ax2.set_ylabel('Specific Energy Change [km^2/s^2]')
plt.show()

# Problem 4
dy1 = op.diff_eq(0, y1, mu, j2=True)
print('Problem 4')
print(dy1, '\n')

# Problem 4.1
dy = op.diff_eq(0, y, mu, j2=True)
print('4.1')
print(dy, '\n')

# Problem 5
j2_value = 0.00108248
R = 6378.1363
energy0 = t.calc_specific_energy(mu, r, v,j2=True)
coes = conv.cart2coes(r, v, mu)
omega_dot = -((3/2) * ((np.sqrt(mu)*j2_value * R**2) / ((1-coes[1]**2)**2 * coes[0]**(7/2)))*np.cos(coes[2]))
print('5.1')
tspan = t.calc_period(t.calc_sma_rv(r, v, mu), mu) * 4
y_result = op.propagate_orbit(y, tspan, 100, mu, j2=True)

print('Final Position: ', y_result[-1][1:4], 'km')
print('Final Velocity: ', y_result[-1][4:], 'km/s\n')

ts = []
a = []
omegas = []
energies = []
period = t.calc_period(t.calc_sma_rv(r, v, mu), mu)

for i in range(len(y_result)):
    ts.append(y_result[i][0] / period)
    energies.append(energy0 - t.calc_specific_energy(mu, y_result[i][1:4], y_result[i][4:], j2=True))
    a.append(conv.cart2coes(y_result[i][1:4], y_result[i][4:], mu)[0])
    omegas.append(conv.cart2coes(y_result[i][1:4], y_result[i][4:], mu)[3])


fig4, ax4 = plt.subplots()
ax4.plot(ts, a, color="#931ed6")
ax4.set_title('Semi-mayor Axis vs. Time')
ax4.set_xlabel('time [h]')
ax4.set_ylabel('Semi-mayor Axis [km]')
plt.show()

fig5, ax5 = plt.subplots()
ax5.plot(ts, omegas, color="#16f77f")
ax5.set_title('RAAN vs Time')
ax5.set_xlabel('Time [h]')
ax5.set_ylabel('RAAN [rad]')
plt.show()

print('5.4')
print('omega_dot: ', omega_dot, 'rad/s')
print('omega_dot * t = ', omega_dot*tspan, 'rad')
print('change in omega = ', abs(omegas[-1]-omegas[0]), 'rad')

fig6, ax6 = plt.subplots()
ax6.plot(ts, energies, color="#1656f7")
ax6.set_title('Specific Energy change vs. time')
ax6.set_xlabel('Time [h]')
ax6.set_ylabel('Specific Energy Change [km^2/s^2]')
plt.show()
