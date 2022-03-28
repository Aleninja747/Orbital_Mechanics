import numpy as np
import tools as t
import Conversions as conv
from TimeClassV2 import Time

Date = Time('UTC', year=2022, month=2, day=19, hour=5, minute=59, second=0)
dut1 = -0.110914
xp = 0.016768
yp = 0.357565
dX = 0.000066
dY = -0.000049
X = 438.4313
Y = 4.391260
s = -0.006758
leaps = 37

Date.to_ut1()
mjd_ut1 = Date.julian()
Date.to_tt()


s_prime = -0.000047 * Date.julian_cent()
s_i = xp/3600
print(conv.W_mat(s_prime, xp, yp)[0, 1])

print('1.1) JD_ut1:', mjd_ut1 + 2400000.5, '\n')

era = t.ERA(mjd_ut1)
era = era - (era//(2*np.pi))*2*np.pi
print('1.2)ERA:', era, '[rad]\n')

W_mat = conv.W_mat(s_prime, xp, yp)
print('1.3)\nW Matrix: \n', W_mat)

R_mat = conv.euler_r3(-era)
print('Earth Rotation Matrix: \n', R_mat)

PN = conv.PN_mat(X, dX, Y, dY, s)

print('PN Matrix: \n', PN)

Q_itrf_gcrf = np.matmul(PN, np.matmul(R_mat, W_mat))
print('1.4) \n Q ITRF to GCRF: \n', Q_itrf_gcrf)

# Problem 2
T_ut1 = t.julian_cent(mjd_ut1)
r_itrf = np.array([-742.845, -5463.244, 3196.066])
theta_gmst = t.theta_gmst(T_ut1)
print('\n2.1) Theta GMST: ', theta_gmst, '[deg]\n')

r3_gmst = conv.euler_r3(-np.deg2rad(theta_gmst))
r_gcrf_r3 = np.matmul(r3_gmst, r_itrf)
print('2.2) r_gcrf: ', r_gcrf_r3, '[km] \n')

r_gcrf = np.matmul(Q_itrf_gcrf, r_itrf)
print('2.3) r_gcrf: ', r_gcrf, '[km]\n')

r_error = abs(r_gcrf_r3 - r_gcrf)
print('2.4) Error using R3 rotation:', r_error, '[km]')
print(np.linalg.norm(r_error))
