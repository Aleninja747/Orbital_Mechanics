import Conversions as conv
import numpy as np
import matplotlib.pyplot as plt

# Problem 5
mu = 398600.4415
a_1 = 8433
e = 0.1179
i = 19
RAAN = 310
arg_periapsis = 111
true_anomaly = 0
J2 = 0.00108248
R = 6378.1363
raan_dot = 1.991e-7

ra = a_1*(1+e)+300
rp = a_1*(1-e)
a_2 = (ra+rp)/2
print(a_2)
e_2 = (ra-rp)/(ra+rp)

r1, v1 = conv.coes2cart([a_1, e, i, RAAN, arg_periapsis, true_anomaly], mu, degree=True)
r, v2 = conv.coes2cart([a_2, e_2, i, RAAN, arg_periapsis, true_anomaly], mu, degree=True)
deltav = v2-v1
print("Problem 5")
print("Delta V =", deltav, "Km/s")
vr = np.dot(deltav, r/np.linalg.norm(r))
v_perp = np.sqrt(np.linalg.norm(deltav)**2 - vr**2)
print("Delta V = ", vr, "Vr +", v_perp, "V⊥")

# Problem 6
i_vec = np.linspace(95, 180, 10000)
e_vec = [0.1, 0.3, 0.4, 0.5, 0.9]
fig, ax = plt.subplots()
for e in e_vec:
    start = 0
    a_vec = []
    for i in i_vec:
        a = ((-3/2)*np.sqrt(mu)*J2*R**2/((1-e**2)**2*raan_dot)*np.cos(np.deg2rad(i)))**(2/7)
        if a*(1-e) > R:
            a_vec.append(a)
        else:
            start += 1
    ax.plot(i_vec[start:], a_vec, label=("e="+str(e)))

ax.legend()
ax.set_xlabel("Inclination [º]")
ax.set_ylabel("Semi-major Axis [km]")
ax.set_title("Semi-major Axis vs. Inclination")
plt.show()

# Problem 7

i = np.pi - np.arcsin(np.sqrt(4/5))
e_vec = np.linspace(0, 0.99, 10000)
a_vec = []
for e in e_vec:
    a_vec.append(((-3/2)*np.sqrt(mu)*J2*R**2/((1-e**2)**2*raan_dot)*np.cos(i))**(2/7))

fig, ax = plt.subplots()
ax.plot(a_vec, e_vec, color="#931ed6")
ax.set_xlabel("Semi-major Axis [km]")
ax.set_ylabel("Eccentricity")
ax.set_title("Semi-major Axis vs. Eccentricity")
plt.show()
