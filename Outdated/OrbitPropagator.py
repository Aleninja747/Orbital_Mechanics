import numpy as np
import planetary_data as pd
from scipy.integrate import ode


class OrbitPropagator:
    def __init__(self, r, v, tspan, dt, cb=pd.earth, j2=False):
        self.r0 = r
        self.v0 = v
        self.cb = cb
        self.tspan = tspan
        self.dt = dt
        self.j2 = j2

        self.y0 = self.r0.tolist() + self.v0.tolist()

        self.n_steps = int(np.ceil(self.tspan / self.dt))

        # initialize arrays
        self.y = np.zeros((self.n_steps, 6))
        self.ts = np.zeros((self.n_steps, 1))
        self.rs = np.zeros((self.n_steps, 3))
        self.vs = np.zeros((self.n_steps, 3))

        # initial conditions
        self.ts[0] = 0
        self.y[0] = self.y0
        self.step = 1

        # initiate solver
        self.solver = ode(self.dif_eq)
        self.solver.set_integrator('dopri5', rtol=1e-4, atol=1e-4)
        self.solver.set_initial_value(self.y0, 0)

    def propagate_orbit(self):
        # propagate orbit
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.step] = self.solver.t
            self.y[self.step] = self.solver.y
            self.step += 1

        self.ts = self.ts[:self.step]
        self.rs = self.y[:, :3]
        self.vs = self.y[:, 3:]

    def dif_eq(self, t, y):
        # unpack state
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])

        # norm of the radius vector
        norm_r = np.linalg.norm(r)

        # two bod acceleration
        a = -r * self.cb['mu'] / norm_r ** 3

        # J2 Perturbation
        if self.j2:
            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)
            a += 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius'] ** 2 / norm_r ** 4 * np.array([tx, ty, tz])

        return [vx, vy, vz, a[0], a[1], a[2]]
