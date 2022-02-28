import numpy as np

G_meters = 6.67408e-11

G = G_meters * 10 ** -9

earth = {
    'name': 'Earth',
    'mass': 5.972e24,
    'mu': 5.972e24 * G,
    'radius': 6378.0,
    'J2': 1.082635854e-3,
}