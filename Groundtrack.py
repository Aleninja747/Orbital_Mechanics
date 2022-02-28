import numpy as np
import matplotlib.pyplot as plt
import OrbitPropagatorV2 as op
import Conversions as conv
import cartopy.crs as ccrs
import Perturbations as perts
from TimeClassV2 import Time
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def generate_groundtracks(y0, lst_green0, omega_earth, tspan, dt=100):
    f = 3.353e-3
    perts_dic = perts.perts_dic
    perts_dic['J2'] = True
    init_time = Time('UTC', day=2451545.0)
    ys = np.array(op.propagate_orbit(y0, init_time, tspan, dt, perts_dic=perts_dic))
    ts = ys[:, 0]
    rs = ys[:, 1:4]
    latitudes = []
    longitudes = []
    for i in range(len(ts)):
        lst_green = lst_green0 + omega_earth * ts[i]
        r_ecf = np.matmul(conv.euler_r3(lst_green), rs[i])
        gc_lat = np.arcsin(r_ecf[2]/np.linalg.norm(r_ecf))
        latitudes.append(np.rad2deg(np.arctan2(np.tan(gc_lat), ((1-f)**2))))
        longitudes.append(np.rad2deg(np.arctan2(r_ecf[1], r_ecf[0])))

    return latitudes, longitudes


def calc_sma_rgt(m, n, mu=398600.4415, omega_earth=2*np.pi/86164):
    return np.cbrt(mu*(m/(n*omega_earth))**2)


def plot_groundtrack(latitudes, longitudes):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.set_global()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    plt.plot(longitudes, latitudes, '.r', markersize=3)
    plt.plot(longitudes[0], latitudes[0], '.g', markersize=5)
    plt.plot(longitudes[-1], latitudes[-1], '.b', markersize=5)

    #  Set the correct aspect ratio for the plot

    #  Show the plot
    plt.show()
