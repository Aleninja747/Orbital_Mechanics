import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import Groundtrack as g

AUS_lat = 30.3
AUS_lon = 97.7

orbits_data = os.listdir('./Props')
for path in orbits_data:
    name = path.split()[0]
    props = glob.glob('./Props/'+path+'/*.csv')
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_ylabel('Density [kg/m^3]')
    ax.set_xlabel('Height [km]')
    ax.set_title(name)
    y_twoBody = np.loadtxt('./Props/'+path+'/TwoBody_prop.csv', delimiter=',')
    r_twoBody = y_twoBody[:, 1:4]
    for file in props:
        prop_name = str(file).split('/')[-1].split('_')[0]
        if prop_name == 'MTBP':
            prop_name = 'Moon'
        elif prop_name == 'STBP':
            prop_name = 'Sun'
        if not prop_name == 'TwoBody':
            y = np.loadtxt(file, delimiter=',')
            r = y[:, 1:4]
            t = y[:, 0]
            r_diff = [0]*len(r)
            for i in range(len(r)):
                r_diff[i] = np.linalg.norm(r_twoBody[i] - r[i])
            ax.plot(t, r_diff, label=prop_name)
    ax.legend()
    plt.savefig(name, dpi=300)


burnout_data = os.listdir('./Burnout_coordinates')
for file in burnout_data:
    name = file.split('_')[0]
    data = np.loadtxt('./Burnout_coordinates/'+file, delimiter=',')
    g.plot_groundtrack(data[0, :], data[1, :])
    error = [AUS_lat-data[0, -1], AUS_lon+data[1, -1]]
    print(error)
