
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flags_data.utilities import bin_centres
import flags_data.scaling_relations as sr
import flags_data.individual as ind


relation = 'Mstar/sSFR'


redshift = 14.


# --- open a specific dataset and plot it
d = sr.read(f'{relation}/models/flares')
fig, ax = d.plot_single_z(redshift)


# --- add observations

x,y = 'log10Mstar', 'log10SFR'

obs = ind.Datasets(datasets = 'obs')

for ds, d in obs.d.items():
    print(ds, '-'*30)

    # v = d.get_values(x,y,redshift)
    # if v: id, X, Y, z = v
    # ax.scatter(X,Y)

    v = d.get_values_errs(x,y,redshift, z_tolerance = 4)
    if v:
        print(v)
        id, X, X_err, Y, Y_err, z = v
        ax.errorbar(X,Y, xerr = X_err, yerr = Y_err, fmt='o', elinewidth=1, ms=5, zorder = 2)


# relation = 'Mstar/sSFR'
#
#
# redshift = 7.
#
#
# # --- open a specific dataset and plot it
# d = sr.read(f'{relation}/models/flares')
# fig, ax = d.plot_single_z(redshift)
#
#
# # --- add observations
#
# x,y = 'log10Mstar', 'log10sSFR'
#
# obs = ind.Datasets(datasets = 'obs')
#
# for ds, d in obs.d.items():
#
#     v = d.get_values(x,y,redshift)
#     if v:
#         id, X, Y, z = v
#         ax.scatter(X,Y)
#
#     # v = d.get_values_errs(x,y,redshift, z_tolerance = 4)
#     # if v:
#     #     id, X, X_err, Y, Y_err, z = v
#     #     ax.errorbar(X,Y, xerr = X_err, yerr = Y_err)
#

plt.show()
