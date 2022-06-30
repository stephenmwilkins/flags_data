
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# import flare.plt as fplt

from flags_data.utilities import bin_centres
import flags_data.distribution_functions as df





df_type = 'LUV'


# --- get lists of the available for this df_type
for datasets in ['', 'models', 'models/schechter','models/binned','obs','obs/schechter','obs/binned']:
    print(datasets, df.list_datasets(f'{df_type}/{datasets}'))

# --- print dataset info for collections of datasets
dataset_info = df.DatasetInfo(datasets = df_type)
dataset_info = df.DatasetInfo(datasets = f'{df_type}/models/binned')


# --- create a redshift range plot of available models/observations
dataset_info = df.DatasetInfo(datasets = df_type)
fig, ax = dataset_info.plot_redshift_range()
plt.show()
# fig.savefig(f'figs/{df_type}DF_redshift_range.pdf')


# --- create a redshift luminosity plot of available models/observations
dataset_info = df.DatasetInfo(datasets = f'{df_type}/models/binned')
fig, ax = dataset_info.plot_redshift_log10X_range()
# plt.show()
fig.savefig(f'figs/{df_type}DF_redshift_log10X_range.pdf')



# model = 'flares'

#
#
#
# lf_binned = df.read('../data/DistributionFunctions/UVLF/models/binned/flares')
#
# lf_schechter = df.read('../data/DistributionFunctions/UVLF/models/Schechter/flares')
#
#
# # --- make a simple plot of the the luminosity function
#
# z = 5.0
#
# fig, ax = plots.simple_fig()
#
# # --- plot the binned LF
# ax.step(lf_binned.log10L[z], lf_binned.log10phi[z], where = 'mid')
#
# # --- plot the Schechter LF
# log10L = lf_binned.log10L[z] # --- use the same binned as the binned LF
# log10phi = lf_schechter.L(z).log10phi_binned(log10L) # --- determine phi (density/dex)
# ax.plot(bin_centres(log10L), log10phi)
#
# # --- add axis lavels
# ax.set_xlabel(r'$\rm \log_{10}(L/erg\ s^{-1}\ Hz^{-1}) $')
# ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ dex^{-1}) $')
#
# plt.show()
# fig.clf()
#
#
#
#
#
# # --- this will produce an automated plot for different redshifts and models
# fig, ax = plots.lf_plot(models = [lf_binned, lf_schechter], x_range = [27.5, 30.5], y_range = [-8., 0.0], redshifts = lf_schechter.redshifts)
# plt.show()
