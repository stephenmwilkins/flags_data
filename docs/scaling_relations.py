
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flags_data.utilities import bin_centres
import flags_data.scaling_relations as sr




# --- list available relations
print(sr.list_relations())

# --- get lists of the available datasets
print(sr.list_datasets())
print(sr.list_datasets('Mstar/sSFR'))

# --- print dataset info for collections of datasets
dataset_info = sr.DatasetInfo() # all relations and models
for datasets in ['Mstar', 'Mstar/sSFR', 'Mstar/sSFR/models']:
    print('-'*50)
    dataset_info = sr.DatasetInfo(datasets = datasets)


# --- create a grid of relations and models
# dataset_info = sr.DatasetInfo()
# fig, ax = dataset_info.plot_matrix()
# # plt.show()
# fig.savefig(f'figs/matrix.pdf')


# --- create a redshift range plot of available models/observations
# relation = 'Mstar/sSFR'
# dataset_info = sr.DatasetInfo(datasets = relation)
# fig, ax = dataset_info.plot_redshift_range()
# # plt.show()
# fig.savefig(f'figs/{relation.replace("/","_")}_redshift_range.pdf')
#
#
# fig, ax = dataset_info.plot_redshift_X_range() # --- create a redshift luminosity plot of available models/observations
# # plt.show()
# fig.savefig(f'figs/{relation.replace("/","_")}_redshift_X_range.pdf')


# --- open a specific dataset and plot it
dataset = 'Mstar/sSFR/models/flares'

d = sr.read(dataset)
fig, ax = d.plot_single_z(5.0)
plt.show()

d = sr.read(dataset)
fig, ax = d.plot_z_evo()
plt.show()
