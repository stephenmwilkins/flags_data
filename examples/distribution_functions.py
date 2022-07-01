
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


datasets = 'LUV/models/binned'
z = 10.0

di = df.DatasetInfo(datasets = datasets)
di.get_info()
dataset_z = di.get_datasets_at_z(z)


# --- plot a list of (dataset, z); in this case all at the same redshift.
fig, ax = df.plots.df(dataset_z)
fig.savefig(f"figs/{datasets.replace('/','-')}_{z}.pdf")
