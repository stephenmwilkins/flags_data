
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import flare.plt as fplt

from flags_data.utilities import bin_centres
import flags_data.distribution_functions as df



def make_tables(df_type):


    # --- get lists of the available for this df_type
    for datasets in ['', 'models', 'models/schechter','models/binned','obs','obs/schechter','obs/binned']:
        print(datasets, df.list_datasets(f'{df_type}/{datasets}'))

    # --- print dataset info for collections of datasets
    dataset_info = df.DatasetInfo(datasets = df_type)
    dataset_info = df.DatasetInfo(datasets = f'{df_type}/models/binned')


def make_range_plots(df_type):

    # --- create a redshift range plot of available models/observations
    dataset_info = df.DatasetInfo(datasets = df_type)
    fig, ax = dataset_info.plot_redshift_range()
    fig.savefig(f'figs/{df_type}DF_redshift_range.png')

    # --- create a redshift luminosity plot of available models/observations
    dataset_info = df.DatasetInfo(datasets = f'{df_type}/models/binned')
    fig, ax = dataset_info.plot_redshift_log10X_range()
    fig.savefig(f'figs/{df_type}DF_redshift_log10X_range.pdf')
    fig.savefig(f'figs/{df_type}DF_redshift_log10X_range.png')



if __name__ == "__main__":

    for df_type in ['LUV','Mstar','SFR']: #,'SFR'
        make_range_plots(df_type)
