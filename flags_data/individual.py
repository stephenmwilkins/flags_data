

import os

# --- io modules
from astropy.table import Table
import astropy.units as units

import numpy as np
from numpy.random import randn

# --- plotting modules/functions
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import cmasher as cmr # provides wider range of cmaps, see https://cmasher.readthedocs.io

from .utilities import log10Lnu_to_M, M_to_log10Lnu, bin_centres, simple_fig, label

this_dir = os.path.dirname(os.path.abspath(__file__))

data_dir = f'{this_dir}/data/Individual'


# --- default units

default_units = {}
default_units['Mstar'] = units.Unit('dex(solMass)')
default_units['SFR'] = units.Unit('dex(solMass yr^-1)')
default_units['sSFR'] = units.Unit('dex(Gyr^-1)')
default_units['Zstar'] = None
default_units['Zstar_young'] = None
default_units['Zgas_young'] = None
default_units['Zgas'] = None
default_units['age'] = units.Unit('Myr')
default_units['Rstar'] = units.Unit('kpc')
default_units['LUV'] = units.Unit('dex(erg s^-1 Hz^-1)')
default_units['beta'] = None
default_units['RUV'] = None


def list_datasets(dataset_type = '', data_dir = data_dir):
    l = [os.path.join(dp, f.split('.')[0]) for dp, dn, fn in os.walk(os.path.expanduser(f'{data_dir}/{dataset_type}')) for f in fn if f.endswith('.ecsv')]
    l = [l_[len(data_dir)+1:] for l_ in l]
    return l


class Dataset:

    def __init__(self, dataset, data_dir = data_dir, verbose = False):

        t = Table.read(f'{data_dir}/{dataset}.ecsv')
        self.t = t

        self.name = t.meta['name']

        if verbose: print(self.name, '-'*20)

        if 'references' in t.meta:
            self.references = t.meta['references']
        else:
            self.references = None

        self.quantities = t.meta['quantities']
        self.observatory = t.meta['observatory']

        self.data = {}


        for c in t.colnames:
            self.data[c] = t[c].data

        # if verbose: print(self.quantities)

        for q in self.quantities:

            self.data[q] = self.t[q].data

            # --- add errors if not available
            for k in ['low','upp']:
                if q+'_err_'+k not in t.colnames:
                    self.data[q+'_err_'+k] = np.zeros(len(self.data[q]))
                else:
                    self.data[q+'_err_'+k] = self.t[q+'_err_'+k].data

            # if verbose: print(list(self.data.keys()))

            if len(q)>5:
                if q[:5] == 'log10':
                    self.data[q[5:]] = 10**self.data[q]
                    self.data[q[5:]+'_err_low'] = 10**(self.data[q]) - 10**(self.data[q]-self.data[q+'_err_low'])
                    self.data[q[5:]+'_err_upp'] = 10**(self.data[q]+self.data[q+'_err_low']) - 10**(self.data[q])

            else:
                if q.split('_')[0] not in ['M_UV', 'A_V', 'beta']:
                    self.data['log10'+q] = np.log10(self.data[q])
                    self.data['log10'+q+'_err_low'] = np.log10(self.data[q]) - np.log10(self.data[q]-self.data[q+'_err_low'])
                    self.data['log10'+q+'_err_upp'] = np.log10(self.data[q]+self.data[q+'_err_upp']) - np.log10(self.data[q])

        if 'sSFR' not in self.data.keys():
            if 'SFR' in self.data.keys() and 'Mstar' in self.data.keys():
                self.data['log10sSFR'] = self.data['log10SFR'] - self.data['log10Mstar'] + 9.0
                self.data['log10sSFR_err_low'] = np.sqrt(self.data['log10SFR_err_low']**2 + self.data['log10Mstar_err_low']**2)


                self.data['log10sSFR_err_upp'] = np.sqrt(self.data['log10SFR_err_upp']**2 + self.data['log10Mstar_err_upp']**2)

                # print('low errors:', self.data['log10SFR_err_low'], self.data['log10Mstar_err_low'], self.data['log10sSFR_err_low'])
                # print('high errors:', self.data['log10SFR_err_upp'], self.data['log10Mstar_err_upp'], self.data['log10sSFR_err_upp'])

        # if verbose:
        #     for k,v in self.data.items():
        #         print(k)

    def get_values(self, x, y, z, z_tolerance = 1.0):
        print(x, y, z)

        if (x in self.data.keys()) and (y in self.data.keys()):
            s = np.fabs(self.data['z']-z)<z_tolerance
            if np.sum(s)>0:
                return self.data['id'][s], self.data[x][s], self.data[y][s], self.data['z'][s]


    def get_values_errs(self, x, y, z, z_tolerance = 0.5):

        if (x in self.data.keys()) and (y in self.data.keys()):
            s = np.fabs(self.data['z']-z)<z_tolerance
            if np.sum(s)>0:

                X = self.data[x][s]
                X_err = [self.data[x+'_err_low'][s], self.data[x+'_err_upp'][s]]

                Y = self.data[y][s]
                Y_err = [self.data[y+'_err_low'][s], self.data[y+'_err_upp'][s]]

                print(y+'_err_upp', self.data[y+'_err_upp'])
                print(y+'_err_low', self.data[y+'_err_low'])
                print('Y_err', Y_err)

                return self.data['id'][s], X, X_err, Y, Y_err, self.data['z'][s]



class Datasets:

    def __init__(self, datasets = '', data_dir = data_dir):

        self.dataset_list = list_datasets(datasets)

        self.d = {ds: Dataset(ds, data_dir = data_dir) for ds in self.dataset_list}

class Datasets_:

    def __init__(self, dataset_list = '', data_dir = data_dir):

        self.dataset_list = dataset_list

        self.d = {ds: Dataset(ds, data_dir = data_dir) for ds in self.dataset_list}
