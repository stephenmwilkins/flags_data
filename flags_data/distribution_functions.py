

import os

# --- io modules
from astropy.table import Table

import numpy as np
from numpy.random import randn

# --- stats modules
from scipy.stats import linregress
import scipy.integrate as cp
import scipy.interpolate as cpi
import scipy.special as cps
from scipy.stats import rv_histogram
from mpmath import gammainc

# --- plotting modules/functions
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import cmasher as cmr # provides wider range of cmaps, see https://cmasher.readthedocs.io

from .utilities import log10Lnu_to_M, M_to_log10Lnu, bin_centres, simple_fig


this_dir = os.path.dirname(os.path.abspath(__file__))



def read(model_name, df_type = 'UVLF', data_dir = f'{this_dir}/data/DistributionFunctions', interp_scheme = 'linear'):

    t = Table.read(f'{data_dir}/{df_type}/{model_name}.ecsv')

    if t.meta['type'] == 'Schechter':
        return Schechter(t, scheme = interp_scheme)

    if t.meta['type'] == 'binned':
        return Binned(t)


def list_data(path):
    return list(map(lambda x: x.replace('.ecsv', ''), os.listdir(path)))

class ListDatasets:

    def __init__(self, df_type = 'UVLF', data_dir = f'{this_dir}/data/DistributionFunctions'):
        self.df_type = df_type
        self.data_dir = data_dir


    def models_schechter(self):
        print(f'{self.data_dir}/{self.df_type}/models/schechter/')
        return list_data(f'{self.data_dir}/{self.df_type}/models/schechter/')

    def models_binned(self):
        return list_data(f'{self.data_dir}/{self.df_type}/models/binned/')

    def models(self):
        return [f'binned/{x}' for x in self.models_binned()] + [f'schechter/{x}' for x in self.models_schechter()]

    def obs_schechter(self):
        return list_data(f'{self.data_dir}/{self.df_type}/obs/schechter/')

    def obs_binned(self):
        return list_data(f'{self.data_dir}/{self.df_type}/obs/binned/')

    def obs(self):
        return [f'binned/{x}' for x in self.obs_binned()] + [f'schechter/{x}' for x in self.obs_schechter()]

    def all(self):
        return [f'models/{x}' for x in self.models()] + [f'obs/{x}' for x in self.obs()]

    def list(self, datasets, full = False):
        if full:
            return [f'{datasets}/{x}' for x in getattr(self, datasets.replace('/', '_'))()]
        else:
            return getattr(self, datasets.replace('/', '_'))()




class DatasetInfo:

    def __init__(self, df_type = 'UVLF', data_dir = f'{this_dir}/data/DistributionFunctions', datasets = 'all'):

        self.data_dir = data_dir
        self.df_type = df_type

        if datasets != 'all':
            self.models = ListDatasets(df_type = df_type, data_dir = data_dir).list(datasets, full = True)
        else:
            self.models = ListDatasets(df_type = df_type, data_dir = data_dir).list(datasets)

        for model_name in self.models:

            m = read(model_name, data_dir = data_dir, df_type = df_type)

            if m.lf_type == 'schechter':
                print(model_name, m.redshifts)

            if m.lf_type == 'binned':

                log10_limits = [np.min(m.log10L[m.redshifts[0]]), np.max(m.log10L[m.redshifts[0]])]

                print(model_name, m.redshifts, log10_limits)


    def plot_redshift_luminosity_range(self, cmap = 'cmr.guppy'):

        fig, ax = simple_fig(fig_size = (4.5, 3.5))

        for model_name, color in zip(self.models, cmr.take_cmap_colors(cmap, len(self.models))):

            m = read(model_name, data_dir = self.data_dir, df_type = self.df_type)

            if m.lf_type == 'binned':

                x, y1, y2 = [], [], []

                for z in m.redshifts:
                    x.append(z)
                    s = (m.phi[z]>0)&(~np.isnan(m.phi[z]))
                    y1.append(np.min(m.log10L[z][s]))
                    y2.append(np.max(m.log10L[z][s]))

                print(model_name, x,y1,y2)

                ax.fill_between(x,y1,y2, color = color, alpha =0.2)
                ax.plot(x + x[::-1] + [x[0]], y1 + y2[::-1] + [y1[0]], lw = 1, color = color, zorder = 3)

                ax.text(x[0] + 0.1, y1[0]+0.4, rf'$\rm {m.name}$', fontsize = 10,  rotation = 90., color = color)

        ax.set_ylabel(r'$\rm \log_{10}(L/erg\ s^{-1}\ Hz^{-1}) $')
        ax.set_xlabel(r'$\rm z$')

        return fig, ax





# def get_models_at_redshift(data_dir, z = None, model_types = ['schechter', 'binned'], verbose = False):
#
#     models = {}
#
#     if 'schechter' in model_types:
#         for model in list_models_schechter():
#
#             model_name = f'models/schechter/{model}'
#             m = read(model_name)
#
#             for redshift in m.redshifts:
#
#                 redshift = int(redshift)
#
#                 if redshift in models.keys():
#                     models[redshift].append(model_name)
#                 else:
#                     models[redshift] = [model_name]
#
#     if 'binned' in model_types:
#         for model in list_models_binned():
#
#             model_name = f'models/binned/{model}'
#             m = read(model_name)
#
#             for redshift in m.redshifts:
#
#                 redshift = int(redshift)
#
#                 if redshift in models.keys():
#                     models[redshift].append(model_name)
#                 else:
#                     models[redshift] = [model_name]
#
#     if verbose:
#         for k, v in models.items():
#             print(k, v)
#
#     if z:
#         return(models[z])
#     else:
#         return(models)









class Schechter:

    def __init__(self, t, scheme = 'linear'):

        self.lf_type = 'schechter'
        self.name = t.meta['name']
        if 'ref' in t.meta: self.ref = t.meta['ref']
        self.t = t
        self.redshifts = self.t['redshift'].data

        if 'alpha' in self.t.colnames:
            self.alpha = self.t['alpha'].data
        else:
            print('WARNING: No faint-end slope set')

        if 'phi*' in self.t.colnames:
            self.phi_star = self.t['log10phi*'].data
            self.log10phi_star = np.log10(self.phi_star)
        elif 'log10phi*' in self.t.colnames:
            self.log10phi_star = self.t['log10phi*'].data
            self.phi_star = 10**(self.log10phi_star)
        else:
            print('WARNING: No characteristic density set')

        if 'M*' in self.t.colnames:
            self.M_star = self.t['M*'].data
            self.log10L_star = M_to_log10Lnu(self.M_star)
        elif 'log10L*' in self.t.colnames:
            self.log10L_star = self.t['log10L*'].data
            self.M_star = log10Lnu_to_M(self.log10L_star)
        else:
            print('WARNING: No characteristic luminosity or magnitude set')


        # --- create dictionary of parameters for each redshift
        self.p = {}

        for i,z in enumerate(self.redshifts):
            self.p[z] = {'M*': self.M_star[i], 'log10L*': self.log10L_star[i], 'log10phi*': self.log10phi_star[i], 'alpha': self.alpha[i]}

    def L(self, z):
        return self.luminosity(self.p[z])


    class luminosity:

        """ Get luminosity function, phi, etc. in terms of luminosity (not magnitude) """

        def __init__(self, p):

            self.alpha = p['alpha']
            self.phistar = 10**p['log10phi*']
            self.Lstar = 10**p['log10L*']

        def phif(self, x):
            return x ** (self.alpha) * np.exp(-x)

        def phi(self, L):

            return self.phistar * self.phif(L/self.Lstar)

        def density(self, L):

            """ get the density down to some limit """

            return self.phistar * self.Lstar * np.float(gammainc(self.alpha + 2, L/self.Lstar))


        def phi_binned(self, log10L):

            """ integrate the LF between the bin edges to get the number density of galaxies in the bin """

            y = np.zeros(len(log10L)-1)

            for i, (a,b) in enumerate(zip(log10L[:-1], log10L[1:])):
                y[i] = self.phistar * cp.quad(self.phif, 10**a/self.Lstar, 10**b/self.Lstar)[0]

            return y/(log10L[1]-log10L[0]) # divide by the bin-width to yield density/dex

        def log10phi_binned(self, log10L):

            """ return log10 of the above """

            return np.log10(self.phi_binned(log10L))









class Binned:

    def __init__(self, t):

        self.lf_type = 'binned'
        self.name = t.meta['name']
        if 'ref' in t.meta: self.ref = t.meta['ref']
        self.t = t

        self.M = {}
        self.log10L = {}

        self.phi = {}
        self.log10phi = {}
        self.log10phi_dex = {}
        self.log10phi_mag = {}

        # --- if data is redshift, log10L, phi ... this is most useful I think
        self.redshifts = list(set(self.t['redshift'].data))
        self.redshifts.sort()
        # --- make sure that redshift list is monotonically increasing
        if self.redshifts[0]>self.redshifts[1]:
            self.redshifts = self.redshifts[::-1]

        # --- extract luminosities or magnitudes and number densities, making conversions where necessary.
        for z in self.redshifts:

            if 'log10L' in self.t.colnames:
                self.log10L[z] = self.t['log10L'][self.t['redshift']==z].data

                self.M[z] = log10Lnu_to_M(self.log10L[z])

                if 'log10phi' in self.t.colnames:
                    self.log10phi_dex[z] = self.t['log10phi'][self.t['redshift']==z].data
                elif 'phi' in self.t.colnames:
                    self.log10phi_dex[z] = np.log10(self.t['phi'][self.t['redshift']==z].data)
                else:
                    print('WARNING: one column should be log10phi or phi')

                self.log10phi_mag[z] = self.log10phi_dex[z] + np.log10(0.4)

            elif 'M' in self.t.colnames:
                self.M[z] = self.t['M'][self.t['redshift']==z].data
                self.log10L[z] = M_to_log10Lnu(self.M[z])
                # print(z, len(self.log10L[z]))
                if 'log10phi' in self.t.colnames:
                    self.log10phi_mag[z] = self.t['log10phi'][self.t['redshift']==z].data
                elif 'phi' in self.t.colnames:
                    self.log10phi_mag[z] = np.log10(self.t['phi'][self.t['redshift']==z].data)
                else:
                    print('WARNING: one column should be log10phi or phi')

                self.log10phi_dex[z] = self.log10phi_mag[z] - np.log10(0.4)

            else:
                print('no luminosity/magnitude column found, use log10L or M')

            self.log10phi[z] = self.log10phi_dex[z]
            self.phi[z] = 10**self.log10phi[z]











def lf_plot(models = None, redshifts = None, cmap = 'cmr.guppy', lum_type = 'Lnu', x_range = [27., 30.], y_range = [-7., 1.]):

    fig, ax = simple_fig()

    colors = cmr.take_cmap_colors(cmap, len(redshifts))

    for m in models:

        for z, color in zip(redshifts, colors):

            if m.lf_type == 'Binned':
                ax.step(m.log10L[z], m.log10phi[z], where = 'mid', color = color, label = rf'$\rm z={z} $')

            if m.lf_type == 'Schechter':
                log10L = np.arange(*x_range, 0.1)
                ax.plot(bin_centres(log10L), m.L(z).log10phi_binned(log10L), color = color)


    ax.legend(fontsize = 8)
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)

    if lum_type == 'Lnu':
        ax.set_xlabel(r'$\rm \log_{10}(L/erg\ s^{-1}\ Hz^{-1}) $')
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ dex^{-1}) $')

    if lum_type == 'M':
        ax.set_xlabel(r'$\rm \log_{10}(M) $')
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ mag^{-1}) $')

    return fig, ax
