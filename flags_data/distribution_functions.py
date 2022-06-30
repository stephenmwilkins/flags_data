

import os

# --- io modules
from astropy.table import Table
from astropy import units

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

from .utilities import log10Lnu_to_M, M_to_log10Lnu, bin_centres, simple_fig, label


this_dir = os.path.dirname(os.path.abspath(__file__))

import ads

ads.config.token = 'qm1AtsIgKukl0jqMjYaEa2LHK9am6gQka6opvce1'


def read(dataset, data_dir = f'{this_dir}/data/DistributionFunctions', interp_scheme = 'linear'):

    t = Table.read(f'{data_dir}/{dataset}.ecsv')

    if t.meta['type'] == 'Schechter':
        d = Schechter(t, scheme = interp_scheme)

    if t.meta['type'] == 'binned':
        d = Binned(t)

    if 'references' in t.meta:
        d.references = t.meta['references']
    else:
        d.references = None

    return d


# def list_datasets(datasets, data_dir = f'{this_dir}/data/DistributionFunctions'):
#     return list(map(lambda x: x.replace('.ecsv', ''), os.listdir(f'{data_dir}/{datasets}')))

def list_datasets(datasets, data_dir = f'{this_dir}/data/DistributionFunctions'):
    l = [os.path.join(dp, f.split('.')[0]) for dp, dn, fn in os.walk(os.path.expanduser(f'{data_dir}/{datasets}')) for f in fn if f.endswith('.ecsv')]
    l = [l_[len(data_dir)+1:] for l_ in l]
    return l





class DatasetInfo:

    def __init__(self, datasets = 'LUV', data_dir = f'{this_dir}/data/DistributionFunctions'):

        self.data_dir = data_dir

        self.df_type = datasets.split('/')[0]

        self.datasets = datasets
        self.models = list_datasets(datasets)

        self.short_models = [model.split('/')[-1]+'/'+model.split('/')[-2] for model in self.models]

        self.long_from_short = dict(zip(self.short_models, self.models))


    def get_info(self):

        for model_name in self.models:

            print(model_name)
            m = read(model_name, data_dir = data_dir)

            if m.df_type == 'schechter':
                print(model_name, m.redshifts)

            if m.df_type == 'binned':
                log10_limits = [np.min(m.log10X[m.redshifts[0]]), np.max(m.log10X[m.redshifts[0]])]
                print(model_name, m.redshifts, log10_limits)



    def plot_redshift_range(self, cmap = 'cmr.guppy'):

        fig = plt.figure(figsize = (5, 4.5))

        left  = 0.4
        height = 0.8
        bottom = 0.15
        width = 0.55

        ax = fig.add_axes((left, bottom, width, height))


        colors = cmr.take_cmap_colors(cmap, len(self.models))

        z_min, z_extent, labels = [], [], []

        for i, short_model_name in enumerate(sorted(self.short_models)):

            print(short_model_name)

            m = read(self.long_from_short[short_model_name], data_dir = self.data_dir)

            label = rf'$\rm \mathbf{{ {m.name} }}\ [{m.df_type}]$'
            ax.text(-0.5, i + 0.2, label, fontsize = 8, ha='right')

            if m.references:

                articles = [list(ads.SearchQuery(bibcode=bibcode))[0] for bibcode in m.references]
                refs = [f"{article.first_author.split(',')[0]}+{article.year}"  for article in articles]

                ax.text(-0.5, i - 0.2, ', '.join(refs), fontsize = 6, ha='right', color = '0.5')

            z_min.append(np.min(m.redshifts))
            z_extent.append(np.max(m.redshifts)-np.min(m.redshifts))

        ax.barh(np.arange(len(self.models)), z_extent, left = z_min, color = colors, align='center')

        ax.set_yticks([])
        ax.set_xlabel(r'$\rm z$')

        return fig, ax


    def plot_redshift_log10X_range(self, cmap = 'cmr.guppy'):

        fig, ax = simple_fig(fig_size = (4.5, 3.5))

        for model_name, color in zip(self.models, cmr.take_cmap_colors(cmap, len(self.models))):

            m = read(model_name, data_dir = self.data_dir)

            if m.df_type == 'binned':

                x, y1, y2 = [], [], []

                for z in m.redshifts:
                    x.append(z)
                    s = (~np.isnan(m.log10phi_dex[z]))&(~np.isinf(m.log10phi_dex[z]))
                    y1.append(np.min(m.log10X[z][s]))
                    y2.append(np.max(m.log10X[z][s]))

                ax.fill_between(x,y1,y2, color = color, alpha =0.2)
                ax.plot(x + x[::-1] + [x[0]], y1 + y2[::-1] + [y1[0]], lw = 1, color = color, zorder = 3)

                ax.text(x[0] + 0.1, y1[0]+0.4, rf'$\rm {m.name}$', fontsize = 10,  rotation = 90., color = color)

        ax.set_ylabel(label(m.log10x, m.log10x_unit))
        ax.set_xlabel(r'$\rm z$')

        return fig, ax












class Schechter:

    def __init__(self, t, scheme = 'linear'):

        self.df_type = 'schechter'
        self.name = t.meta['name']
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
        elif 'log10M*' in self.t.colnames:
            self.log10Mstar_star = self.t['log10M*'].data
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

        self.df_type = 'binned'
        self.name = t.meta['name']
        self.t = t


        self.log10X = {}

        # original x quantity and units
        self.x = t.meta['x']
        self.x_unit = t[self.x].unit

        if self.x in ['M', 'log10L']: self.M = {}

        if self.x == 'M':
            self.log10x = 'log10L'
            self.log10x_unit = units.Unit('dex(erg s^-1 Hz^-1)')
        else:
            self.log10x = self.x
            self.log10x_unit = self.x_unit

        print(self.name, self.log10x_unit)


        # original y quantity
        self.y = t.meta['y']
        self.y_unit = t[self.y].unit

        self.log10phi_dex = {}
        self.log10phi_mag = {}

        # --- if data is redshift, log10L, phi ... this is most useful I think
        self.redshifts = list(set(self.t['z'].data))
        self.redshifts.sort()
        # --- make sure that redshift list is monotonically increasing
        if self.redshifts[0]>self.redshifts[1]:
            self.redshifts = self.redshifts[::-1]

        # --- extract luminosities or magnitudes and number densities, making conversions where necessary.
        for z in self.redshifts:

            s = self.t['z'] == z

            if self.x in ['log10L','log10Mstar','log10SFR']:
                self.log10X[z] = self.t[self.x][s].data
                if self.x in ['M', 'log10L']: self.M[z] = log10Lnu_to_M(self.log10X[z])
            elif self.x == 'M':
                self.M[z] = self.t['M'][s].data
                self.log10X[z] = M_to_log10Lnu(self.M[z])
            else:
                print('WARNING [x]')

            if self.y == 'phi':
                phi = self.t['phi'][s].data
                log10phi = np.log10(phi)
            elif self.y == 'log10phi':
                log10phi = self.t['log10phi'][s].data
                phi = 10**log10phi
            else:
                print('WARNING [phi]')

            if self.y_unit == units.Unit('1 / (dex Mpc3)') or self.y_unit == units.Unit('dex(1 / (dex Mpc3))'):
                self.log10phi_dex[z] = log10phi
                self.log10phi_mag[z] = log10phi + np.log10(0.4)
            elif self.y_unit == units.Unit('1 / (mag Mpc3)') or self.y_unit == units.Unit('dex(1 / (mag Mpc3))'):
                self.log10phi_mag[z] = log10phi
                self.log10phi_dex[z] = log10phi - np.log10(0.4)
            else:
                print('WARNING [unit]')








class plots:

    def lf(models = None, redshifts = None, cmap = 'cmr.guppy', lum_type = 'Lnu', x_range = [27., 30.], y_range = [-7., 1.]):

        fig, ax = simple_fig()

        colors = cmr.take_cmap_colors(cmap, len(redshifts))

        for m in models:

            for z, color in zip(redshifts, colors):

                if m.df_type == 'Binned':
                    ax.step(m.log10L[z], m.log10phi[z], where = 'mid', color = color, label = rf'$\rm z={z} $')

                if m.df_type == 'Schechter':
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
