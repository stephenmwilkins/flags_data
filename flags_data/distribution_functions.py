

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

from .utilities import log10Lnu_to_M, M_to_log10Lnu, bin_centres, simple_fig, label, label_


this_dir = os.path.dirname(os.path.abspath(__file__))



data_dir = f'{this_dir}/data/DistributionFunctions'

x_ranges = {}
x_ranges['LUV'] = [27.51, 30.]
x_ranges['SFR'] = [-0.5, 3.]
x_ranges['Mstar'] = [7.5, 11.5]



def read(dataset, data_dir = data_dir, interp_scheme = 'linear'):

    print(f'Reading {dataset}')

    t = Table.read(f'{data_dir}/{dataset}.ecsv')

    if t.meta['type'] == 'Schechter':
        d = Schechter(t, scheme = interp_scheme)

    if t.meta['type'] == 'binned':
        d = Binned(t)

    if 'references' in t.meta:
        d.references = t.meta['references']
    else:
        d.references = None

    d.id = dataset.split('/')[-1]

    if 'intrinsic' in t.meta:
        d.intrinsic = t.meta['intrinsic']

    return d


# def list_datasets(datasets, data_dir = f'{this_dir}/data/DistributionFunctions'):
#     return list(map(lambda x: x.replace('.ecsv', ''), os.listdir(f'{data_dir}/{datasets}')))

def list_datasets(datasets, data_dir = data_dir, remove_repeats = False):
    l = [os.path.join(dp, f.split('.')[0]) for dp, dn, fn in os.walk(os.path.expanduser(f'{data_dir}/{datasets}')) for f in fn if f.endswith('.ecsv')]
    dataset_list = [l_[len(data_dir)+1:] for l_ in l]

    if remove_repeats:
        dataset_list_ = [x.split('/') for x in dataset_list]
        dataset_list_ = sorted(dataset_list_, key = lambda x: x[-1])

        studies = np.array([x[-1] for x in dataset_list_])

        # --- get rid of the Schechter parameters where we have the binned version
        if remove_repeats:
            for x in set(studies):
                if np.sum(studies==x)>1:
                    x_ = [datasets.split('/')[0], 'models', 'schechter', x] # eughh
                    dataset_list_.remove(x_)

        dataset_list = ['/'.join(x) for x in dataset_list_]



    return dataset_list





class DatasetInfo:

    def __init__(self, datasets = 'LUV', data_dir = data_dir, remove_repeats = True):

        self.data_dir = data_dir

        if type(datasets) == str:
            self.df_type = datasets.split('/')[0]
            self.datasets_ = datasets
            self.dataset_list = list_datasets(datasets)
        else:
            self.dataset_list = datasets


        self.datasets = [x.split('/') for x in self.dataset_list]
        self.datasets = sorted(self.datasets, key = lambda x: x[-1])

        self.studies = np.array([x[-1] for x in self.datasets])

        # --- get rid of the Schechter parameters where we have the binned version
        if remove_repeats:
            for x in set(self.studies):
                if np.sum(self.studies==x)>1:
                    x_ = [datasets.split('/')[0], 'models', 'schechter', x] # eughh
                    self.datasets.remove(x_)

        self.dataset_list = ['/'.join(x) for x in self.datasets]
        self.n = len(self.dataset_list)

        # --- read all datasets
        self.dataset = {}
        self.redshifts = {}
        self.log10X_range = {}
        self.names = {}
        for dataset_name in self.dataset_list:

            m = read(dataset_name, data_dir = data_dir)

            self.dataset[dataset_name] = m

            # --- get dataset name
            self.names[dataset_name] = m.name

            # --- get redshifts of each models
            self.redshifts[dataset_name] = m.redshifts

            # --- get log10X range of binned models
            if dataset_name.split('/')[-1] == 'binned':
                self.log10X_range[dataset_name] = [np.min(m.log10X[m.redshifts[0]]), np.max(m.log10X[m.redshifts[0]])]


    def get_info(self):
        for dataset_name in self.dataset_list:
            if dataset_name.split('/')[-1] == 'binned':
                print(dataset_name, self.redshifts[dataset_name], self.log10X_range[dataset_name])
            if dataset_name.split('/')[-1] == 'schechter':
                print(dataset_name, self.redshifts[dataset_name])


    # def get_datasets_at_z(self, z, z_tolerance = 0.2, remove_repeats = True):
    #     datasets_z = []
    #     datasets_ = []
    #     z_ = []
    #
    #     for dataset_name in self.dataset_list:
    #         for z_ in self.redshifts[dataset_name]:
    #             if np.fabs(z-z_)<z_tolerance:
    #                 datasets_z.append((dataset_name, z_))


    def get_datasets_at_z(self, z, z_tolerance = 0.2, remove_repeats = True):
        datasets_z = []

        if remove_repeats:
            redshift = {}
            for dataset_name in self.dataset_list:
                for z_ in self.redshifts[dataset_name]:
                    if np.fabs(z-z_)<z_tolerance:
                            if dataset_name not in redshift.keys():
                                redshift[dataset_name] = z_
                            else:
                                if np.fabs(z-z_)<np.fabs(z-redshift[dataset_name]):
                                    redshift[dataset_name] = z_
            datasets_z = [(k,v) for k,v in redshift.items()]

        else:
            for dataset_name in self.dataset_list:
                for z_ in self.redshifts[dataset_name]:
                    if np.fabs(z-z_)<z_tolerance:
                        datasets_z.append((dataset_name, z_))


        return(datasets_z)

    def plot_redshift_range(self, cmap = 'cmr.guppy'):


        xtotal_ = 7.  # in "

        bottom_ = 0.35  # in "
        height_ = 0.15 * self.n  # in "
        top_ = 0.1  # in "
        ytotal_ = bottom_ + height_ + top_  # in "

        left  = 0.4
        width = 0.55

        bottom = bottom_/ytotal_
        height = height_/ytotal_

        fig = plt.figure(figsize = (xtotal_, ytotal_))

        ax = fig.add_axes((left, bottom, width, height))

        colors = cmr.take_cmap_colors(cmap, self.n)

        for i_, ((df, modobs, df_type, study), color) in enumerate(zip(self.datasets, colors)):

            i = self.n-i_-1

            dataset_name = f'{df}/{modobs}/{df_type}/{study}'
            m = self.dataset[dataset_name]

            nm = m.name.replace(' ','\ ')

            label = rf'$\rm \mathbf{{ {nm} }}\ [{m.df_type}]$'
            print(i, df_type, study, nm)
            ax.text(3.3, i, label, fontsize = 8, ha='right', va = 'center')

            ax.barh(i, np.max(m.redshifts)-np.min(m.redshifts), left = np.min(m.redshifts), color = color, align='center', alpha = 0.5, height = 1)
            for z in m.redshifts:
                ax.plot([z,z],[i-0.5,i+0.5], color = color, lw=1, solid_capstyle='butt')

        ax.set_xlim([3.5, 15.5])
        ax.set_ylim([-0.5, self.n - 0.5])
        ax.set_yticks([])
        ax.set_xlabel(r'$\rm z$')

        return fig, ax

    def plot_redshift_log10X_range(self, cmap = 'cmr.guppy'):

        n_datasets = len(self.dataset_list)

        fig, ax = simple_fig(fig_size = (4.5, 3.5))

        colors = cmr.take_cmap_colors(cmap, self.n)

        for i_, ((df, modobs, df_type, study), color) in enumerate(zip(self.datasets, colors)):

            i = self.n-i_-1

            dataset_name = f'{df}/{modobs}/{df_type}/{study}'
            m = self.dataset[dataset_name]

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

        ax.set_xlim([3.5, 15.5])
        # ax.set_ylim(x_ranges[self.df_type])
        # ax.set_ylabel(label(m.log10x, m.log10x_unit))
        ax.set_ylabel(label_(self.df_type))
        ax.set_xlabel(r'$\rm z$')

        return fig, ax


    def plot_dfs(self, redshifts = np.arange(5, 17, 1), cmap = 'cmr.guppy', x_range = None, y_range = [-6.99, -0.51]):

        if not x_range:
            x_range = x_ranges[self.df_type]


        fig = plt.figure(figsize = (7,5))

        left = 0.1
        right = 0.75
        bottom = 0.1
        top = 0.95

        gs = fig.add_gridspec(4, 3, left = left, bottom = bottom, right = right, top = top, hspace=0, wspace=0)
        axes = gs.subplots(sharex=True, sharey=True)

        colors = dict(zip(self.dataset_list, cmr.take_cmap_colors(cmap, self.n)))

        lss = dict(zip(self.dataset_list, ['-','--','-.',':']*10))

        print(self.names)

        # --- create legend
        lax = fig.add_axes([right, bottom, 0.3, top-bottom])
        lax.axis('off')
        handles = [Line2D([0], [0], color = colors[ds], lw=2, ls=lss[ds], label=self.names[ds]) for ds in self.dataset_list]
        lax.legend(handles=handles, loc='center left', title = self.datasets_, fontsize = 8)

        for ax, z in zip(axes.flatten(), redshifts):

            ax.label_outer()

            ax.text(0.05, 0.9, rf'$\rm z={z:.0f}$', color = '0.3', transform = ax.transAxes)

            dataset_z = self.get_datasets_at_z(z, z_tolerance = 0.51) # --- get list of datasets at this redshift



            for dataset_name, z in dataset_z:


                _, obsmod, df_type, study = dataset_name.split('/')

                dataset = self.dataset[dataset_name]

                color = colors[dataset_name]
                ls  = lss[dataset_name]

                if df_type == 'binned' and obsmod == 'models':
                    ax.plot(dataset.log10X[z], dataset.log10phi_dex[z], color = color, ls = ls)

                if df_type == 'binned' and obsmod == 'obs':
                    ax.scatter(dataset.log10X[z], dataset.log10phi[z], color = 'k', label = rf'$\rm {dataset.name} $', s = 3)
                    ax.errorbar(dataset.log10X[z], dataset.log10phi[z], xerr = dataset.log10X_binw[z]/2., yerr = dataset.log10phi_err[z], fmt='o', elinewidth=1, ms=4, c='k', mec='k', zorder = 2)

                if df_type == 'schechter':
                    log10L = np.arange(*x_range, 0.01)
                    ax.plot(bin_centres(log10L), dataset.L(z).log10phi_binned(log10L), color = color, ls = ls)



            ax.set_xlim(x_range)
            ax.set_ylim(y_range)
            ax.set_xticks(np.arange(*np.round(x_range, 0), 1))

            # ax.legend(fontsize = 8)

        axes[2,1].set_xlabel(label_(self.df_type), fontsize = 10)
        axes[1,0].set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ dex^{-1}) $', fontsize = 10)

        return fig, axes




class Schechter:

    def __init__(self, t, scheme = 'linear'):

        self.df_type = 'schechter'
        self.name = t.meta['name']

        nm = self.name.replace(' ', '\ ')
        self.label = rf'$\rm {nm}$'

        if 'shortname' in t.meta.keys():
            self.slabel = rf'$\rm {t.meta["shortname"]}$'
        else:
            self.slabel = self.label


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

        nm = self.name.replace(' ', r'\ ')
        self.label = rf'$\rm {nm}$'

        if 'shortname' in t.meta.keys():
            print(self.name)
            self.slabel = rf'$\rm {t.meta["shortname"]}$'
        else:
            self.slabel = self.label

        # original x quantity and units
        self.x = t.meta['x']
        self.x_unit = t[self.x].unit

        if self.x == 'M':
            self.log10x = 'log10L'
            self.log10x_unit = units.Unit('dex(erg s^-1 Hz^-1)')
        else:
            self.log10x = self.x
            self.log10x_unit = self.x_unit

        # original y quantity
        self.y = t.meta['y']
        self.y_unit = t[self.y].unit


        if 'log10phi_err_low' in t.colnames or 'phi_err_low' in t.colnames:
            uncertainties = True
        else:
            uncertainties = False


        self.log10X = {}
        self.log10X_binw = {}

        self.phi_dex = {}
        if uncertainties:
            self.phi_dex_err_upp = {}
            self.phi_dex_err_low = {}
            self.phi_dex_err = {}

        self.log10phi_dex = {}
        if uncertainties:
            self.log10phi_dex_err_upp = {}
            self.log10phi_dex_err_low = {}
            self.log10phi_dex_err = {}

        if self.x in ['M', 'log10L']:

            self.M = {}
            self.M_binw = {}

            self.phi_mag = {}
            if uncertainties:
                self.phi_mag_err_upp = {}
                self.phi_mag_err_low = {}
                self.phi_mag_err = {}

            self.log10phi_mag = {}
            if uncertainties:
                self.log10phi_mag_err_upp = {}
                self.log10phi_mag_err_low = {}
                self.log10phi_mag_err = {}



        # --- if data is redshift, log10L, phi ... this is most useful I think
        self.redshifts = list(set(self.t['z'].data))
        self.redshifts.sort()
        # --- make sure that redshift list is monotonically increasing
        if len(self.redshifts)>1:
            if self.redshifts[0]>self.redshifts[1]:
                self.redshifts = self.redshifts[::-1]

        # --- extract luminosities or magnitudes and number densities, making conversions where necessary.
        for z in self.redshifts:

            s = self.t['z'] == z

            # --- calculate log10X and M (if necessary)
            if self.x in ['log10L','log10Mstar','log10SFR']:
                self.log10X[z] = self.t[self.x][s].data

                if 'delta'+self.x in self.t.colnames:
                    self.log10X_binw[z] = self.t['delta'+self.x][s].data
                else:
                    self.log10X_binw[z] = (self.log10X[z][1]-self.log10X[z][0])*np.ones(len(self.log10X[z]))

                if self.x in ['M', 'log10L']:
                    self.M[z] = log10Lnu_to_M(self.log10X[z])
                    self.M_binw[z] = self.log10X_binw[z]*2.5


            elif self.x == 'M':
                self.M[z] = self.t['M'][s].data
                self.log10X[z] = M_to_log10Lnu(self.M[z])

                if 'delta'+self.x in self.t.colnames:
                    self.M_binw[z] = self.t['delta'+self.x][s].data
                else:
                    self.M_binw[z] = (self.M[z][1]-self.M[z][0])*np.ones(len(self.M[z]))

                self.log10X_binw[z] = self.M_binw[z]/2.5


            else:
                print('WARNING [x]')

            # --- calculate log10X and M (if necessary)

            if self.y == 'phi':
                phi = self.t['phi'][s].data
                log10phi = np.log10(phi)

                if uncertainties:
                    phi_err_low = self.t['phi_err_low'][s].data
                    phi_err_upp = self.t['phi_err_upp'][s].data
                    log10phi_err_low = np.log10(phi)-np.log10(phi-phi_err_low)
                    log10phi_err_upp = np.log10(phi+phi_err_upp)-np.log10(phi)
                    self.log10phi_dex_err_low[z] = log10phi_err_low
                    self.log10phi_dex_err_upp[z] = log10phi_err_upp
                    self.log10phi_dex_err[z] = [log10phi_err_low, log10phi_err_upp]
                    self.log10phi_mag_err_low[z] = log10phi_err_low
                    self.log10phi_mag_err_upp[z] = log10phi_err_upp
                    self.log10phi_dex_err[z] = [self.log10phi_mag_err_low[z], self.log10phi_mag_err_upp[z]]

            elif self.y == 'log10phi':
                log10phi = self.t['log10phi'][s].data
                phi = 10**log10phi

                if uncertainties:
                    self.log10phi_dex_err_low[z] = self.t['log10phi_err_low'][s].data
                    self.log10phi_dex_err_upp[z] = self.t['log10phi_err_upp'][s].data
                    self.log10phi_dex_err[z] = [self.log10phi_dex_err_low[z], self.log10phi_dex_err_upp[z]]
                    self.log10phi_mag_err_low[z] = self.t['log10phi_err_low'][s].data
                    self.log10phi_mag_err_upp[z] = self.t['log10phi_err_upp'][s].data
                    self.log10phi_dex_err[z] = [self.log10phi_mag_err_low[z], self.log10phi_mag_err_upp[z]]

                    # -- should add the inverse (e.g. phi_mag_err)

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


            if len(self.log10X[z])>1:

                if self.log10X[z][1]<self.log10X[z][0]:

                    self.M[z] = self.M[z][::-1]
                    self.log10X[z] = self.log10X[z][::-1]
                    self.log10phi_dex[z] = self.log10phi_dex[z][::-1]

                    if uncertainties:
                        self.log10phi_dex_err_upp[z] = self.log10phi_dex_err_upp[z][::-1]
                        self.log10phi_dex_err_low[z] = self.log10phi_dex_err_low[z][::-1]
                        self.log10phi_dex_err[z][0] = self.log10phi_dex_err[z][0][::-1]
                        self.log10phi_dex_err[z][1] = self.log10phi_dex_err[z][1][::-1]


        # --- short hand
        self.log10phi = self.log10phi_dex
        if uncertainties: self.log10phi_err = self.log10phi_dex_err







class Plots:

    def df(dataset_z_, cmap = 'cmr.guppy', x_range = [27., 30.], y_range = [-7., 1.], data_dir = data_dir):

        if type(dataset_z_) is not list: dataset_z_ = [dataset_z_]

        df_type = dataset_z_[0][0].split('/')[0]





        fig, ax = simple_fig()

        colors = cmr.take_cmap_colors(cmap, len(dataset_z_))
        lss = ['-','--','-.',':']*5

        for dataset_z, color, ls in zip(dataset_z_, colors, lss):

            dataset, z = dataset_z

            m = read(dataset, data_dir = data_dir)

            if m.df_type == 'binned':
                ax.plot(m.log10X[z], m.log10phi_dex[z], color = color, label = rf'$\rm {m.name} $', ls = ls)

            if m.df_type == 'schechter':
                log10L = np.arange(*x_range, 0.01)
                ax.plot(bin_centres(log10L), m.L(z).log10phi_binned(log10L), color = color, label = rf'$\rm {m.name} $', ls = ls)


        ax.legend(fontsize = 8)
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)

        ax.set_xlabel(label_(df_type))
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ dex^{-1}) $')

        return fig, ax
