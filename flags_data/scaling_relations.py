

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



data_dir = f'{this_dir}/data/ScalingRelations'


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




class read:

    def __init__(self, dataset, data_dir = data_dir):

        t = Table.read(f'{data_dir}/{dataset}.ecsv')
        self.t = t

        print(dataset)

        self.name = t.meta['name']

        if 'references' in t.meta:
            self.references = t.meta['references']
        else:
            self.references = None

        self.dataset_type = t.meta['dataset_type']


        self.x, self.y, self.om, self.study = dataset.split('/')

        # --- grab default units
        self.x_unit = default_units[self.x]
        self.y_unit = default_units[self.y]


        # --- need to  handle unit conversions properly
        self._x = t.meta['x']
        self._y = t.meta['y']

        # self._x_unit = t[self.x].unit
        # self._y_unit = None


        self.redshifts = list(set(t['z'].data))

        self.X = {}
        self.Y = {}
        if 'N' in t.colnames: self.N = {}

        for z in self.redshifts:
            s = (t['z'].data == z)
            self.X[z] = t[self._x].data[s]
            self.Y[z] = {}
            if 'N' in t.colnames: self.N[z] = t['N'].data[s]

            if f'{self._y}_P50.0' not in self.t.colnames:
                self.Y[z][50.0] = self.t[f'{self._y}'].data[s]

            for p in [2.2, 15.8, 50.0, 84.2, 97.8]:
                if f'{self._y}_P{p}' in self.t.colnames:
                    self.Y[z][p] = self.t[f'{self._y}_P{p}'].data[s]
                    if not self.y_unit: self.y_unit = self.t[f'{self._y}_P{p}'].unit


    def plot_single_z(self, z, color = 'k'):

        fig, ax = simple_fig()

        if 2.2 in self.Y[z].keys():
            ax.fill_between(self.X[z], self.Y[z][2.2], self.Y[z][97.8], alpha = 0.1, color = color)
        if 15.8 in self.Y[z].keys():
            ax.fill_between(self.X[z], self.Y[z][15.8], self.Y[z][84.2], alpha = 0.2, color = color)

        ax.plot(self.X[z], self.Y[z][50.0], color = color)

        ax.set_ylabel(label(self.y, self.y_unit))
        ax.set_xlabel(label(self.x, self.x_unit))

        return fig, ax

    def plot_z_evo(self, cmap = 'cmr.guppy'):

        fig, ax = simple_fig()

        colors = cmr.take_cmap_colors(cmap, len(self.redshifts))

        for z, color in zip(self.redshifts, colors):
            ax.plot(self.X[z], self.Y[z][50.0], color = color, label = rf'$\rm z={z}$')

        ax.legend(fontsize = 8)
        ax.set_ylabel(label(self.y, self.y_unit))
        ax.set_xlabel(label(self.x, self.x_unit))

        return fig, ax



# def list_datasets(datasets, data_dir = f'{this_dir}/data/DistributionFunctions'):
#     return list(map(lambda x: x.replace('.ecsv', ''), os.listdir(f'{data_dir}/{datasets}')))


def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f

def list_relations(data_dir = data_dir):
    relations = {}
    for x in listdir_nohidden(f'{data_dir}'):
        if os.path.isdir(os.path.join(data_dir, x)):
            relations[x] = listdir_nohidden(f'{data_dir}/{x}')
    return relations


def list_datasets(relation = '', data_dir = data_dir):
    l = [os.path.join(dp, f.split('.')[0]) for dp, dn, fn in os.walk(os.path.expanduser(f'{data_dir}/{relation}')) for f in fn if f.endswith('.ecsv')]
    l = [l_[len(data_dir)+1:] for l_ in l]
    return l







class Datasets:

    def __init__(self, datasets = '', data_dir = f'{this_dir}/data/ScalingRelations'):

        self.data_dir = data_dir
        self.datasets = datasets
        self.dataset_list = list_datasets(datasets)
        self.n = len(self.dataset_list)
        self.relation = '/'.join(datasets.split('/')[:2])
        print(self.relation)
        self.x, self.y = self.relation.split('/')
        self.x_unit = default_units[self.x]
        self.y_unit = default_units[self.y]

        # print(f'{self.x}/{self.x_unit:latex}', f'{self.y}/{self.y_unit:latex}')


        self.type_studies = list(set(['/'.join(x.split('/')[-2:]) for x in self.dataset_list]))
        self.studies = list(set([x.split('/')[-1] for x in self.dataset_list]))

        self.ts_from_s = dict(zip(self.studies, self.type_studies))

        # --- read all datasets
        self.dataset = {}
        self.redshifts = {}
        self.log10X_range = {}
        self.names = {}

        for dataset_name in self.dataset_list:

            ds = read(dataset_name, data_dir = data_dir)

            self.dataset[dataset_name] = ds

            # --- get dataset name
            self.names[dataset_name] = ds.name

            # --- get redshifts of each models
            self.redshifts[dataset_name] = ds.redshifts

            # --- get log10X range of binned models
            if dataset_name.split('/')[-1] == 'binned':
                self.log10X_range[dataset_name] = [np.min(ds.X[ds.redshifts[0]]), np.max(ds.X[ds.redshifts[0]])]







    def get_datasets_at_z(self, z, z_tolerance = 0.2):
        datasets_z = []
        for dataset_name in self.dataset_list:
            for z_ in self.redshifts[dataset_name]:
                if np.fabs(z-z_)<z_tolerance:
                    datasets_z.append((dataset_name, z_))
        return(datasets_z)


    def plot_redshift_range(self, cmap = 'cmr.guppy', add_references = False):

        xtotal_ = 7.  # in "

        bottom_ = 0.45  # in "
        height_ = 0.3 * self.n  # in "
        top_ = 0.1  # in "
        ytotal_ = bottom_ + height_ + top_  # in "

        left  = 0.4
        width = 0.55

        bottom = bottom_/ytotal_
        height = height_/ytotal_

        fig = plt.figure(figsize = (xtotal_, ytotal_))

        ax = fig.add_axes((left, bottom, width, height))



        colors = cmr.take_cmap_colors(cmap, self.n)

        z_min, z_extent, labels = [], [], []

        for i, ts in enumerate(sorted(self.type_studies)):

            ds = self.dataset[f'{self.relation}/{ts}']

            label = rf'$\rm \mathbf{{ {ds.name} }}$'
            # ax.text(3.0, i + 0.2, label, fontsize = 8, ha='right')
            ax.text(3.3, i, label, fontsize = 8, ha='right', va = 'center')

            # --- add references
            if add_references and ds.references:
                articles = [list(ads.SearchQuery(bibcode=bibcode))[0] for bibcode in ds.references]
                refs = [f"{article.first_author.split(',')[0]}+{article.year}"  for article in articles]
                ax.text(3.0, i - 0.2, ', '.join(refs), fontsize = 6, ha='right', color = '0.5')

            z_min.append(np.min(ds.redshifts))
            z_extent.append(np.max(ds.redshifts)-np.min(ds.redshifts))


        ax.barh(np.arange(self.n), z_extent, left = z_min, color = colors, align='center')

        ax.set_xlim([3.5, 15.5])
        ax.set_ylim([-0.75, self.n - 0.25])
        ax.set_yticks([])
        ax.set_xlabel(r'$\rm z$')

        return fig, ax


    def plot_redshift_X_range(self, cmap = 'cmr.guppy', threshold = 1):

        fig, ax = simple_fig(fig_size = (4.5, 3.5))

        for ds_name, color in zip(self.dataset_list, cmr.take_cmap_colors(cmap, self.n)):

            ds = self.dataset[ds_name]

            x, y1, y2 = [], [], []

            for z in ds.redshifts:
                x.append(z)
                # s = (m.N[z]>threshold)
                # y1.append(np.min(m.X[z][s]))
                # y2.append(np.max(m.X[z][s]))
                y1.append(np.min(ds.X[z]))
                y2.append(np.max(ds.X[z]))

            ax.fill_between(x,y1,y2, color = color, alpha =0.2)
            ax.plot(x + x[::-1] + [x[0]], y1 + y2[::-1] + [y1[0]], lw = 1, color = color, zorder = 3)

            ax.text(x[0] + 0.1, y1[0]+0.4, rf'$\rm {ds.name}$', fontsize = 10,  rotation = 90., color = color)

        ax.set_ylabel(label(ds.x, ds.x_unit))
        ax.set_xlabel(r'$\rm z$')

        return fig, ax



    def plot_srs(self, redshifts = np.arange(5, 14, 1), cmap = 'cmr.guppy'):

        """ plot all datasets at a range of redshifts """


        fig = plt.figure(figsize = (7,5))

        left = 0.1
        right = 0.75
        bottom = 0.1
        top = 0.95

        gs = fig.add_gridspec(3, 3, left = left, bottom = bottom, right = right, top = top, hspace=0, wspace=0)
        axes = gs.subplots(sharex=True, sharey=True)


        colors = dict(zip(self.dataset_list, cmr.take_cmap_colors(cmap, self.n)))
        lss = dict(zip(self.dataset_list, ['-','--','-.',':']*5))

        # --- create legend
        lax = fig.add_axes([right, bottom, 0.3, top-bottom])
        lax.axis('off')
        handles = [Line2D([0], [0], color = colors[ds], lw=2, ls=lss[ds], label=self.names[ds]) for ds in self.dataset_list]
        lax.legend(handles=handles, loc='center left', title = self.datasets, fontsize = 8)

        for ax, z in zip(axes.flatten(), redshifts):

            ax.label_outer()

            ax.text(0.05, 0.9, rf'$\rm z={z:.0f}$', color = '0.3', transform = ax.transAxes)

            dataset_z = self.get_datasets_at_z(z, z_tolerance = 0.1) # --- get list of datasets at this redshift

            for dataset_name, z in dataset_z:

                ds = self.dataset[dataset_name]

                color = colors[dataset_name]
                ls  = lss[dataset_name]

                ax.plot(ds.X[z], ds.Y[z][50.], color = color, label = rf'$\rm {ds.name} $', ls = ls)


            # ax.set_xlim(x_range)
            # ax.set_ylim(y_range)
            # ax.set_xticks(np.arange(*np.round(x_range, 0), 1))

            # ax.legend(fontsize = 8)

        axes[2,1].set_xlabel(label(self.x, self.x_unit), fontsize = 10)
        axes[1,0].set_ylabel(label(self.y, self.y_unit), fontsize = 10)

        return fig, axes
