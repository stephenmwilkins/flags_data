

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

import ads

ads.config.token = 'qm1AtsIgKukl0jqMjYaEa2LHK9am6gQka6opvce1'

data_dir = f'{this_dir}/data/ScalingRelations'




class read:

    def __init__(self, dataset, data_dir = data_dir):

        t = Table.read(f'{data_dir}/{dataset}.ecsv')
        self.t = t

        self.name = t.meta['name']

        if 'references' in t.meta:
            self.references = t.meta['references']
        else:
            self.references = None

        self.x = t.meta['x']
        self.y = t.meta['y']

        self.x_unit = t[self.x].unit
        self.y_unit = None


        self.redshifts = list(set(t['z'].data))

        self.X = {}
        self.Y = {}
        if 'N' in t.colnames: self.N = {}
        for z in self.redshifts:
            s = (t['z'].data == z)
            self.X[z] = t[self.x].data[s]
            self.Y[z] = {}
            if 'N' in t.colnames: self.N[z] = t['N'].data[s]

            if f'{self.y}_P50.0' not in self.t.colnames:
                self.Y[z][50.0] = self.t[f'{self.y}'].data[s]

            for p in [2.2, 15.8, 50.0, 84.2, 97.8]:
                if f'{self.y}_P{p}' in self.t.colnames:
                    self.Y[z][p] = self.t[f'{self.y}_P{p}'].data[s]
                    if not self.y_unit: self.y_unit = self.t[f'{self.y}_P{p}'].unit


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


def list_relations(data_dir = data_dir):
    relations = {}
    for x in os.listdir(f'{data_dir}'):
        if os.path.isdir(os.path.join(data_dir, x)) and x != '.DS_Store':
            relations[x] = os.listdir(f'{data_dir}/{x}')
    return relations


def list_datasets(relation = '', data_dir = data_dir):
    l = [os.path.join(dp, f.split('.')[0]) for dp, dn, fn in os.walk(os.path.expanduser(f'{data_dir}/{relation}')) for f in fn if f.endswith('.ecsv')]
    l = [l_[len(data_dir)+1:] for l_ in l]
    return l





class DatasetInfo:

    def __init__(self, datasets = '', data_dir = f'{this_dir}/data/ScalingRelations'):

        self.data_dir = data_dir
        self.datasets = datasets
        self.dataset_list = list_datasets(datasets)
        self.relations = list(set(['/'.join(x.split('/')[:-2]) for x in self.dataset_list]))
        self.type_studies = list(set(['/'.join(x.split('/')[-2:]) for x in self.dataset_list]))
        self.studies = list(set([x.split('/')[-1] for x in self.dataset_list]))

        self.ts_from_s = dict(zip(self.studies, self.type_studies))

        for dset in self.dataset_list:

            m = read(dset, data_dir = data_dir)

            log10_limits = [np.min(m.X[m.redshifts[0]]), np.max(m.X[m.redshifts[0]])]

            print(dset, m.redshifts, log10_limits)


    def plot_matrix(self, cmap = 'cmr.horizon_r'):

        fig = plt.figure(figsize = (5, 4.5))

        left  = 0.4
        height = 0.8
        bottom = 0.15
        width = 0.55

        ax = fig.add_axes((left, bottom, width, height))

        M = np.zeros([len(self.relations), len(self.studies)])

        study_name = {}

        for i, study in enumerate(sorted(self.studies)):
            for j, relation in enumerate(self.relations):
                dset = f'{relation}/{self.ts_from_s[study]}'
                if dset in self.dataset_list:
                    m = read(dset, data_dir = self.data_dir)
                    study_name[study] = rf'$\rm \mathbf{{ {m.name} }}$'
                    M[j, i] = (i+1)/len(self.studies)

        ax.imshow(M, cmap = cmap, origin = 'lower')

        ax.set_yticks(np.arange(len(self.relations)), labels = self.relations)
        ax.set_xticks(np.arange(len(self.studies)), labels = [study_name[study] for study in sorted(self.studies)])

        # ax.set_xlabel(r'$\rm z$')

        return fig, ax



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

            m = read(self.long_from_short[short_model_name], data_dir = self.data_dir)

            label = rf'$\rm \mathbf{{ {m.name} }}$'
            ax.text(3.0, i + 0.2, label, fontsize = 8, ha='right')

            if m.references:

                articles = [list(ads.SearchQuery(bibcode=bibcode))[0] for bibcode in m.references]

                refs = [f"{article.first_author.split(',')[0]}+{article.year}"  for article in articles]

                ax.text(3.0, i - 0.2, ', '.join(refs), fontsize = 6, ha='right', color = '0.5')

            z_min.append(np.min(m.redshifts))
            z_extent.append(np.max(m.redshifts)-np.min(m.redshifts))

        ax.barh(np.arange(len(self.models)), z_extent, left = z_min, color = colors, align='center')

        ax.set_xlim([3.5, 15.5])
        ax.set_ylim([-0.5, np.max([5, len(self.models)])])
        ax.set_yticks([])
        ax.set_xlabel(r'$\rm z$')

        return fig, ax


    def plot_redshift_X_range(self, cmap = 'cmr.guppy', threshold = 1):

        fig, ax = simple_fig(fig_size = (4.5, 3.5))

        for model_name, color in zip(self.models, cmr.take_cmap_colors(cmap, len(self.models))):

            m = read(model_name, data_dir = self.data_dir)

            x, y1, y2 = [], [], []

            for z in m.redshifts:
                x.append(z)
                # s = (m.N[z]>threshold)
                # y1.append(np.min(m.X[z][s]))
                # y2.append(np.max(m.X[z][s]))
                y1.append(np.min(m.X[z]))
                y2.append(np.max(m.X[z]))

            ax.fill_between(x,y1,y2, color = color, alpha =0.2)
            ax.plot(x + x[::-1] + [x[0]], y1 + y2[::-1] + [y1[0]], lw = 1, color = color, zorder = 3)

            ax.text(x[0] + 0.1, y1[0]+0.4, rf'$\rm {m.name}$', fontsize = 10,  rotation = 90., color = color)

        ax.set_ylabel(label(m.x, m.x_unit))
        ax.set_xlabel(r'$\rm z$')

        return fig, ax
