
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flags_data.utilities import bin_centres
import flags_data.scaling_relations as sr

import flare.plt as fplt



if __name__ == "__main__":

    # --- list available relations
    relations = sr.list_relations()

    # --- create a grid of relations and models
    # di = sr.Datasets()
    # fig, ax = di.plot_matrix()
    # fig.savefig(f'figs/sr/matrix.png')

    for x, Y in relations.items():
        for y in Y:
            r = f'{x}/{y}'
            # os.mkdir(f'figs/sr/{x}-{y}')
            print(r, sr.list_datasets(f'{x}/{y}'))

            ds = sr.Datasets(datasets = r)

            if ds.n > 0:

                fig, ax = ds.plot_redshift_range()
                fig.savefig(f'figs/sr/{x}-{y}/z_r.png')

                fig, ax = ds.plot_redshift_X_range()
                fig.savefig(f'figs/sr/{x}-{y}/z_X_r.png')

                fig, ax = ds.plot_srs()
                fig.savefig(f'figs/sr/{x}-{y}/sr.png')

                # --- make README

                # lines = []
                # lines.append(f'# {x}-{y}\n')
                # lines.append(f'![](../figs/sr/{x}-{y}/z_r.png)\n')
                # lines.append(f'![](../figs/sr/{x}-{y}/z_X_r.png)\n')
                # lines.append(f'![](../figs/sr/{x}-{y}/sr.png)\n')
                #
                # f = open(f'sr/{x}-{y}.md','w+').writelines(lines)
