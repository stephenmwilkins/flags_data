
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flags_data.utilities import bin_centres
import flags_data.scaling_relations as sr

import flare.plt as fplt

import ads

ads.config.token = 'qm1AtsIgKukl0jqMjYaEa2LHK9am6gQka6opvce1'


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

                # fig, ax = ds.plot_redshift_range()
                # fig.savefig(f'figs/sr/{x}-{y}/z_r.png')
                #
                # fig, ax = ds.plot_redshift_X_range()
                # fig.savefig(f'figs/sr/{x}-{y}/z_X_r.png')
                #
                # fig, ax = ds.plot_srs()
                # fig.savefig(f'figs/sr/{x}-{y}/sr.png')

                # --- make README

                lines = []
                lines.append(f'# {x}-{y}\n')
                lines.append(f'\n')

                # --- add references

                for dset_name, dset in ds.dataset.items():

                    if dset.references:
                        articles = [list(ads.SearchQuery(bibcode=bibcode))[0] for bibcode in dset.references]
                        refs = [f"{article.first_author.split(',')[0]}+{article.year}"  for article in articles]
                        links = [f'[{ref}](https://ui.adsabs.harvard.edu/abs/{bibcode}/abstract)' for ref, bibcode in zip(refs, dset.references)]
                        lines.append(f'| {dset.name} | {", ".join(links)} |\n')


                lines.append('\n')
                lines.append(f'![](../figs/sr/{x}-{y}/z_r.png)\n')
                lines.append(f'![](../figs/sr/{x}-{y}/z_X_r.png)\n')
                lines.append(f'![](../figs/sr/{x}-{y}/sr.png)\n')

                f = open(f'sr/{x}-{y}.md','w+').writelines(lines)
