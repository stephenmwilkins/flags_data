
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'Furlong_2015'



redshifts = [4,5,6,7]



tables = []

for z in redshifts:

    # print(z, f)

    data = np.loadtxt(f'original_data/{indir}/smf_z{z}.dat').T

    X = data[0]
    Y = data[1]
    s = ~np.isnan(Y)

    if np.sum(s)>0:
        t = Table()
        t.add_column(Column(data = z*np.ones(np.sum(s)), name = 'z'))
        t.add_column(Column(data = X[s], name = 'log10Mstar', unit = 'dex(Msun)'))
        t.add_column(Column(data = Y[s], name = 'log10phi', unit = 'dex(Mpc^-3 dex^-1)'))

        tables.append(t)

table = vstack(tables)

table.meta['x'] = 'log10Mstar'
table.meta['y'] = 'log10phi'
table.meta['name'] = 'EAGLE'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2015MNRAS.450.4486F','2015MNRAS.446..521S']

out_name = 'eagle'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
