
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'THESAN/UVLF'

redshifts = [6,7,8,9,10,11,12,13,14]

tables = []

for z in redshifts:

    # print(z, f)

    data = np.loadtxt(f'../../original_data/{indir}/z_{z}.txt').T

    M = data[0]
    phi = data[1]
    s = phi>0

    if np.sum(s)>0:
        t = Table()
        t.add_column(Column(data = z*np.ones(np.sum(s)), name = 'z'))
        t.add_column(Column(data = M[s], name = 'M'))
        t.add_column(Column(data = phi[s], name = 'phi', unit = 'Mpc^-3 mag^-1'))

        tables.append(t)

table = vstack(tables)

table.meta['x'] = 'M'
table.meta['y'] = 'phi'
table.meta['name'] = 'THESAN'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = []

out_name = 'thesan'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
