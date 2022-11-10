
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u


indir = 'scsam'

redshifts = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

tables = []

for z_ in redshifts:

    # print(z, f)

    data = np.loadtxt(f'original_data/{indir}/UVLF_nd_z{z_}.dat').T

    M = data[0]
    phi = data[1]
    z = np.ones(len(M))*z_

    t = Table()
    t.add_column(Column(data=z, name='z'))
    t.add_column(Column(data=M, name='M'))
    t.add_column(Column(data=phi, name='phi', unit='Mpc^-3 mag^-1'))
    tables.append(t)

table = vstack(tables)

table.meta['x'] = 'M'
table.meta['y'] = 'log10phi'
table.meta['name'] = 'Santa Cruz SAM'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = []

out_name = 'scsam'

table.write(f'models/binned/{out_name}.ecsv', format='ascii.ecsv', overwrite=True)
