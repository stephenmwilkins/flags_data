
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'behroozi_silk_2015/lf_nodust'

files = os.listdir(f'../LUV/original_data/{indir}')

files = sorted(files)

redshifts = [10, 12.5, 14, 15]

tables = []

for z, f in zip(redshifts, files):

    # print(z, f)

    data = np.loadtxt(f'../LUV/original_data/{indir}/{f}').T

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
table.meta['name'] = 'Behroozi and Silk (2015)'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2015ApJ...799...32B']

out_name = 'behroozi_silk2015'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
