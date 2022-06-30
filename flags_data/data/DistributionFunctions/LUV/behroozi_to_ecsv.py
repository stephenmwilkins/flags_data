
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'behroozi_2019'

files = os.listdir(f'original_data/{indir}/')

files = sorted(files)

scale_factors = np.array([np.float('0.'+fn.split('.')[1]) for fn in files])

redshifts = 1/(scale_factors) - 1

# print(scale_factors)
# print(redshifts)


tables = []

for z, f in zip(redshifts, files):

    # print(z, f)

    data = np.loadtxt(f'original_data/{indir}/{f}').T

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
table.meta['name'] = 'Universe\ Machine'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2020MNRAS.499.5702B']

out_name = 'universe_machine'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
