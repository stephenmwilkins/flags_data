
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'Behroozi_2019'



redshifts = [4,5,6,7,8,9,10]



tables = []

for z in redshifts:

    # print(z, f)

    data = np.loadtxt(f'original_data/{indir}/behroozi_smf_z{z}.dat').T

    X = data[0]
    Y = data[1]
    s = ~np.isnan(Y)

    if np.sum(s)>0:
        t = Table()
        t.add_column(Column(data = z*np.ones(np.sum(s)), name = 'z'))
        t.add_column(Column(data = X[s], name = 'log10Mstar', unit = 'dex(Msun)'))
        t.add_column(Column(data = Y[s], name = 'phi', unit = 'Mpc^-3 dex^-1'))

        tables.append(t)

table = vstack(tables)

table.meta['x'] = 'log10Mstar'
table.meta['y'] = 'phi'
table.meta['name'] = 'Universe Machine'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2019MNRAS.488.3143B']

out_name = 'universe_machine'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
