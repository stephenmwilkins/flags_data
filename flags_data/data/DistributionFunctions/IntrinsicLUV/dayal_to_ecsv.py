
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'dayal'

redshifts = [5,6,7,8,9,10,11,12,13,14,15,16,18,20]

tables = []

for z_ in redshifts:

    # print(z, f)

    data = np.loadtxt(f'../LUV/original_data/{indir}/Delphi_UVLF_nouv_int_z{z_}.dat').T

    z = data[0]
    print(set(data[0]))
    M = data[1]
    phi = data[2]

    t = Table()
    t.add_column(Column(data = z, name = 'z'))
    t.add_column(Column(data = M, name = 'M'))
    t.add_column(Column(data = phi, name = 'log10phi', unit = 'dex(Mpc^-3 mag^-1)'))
    tables.append(t)

table = vstack(tables)

table.meta['x'] = 'M'
table.meta['y'] = 'log10phi'
table.meta['name'] = 'Delphi'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2014MNRAS.445.2545D', '2022MNRAS.512..989D']

out_name = 'delphi'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
