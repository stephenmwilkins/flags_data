
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u



indir = 'Astrid/Distribution/UVLF'

redshifts = [3,4,5,6,7,8,9,10]

tables = []

for z in redshifts:

    # print(z, f)

    data = np.loadtxt(f'../../original_data/{indir}/UVLF_int_z{z}.dat').T

    M = data[0]
    phi = data[1]

    t = Table()
    t.add_column(Column(data = z*np.ones(len(M)), name = 'z'))
    t.add_column(Column(data = M, name = 'M', unit = 'mag'))
    t.add_column(Column(data = phi, name = 'log10phi', unit = 'dex(Mpc^-3 mag^-1)'))
    tables.append(t)

table = vstack(tables)

table.meta['x'] = 'M'
table.meta['y'] = 'log10phi'
table.meta['name'] = 'Astrid'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2022MNRAS.512.3703B']

out_name = 'astrid'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
