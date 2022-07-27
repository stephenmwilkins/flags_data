
import os

import numpy as np
from astropy.table import Table, Column, vstack, join
import astropy.io.ascii as ascii
from astropy import units as u
from astropy.coordinates import SkyCoord


fn = 'yan22'


d = ascii.read(f'{fn}.tex')
table = Table()

for c in d.colnames:
    if c not in ['radec']:
        table.add_column(Column(data =d[c], name = c))


name = np.array(d['radec'].data).astype(str)

ra = [SkyCoord.from_name(nm).ra.deg for nm in name]
dec = [SkyCoord.from_name(nm).dec.deg for nm in name]

table.add_column(Column(data = ra, name = 'ra', unit = 'deg'))
table.add_column(Column(data = dec, name = 'dec', unit = 'deg'))

table.add_column(Column(data = ['SMACS']*len(name), name = 'field'))


table.meta['name'] = f'Yan et al. (2022)'
table.meta['references'] = []
table.meta['quantities'] = []
table.meta['observatory'] = 'Webb'

table.write(f'../obs/{fn}.ecsv', format = 'ascii.ecsv', overwrite=True)
