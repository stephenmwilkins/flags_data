
import os

import numpy as np
from astropy.table import Table, Column, vstack, join
import astropy.io.ascii as ascii
from astropy import units as u
from astropy.coordinates import SkyCoord


fn = 'harikane22'


d = ascii.read(f'{fn}.tex')
table = Table()

for c in d.colnames:
    if c not in ['ra', 'dec']:
        table.add_column(Column(data = d[c], name = c))


c = SkyCoord(d['ra'], d['dec'], unit=(u.hourangle, u.deg))


table.add_column(Column(data = c.ra.deg, name = 'ra', unit = 'deg'))
table.add_column(Column(data = c.dec.deg, name = 'dec', unit = 'deg'))


table.meta['name'] = f'Harikane et al. (2022)'
table.meta['references'] = []
table.meta['quantities'] = []
table.meta['observatory'] = 'Webb'

table.write(f'../obs/{fn}.ecsv', format = 'ascii.ecsv', overwrite=True)
