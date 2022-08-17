
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import astropy.units as u

from flare import photom


fn = 'adams22'

d = ascii.read(f'{fn}.tex')
table = Table()

for c in d.colnames:
    table.add_column(Column(data = d[c], name = c))


z = d['z'].data
M_UV = np.zeros(len(d['id']))


M_UV = d['F200W'] - photom.DM(z)
print(M_UV)

table.add_column(Column(data = M_UV, name = 'M_UV', unit = 'mag'))


table.meta['name'] = f'Adams et al. (2022)'
table.meta['references'] = []
table.meta['quantities'] = []
table.meta['observatory'] = 'Webb'

table.write(f'../obs/{fn}.ecsv', format = 'ascii.ecsv', overwrite=True)
