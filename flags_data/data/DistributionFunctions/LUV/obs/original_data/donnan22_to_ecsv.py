

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u
from flare.LF.literature import UV



nm = 'donnan22'
d = np.loadtxt(f'{nm}.tex').T


z = 10
xname, xunit = 'M', 'mag'
yname, yunit = 'phi', 'Mpc^-3 mag^-1'

ones = np.ones(len(d[0]))

table = Table()
table.add_column(Column(data = d[0], name = 'z'))
table.add_column(Column(data = d[1], name = xname, unit = xunit))
table.add_column(Column(data = d[2], name = 'deltaM', unit = xunit))
table.add_column(Column(data = d[3]*1E-6, name = yname, unit = yunit))
table.add_column(Column(data = d[4]*1E-6, name = 'phi_err_low', unit = yunit))
table.add_column(Column(data = d[4]*1E-6, name = 'phi_err_upp', unit = yunit))

table.meta['x'] = 'M'
table.meta['y'] = 'phi'
table.meta['name'] = 'Donnan et al. (2022)'
table.meta['type'] = 'binned'
table.meta['redshifts'] = list(set(table['z'].data))
table.meta['references'] = ['2022arXiv220712356D']


table.write(f'../binned/{nm}.ecsv', overwrite = True)
