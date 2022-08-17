

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u
from flare.LF.literature import UV


from flare.stats import poisson_confidence_interval

nm = 'donnan22'
d = np.loadtxt(f'{nm}.tex').T


z = 10
xname, xunit = 'M', 'mag'
yname, yunit = 'phi', 'Mpc^-3 mag^-1'

ones = np.ones(len(d[0]))



N = (d[3]/d[4])**2
print(N)
N = np.round(N)

print(poisson_confidence_interval(1.0))


phi_err_low = []
phi_err_upp = []

for phi, n in zip(d[3]*1E-6, N):

    ci = poisson_confidence_interval(n)

    phi_err_low.append((n-ci[0])*phi)
    phi_err_upp.append((ci[1]-n)*phi)



table = Table()
table.add_column(Column(data = d[0], name = 'z'))
table.add_column(Column(data = d[1], name = xname, unit = xunit))
table.add_column(Column(data = d[2], name = 'deltaM', unit = xunit))
table.add_column(Column(data = N, name = 'N'))
table.add_column(Column(data = d[3]*1E-6, name = yname, unit = yunit))
# table.add_column(Column(data = d[4]*1E-6, name = 'phi_err_low', unit = yunit))
# table.add_column(Column(data = d[4]*1E-6, name = 'phi_err_upp', unit = yunit))
table.add_column(Column(data = phi_err_low, name = 'phi_err_low', unit = yunit))
table.add_column(Column(data = phi_err_upp, name = 'phi_err_upp', unit = yunit))

table.meta['x'] = 'M'
table.meta['y'] = 'phi'
table.meta['name'] = 'Donnan et al. (2022)'
table.meta['type'] = 'binned'
table.meta['redshifts'] = list(set(table['z'].data))
table.meta['references'] = ['2022arXiv220712356D']


table.write(f'../binned/{nm}.ecsv', overwrite = True)
