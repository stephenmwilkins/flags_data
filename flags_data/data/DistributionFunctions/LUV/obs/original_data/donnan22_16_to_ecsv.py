

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u
from flare.LF.literature import UV


from astropy.cosmology import Planck18 as cosmo

z1, z2 = 10, 13
z1, z2 = 15, 20

vol = cosmo.comoving_volume(z2) - cosmo.comoving_volume(z1)
print(vol)

area = 49.*u.arcmin*u.arcmin # naidu
area = 45.*u.arcmin*u.arcmin # donnen

vol_survey = vol*area.to('sr')/(4*np.pi*u.sr)

# print(vol_survey)
# print(vol_survey/1E5)

nm = 'donnan22_16'



xname, xunit = 'M', 'mag'
yname, yunit = 'phi', 'Mpc^-3 mag^-1'



phi = 1./vol_survey.to('Mpc^3').value


table = Table()
table.add_column(Column(data = [16], name = 'z'))
table.add_column(Column(data = [-21.5], name = xname, unit = xunit))
table.add_column(Column(data = [1], name = 'deltaM', unit = xunit))
table.add_column(Column(data = [phi], name = yname, unit = yunit))
table.add_column(Column(data = [(1-0.174)*phi], name = 'phi_err_low', unit = yunit))
table.add_column(Column(data = [(3.289-1)*phi], name = 'phi_err_upp', unit = yunit))

table.meta['x'] = 'M'
table.meta['y'] = 'phi'
table.meta['name'] = 'Donnan et al. (2022) (z=15-17)'
table.meta['type'] = 'binned'
table.meta['redshifts'] = list(set(table['z'].data))
table.meta['references'] = ['2022arXiv220712356D']


table.write(f'../binned/{nm}.ecsv', overwrite = True)
