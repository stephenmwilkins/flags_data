
import os

import numpy as np
from astropy.table import Table, Column, vstack, join
import astropy.io.ascii as ascii
from astropy import units as u
from astropy.coordinates import SkyCoord


fn = 'donnan22'





d = ascii.read(f'{fn}_radec.tex')
t1 = Table()

for c in 'id z z_err_upp z_err_low M_UV'.split():
    t1.add_column(Column(data = d[c], name = c))

c = SkyCoord(d['ra'], d['dec'], unit=(u.hourangle, u.deg))

t1.add_column(Column(data = c.ra.deg, name = 'ra', unit = 'deg'))
t1.add_column(Column(data = c.dec.deg, name = 'dec', unit = 'deg'))

t2 = Table()
d = ascii.read(f'{fn}_B2.tex')

for c in 'id field F090W F090W_err_upp F090W_err_low F115W F115W_err_upp F115W_err_low F150W  F150W_err_upp F150W_err_low F200W  F200W_err_upp F200W_err_low F277W  F277W_err_upp F277W_err_low F356W  F356W_err_upp F356W_err_low F410M  F410M_err_upp F410M_err_low F444W'.split():
    t2.add_column(Column(data = d[c], name = c))

print(len(t1))
print(len(t2))

table = join(t1, t2, keys='id')


table.meta['name'] = f'Donnan et al. (2022)'
table.meta['references'] = []
table.meta['quantities'] = ['M_UV']
table.meta['observatory'] = 'Webb'

table.write(f'../obs/{fn}.ecsv', format = 'ascii.ecsv', overwrite=True)
