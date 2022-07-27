
import os

import numpy as np
from astropy.table import Table, Column, vstack, join
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import astropy.units as u



fn = 'chen22'

t1 = ascii.read(f'{fn}_tb1.txt')
t2 = ascii.read(f'{fn}_tb2.txt')

table = join(t1, t2, keys='id')

table.meta['name'] = f'Chen et al. (2022)'
table.meta['references'] = ['2022arXiv220712657C']
table.meta['quantities'] = ['M_UV', 'log10Mstar', 'log10sSFR', 'age']
table.meta['observatory'] = 'Webb'

table.write(f'../obs/{fn}.ecsv', format = 'ascii.ecsv', overwrite=True)
