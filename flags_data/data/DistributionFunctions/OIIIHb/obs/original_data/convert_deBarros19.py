
""" convert the deBarros+ individual sources into a GSMF """

import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column, vstack
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo
from flare.stats import poisson_confidence_interval

study_name = 'debarros19'

data = ascii.read(f'OIII_LF_z8_wSigmaInt_1SigmaContour.txt')


log10L = data['Lline']
log10phi_16 = data['LF16']
log10phi_84 = data['LF84']

log10phi = 0.5*(log10phi_16+log10phi_84)  # guess


table = Table()
table.add_column(Column(data=8.0*np.ones(len(log10L)), name='z'))
table.add_column(Column(data=log10L, name='log10L', unit='dex(erg/s)'))
table.add_column(Column(data=log10phi, name='log10phi', unit='dex(1 / (dex Mpc3))'))


table.meta['x'] = 'log10L'
table.meta['y'] = 'log10phi'
table.meta['name'] = 'de Barros et al. (2019)'
table.meta['type'] = 'binned'
table.meta['redshifts'] = list(set(table['z'].data))
table.meta['references'] = []

table.write(f'../binned/{study_name}.ecsv', overwrite=True)
