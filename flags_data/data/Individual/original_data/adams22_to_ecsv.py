
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import astropy.units as u



hdul = fits.open('SMACS0724_MASTER_Sel-f200W_v2_hz3zp.fits')

hdul.info()

t = Table.read('SMACS0724_MASTER_Sel-f200W_v2_hz3zp.fits')

print(t.colnames)

#
#
# table = Table()
#
# table.meta['name'] = f'{lab}'
# table.meta['quantities'] = ['log10Mstar', 'sSFR', 'age']
# table.meta['references'] = []
#
# table.add_column(Column(data = t1['ID'], name = 'id'))
# table.add_column(Column(data = 7*np.ones(len(t1['ID'])), name = 'z')) # wrong! should cross match with Endsley
# table.add_column(Column(data = t1['logmass_'+mod], name = 'log10Mstar', unit = 'dex(solMass)'))
# table.add_column(Column(data = t1['age_'+mod], name = 'age', unit = 'Myr'))
# table.add_column(Column(data = t2['ssfr_'+modssfr], name = 'sSFR', unit = 'Gyr^-1'))
#
# table.write(f'obs/whitler22_{mod}.ecsv', format = 'ascii.ecsv', overwrite=True)
