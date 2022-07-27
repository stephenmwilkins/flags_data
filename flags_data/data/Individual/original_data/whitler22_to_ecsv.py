
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import astropy.units as u


hdul = fits.open(f'original_data/whitler22_age_mass_medians.fits')
t1 = hdul[1].data # assuming the first extension is a table
print(hdul[1].columns)
hdul = fits.open(f'original_data/whitler22_ssfr_medians.fits')
t2 = hdul[1].data # assuming the first extension is a table
print(hdul[1].columns)
hdul.close()
#


# t1 = ascii.read(f'original_data/whitler22_age_mass_medians.fits', encoding='utf-8')
# # t2 = ascii.read(f'original_data/whitler22_ssfr_medians.fits')
#
# print(t1)
#
# print(t1.colnames)
# print(t2.colnames)

modlabs = []
modlabs.append(('beagle_csfh',  'beagle_csfh','Whitler et al. (2022) - BEAGLE CSFH'))
modlabs.append(('prospect_csfh', 'prospect_csfh','W22 - Prospector CSFH'))
modlabs.append(('prospect_nonpar', 'prospect_nonpar100','W22 - Prospector non-par 100 Myr]'))


for mod, modssfr, lab in modlabs:

    table = Table()

    table.meta['name'] = f'{lab}'
    table.meta['quantities'] = ['log10Mstar', 'sSFR', 'age']
    table.meta['references'] = []

    table.add_column(Column(data = t1['ID'], name = 'id'))
    table.add_column(Column(data = 7*np.ones(len(t1['ID'])), name = 'z')) # wrong! should cross match with Endsley
    table.add_column(Column(data = t1['logmass_'+mod], name = 'log10Mstar', unit = 'dex(solMass)'))
    table.add_column(Column(data = t1['age_'+mod], name = 'age', unit = 'Myr'))
    table.add_column(Column(data = t2['ssfr_'+modssfr], name = 'sSFR', unit = 'Gyr^-1'))

    table.write(f'obs/whitler22_{mod}.ecsv', format = 'ascii.ecsv', overwrite=True)
