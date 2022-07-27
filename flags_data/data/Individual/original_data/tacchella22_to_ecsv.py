
import os

import numpy as np
from astropy.table import Table, Column, vstack
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import astropy.units as u


 # AV_q50 AV_q16 AV_q84 log_stellar_mass_q50 log_stellar_mass_q16 log_stellar_mass_q84 log_sfr_10_q50 log_sfr_10_q16 log_sfr_10_q84 log_ssfr_10_q50 log_ssfr_10_q16 log_ssfr_10_q84 log_sfr_50_q50 log_sfr_50_q16 log_sfr_50_q84 log_ssfr_50_q50 log_ssfr_50_q16 log_ssfr_50_q84 log_sfr_100_q50 log_sfr_100_q16 log_sfr_100_q84 log_ssfr_100_q50 log_ssfr_100_q16 log_ssfr_100_q84 dust2_q50 dust2_q16 dust2_q84 dust_index_q50 dust_index_q16 dust_index_q84 dust1_fraction_q50 dust1_fraction_q16 dust1_fraction_q84 time_50_q50 time_50_q16 time_50_q84 logzsol_q50 logzsol_q16 logzsol_q84 mag_1500_q50 mag_1500_q16 mag_1500_q84 mag_1500_intrinsic_q50 mag_1500_intrinsic_q16 mag_1500_intrinsic_q84 log_nion_q50 log_nion_q16 log_nion_q84 log_xion_q50 log_xion_q16 log_xion_q84


fn = 'tacchella22'

dkey = {}

dkey['z'] = 'redshift'
dkey['log10Mstar'] = 'log_stellar_mass'
dkey['log10sSFR'] = 'log_ssfr_50'
dkey['log10SFR'] = 'log_sfr_50'
# dkey['age'] = 'time_50'

conversion = {k: 0.0 for k in dkey.keys()}
conversion['log10sSFR'] = 9.


data = ascii.read(f'{fn}.txt')

table = Table()

table.meta['name'] = f'Tacchella et al. (2022)'
table.meta['quantities'] = ['log10Mstar', 'log10SFR', 'log10sSFR', 'age']
table.meta['references'] = []
table.meta['observatory'] = 'Webb'


table.add_column(Column(data = data['id'], name = 'id'))

for k,v in dkey.items():
    table.add_column(Column(data = data[v+'_q50'] + conversion[k], name = k))
    table.add_column(Column(data = data[v+'_q50']-data[v+'_q16'], name = k+'_err_low'))
    table.add_column(Column(data = data[v+'_q84']-data[v+'_q50'], name = k+'_err_upp'))

k, v = 'age', 'time_50'
table.add_column(Column(data = 1E3*data[v+'_q50'], name = k))
table.add_column(Column(data = 1E3*(data[v+'_q50']-data[v+'_q16']), name = k+'_err_low'))
table.add_column(Column(data = 1E3*(data[v+'_q84']-data[v+'_q50']), name = k+'_err_upp'))


table.write(f'../obs/{fn}.ecsv', format = 'ascii.ecsv', overwrite=True)
