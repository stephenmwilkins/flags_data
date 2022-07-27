



import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u


redshifts = [5.0, 6.0]

log10Mstar = {}
log10Mstar_err = {}
log10SFR = {}
log10SFR_uppErr = {}
log10SFR_lowErr = {}
log10SFR_err = {}


log10Mstar[5] = np.array([9.0,9.25,9.50,9.75,10.,10.25])
log10SFR[5] = np.array([0.88,1.04,1.12,1.23,1.46,1.62])
log10SFR_uppErr[5] = log10SFR_lowErr[5] = np.array([0.42,0.38,0.41,0.43,0.31,0.37])

log10Mstar[6] = np.array([9.0,9.25,9.50,9.75,10.])
log10SFR[6] = np.array([0.92,1.07,1.27,1.40,1.147])
log10SFR_uppErr[6] = log10SFR_lowErr[6] = np.array([0.19,0.21,0.35,0.26,0.07])

tables = []

for z in redshifts:

    t = Table()
    t.add_column(Column(data = z*np.ones(len(log10Mstar[z])), name = 'z'))
    t.add_column(Column(data = log10Mstar[z], name = 'log10Mstar', unit = 'dex(Msun)'))

    v_ = log10SFR[z] - log10Mstar[z] + 9.0

    t.add_column(Column(data = v_, name = 'log10sSFR', unit = 'dex(1/Gyr)'))
    v = v_+log10SFR_uppErr[z]
    t.add_column(Column(data = v, name = 'log10sSFR_P15.8', unit = 'dex'))
    v = v_-log10SFR_uppErr[z]
    t.add_column(Column(data = v, name = 'log10sSFR_P84.2', unit = 'dex'))
    tables.append(t)

table = vstack(tables)

table.meta['x'] = 'log10Mstar'
table.meta['y'] = 'log10sSFR'
table.meta['name'] = 'Salmon et al. (2015)'
table.meta['redshifts'] = list(set(table['z'].data))
table.meta['references'] = []


filename = table.write(f'obs/salmon15.ecsv', format = 'ascii.ecsv', overwrite=True)
