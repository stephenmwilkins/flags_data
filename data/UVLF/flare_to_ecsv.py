

import numpy as np
from astropy.table import Table, Column
import astropy.units as u
from flare.LF.literature import UV






# --- binned

run_model_binned = False

if run_model_binned:

    for model in ['FLARES_binned', 'FLARES_intrinsic_binned']:

        m = getattr(UV, model)()

        t = Table()

        redshift = np.array([])
        log10L = np.array([])
        phi = np.array([])

        for z in m.redshifts:
            redshift = np.append(redshift, z*np.ones(len(m.log10L))) # only works if only one set of luminosities
            log10L = np.append(np.round(log10L, 2), m.log10L)
            phi = np.append(phi, m.phi[z])

        t.add_column(Column(data = redshift, name = 'redshift', description = 'redshift'))
        t.add_column(Column(data = log10L, name = 'log10L', description = 'log10(luminosity/erg/s/Hz)'))
        t.add_column(Column(data = phi, name = f'phi', description = 'log10(phi/Mpc^-3/dex)'))

        t.meta['name'] = m.name
        t.meta['redshifts'] = m.redshifts
        t.meta['model type'] = m.type
        t.meta['type'] = 'binned'

        t.write(f'models/binned/{m.fname}.ecsv', format = 'ascii.ecsv', overwrite=True)



run_obs_binned = True

if run_obs_binned:

    for model in ['FLARES_binned', 'FLARES_intrinsic_binned']:

        m = getattr(UV, model)()

        t = Table()

        redshift = np.array([])
        log10L = np.array([])
        phi = np.array([])

        for z in m.redshifts:
            redshift = np.append(redshift, z*np.ones(len(m.log10L))) # only works if only one set of luminosities
            log10L = np.append(np.round(log10L, 2), m.log10L)
            phi = np.append(phi, m.phi[z])

        t.add_column(Column(data = redshift, name = 'redshift', description = 'redshift'))
        t.add_column(Column(data = log10L, name = 'log10L', description = 'log10(luminosity/erg/s/Hz)'))
        t.add_column(Column(data = phi, name = f'phi', description = 'log10(phi/Mpc^-3/dex)'))

        t.meta['name'] = m.name
        t.meta['redshifts'] = m.redshifts
        t.meta['model type'] = m.type
        t.meta['type'] = 'binned'

        t.write(f'models/binned/{m.fname}.ecsv', format = 'ascii.ecsv', overwrite=True)



# --- schechter

for model in ['FLARES', 'FLARES_intrinsic', 'Bluetides', 'Ma2019', 'Mason2015', 'Yung2018','TNG_A','TNG_B','TNG_C']:

    m = getattr(UV, model)()
    t = Table()

    t.add_column(Column(data = m.redshifts, name = 'redshift', description = 'redshift'))
    t.add_column(Column(data = m.M_star, name = 'M*', description = 'characteristic absolute magnitude'))
    t.add_column(Column(data = m.log10phi_star, name = 'log10phi*', description = 'log10(number density/Mpc^-3)'))
    t.add_column(Column(data = m.alpha, name = 'alpha', description = 'faint end slope'))

    t.meta['name'] = m.name
    t.meta['model type'] = m.type
    t.meta['type'] = 'Schechter'

    try:
        t.meta['ref'] = m.ref
    except:
        print('no reference')

    try:
        t.meta['ads'] = m.ads
    except:
        print('no ADS')

    try:
        t.meta['arxiv'] = m.arxiv
    except:
        print('no arxiv')

    t.write(f'models/schechter/{m.fname}.ecsv', format = 'ascii.ecsv', overwrite=True)



# t.add_column(Column(data = m.log10L, name = 'log10L', description = 'log10(luminosity/erg/s/Hz)'))
# print(t['log10L'])
#
# for z in m.redshifts:
#     t.add_column(m.phi[z], name = f'phi_{z}')



# q = 10 * u.erg / u.s / u.Hz
# print(np.log10(q))
#
# q = m.log10L * u.dex(u.erg / u.s / u.Hz)
#
# print(q)
# phi = 1.0 * u.Mpc**-3 / u.dex
#
# print(phi)
#
# print(phi.to('Mpc^-3 / mag'))


# t.add_column(m.log10L, name = 'log10L')
