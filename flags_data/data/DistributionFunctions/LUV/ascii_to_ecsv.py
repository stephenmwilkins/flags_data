

import numpy as np
from astropy.table import Table, Column
import astropy.units as u
from flare.LF.literature import UV


models = []

#(name, redshifts, filename pattern, dex^-1, log10)
# models.append(('liu_2016', [5,6,7,8,9,10], lambda z: f'lf_z{z}.dat', False, True))
# models.append(('ocvirk_2016', [6,7,8,10], lambda z: f'lf_z{z}.dat', False, True))
models.append(('yung_2019', [4,5,6,7,8,10,11,12,13,14,15], lambda z: f'UVLF_z{z}.dat', True, False))


for model in models:

    mname, redshifts, f, dex, lg = model

    fname = mname.replace('_','')

    t = Table()

    redshift = np.array([])
    M = np.array([])
    # log10L = np.array([])
    log10phi = np.array([])

    for z in redshifts:

        data = np.loadtxt(f'original_data/{mname}/{f(z)}').T

        print(data)

        redshift = np.append(redshift, z*np.ones(len(data[0]))) # only works if only one set of luminosities
        # log10L = np.append(np.round(log10L, 2), m.log10L)
        M = np.append(M, data[0])

        if lg:
            log10phi = np.append(log10phi, data[1])
        else:
            log10phi = np.append(log10phi, np.log10(data[1]))

    # if number densities are expressed in dex^-1, convert to mag^-1
    if dex: log10phi += np.log10(0.4)


    t.add_column(Column(data = redshift, name = 'redshift', description = 'redshift'))
    t.add_column(Column(data = M, name = 'M', description = 'absolute magnitude'))
    # t.add_column(Column(data = log10L, name = 'log10L', description = 'log10(luminosity/erg/s/Hz)'))
    t.add_column(Column(data = log10phi, name = f'log10phi', description = 'log10(phi/Mpc^-3/mag)'))

    t.meta['name'] = fname
    t.meta['redshifts'] = redshifts
    t.meta['type'] = 'binned'

    t.write(f'models/binned/{fname}.ecsv', format = 'ascii.ecsv', overwrite=True)
