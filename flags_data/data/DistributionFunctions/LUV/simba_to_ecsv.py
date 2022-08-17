
import os
import h5py
import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u




# 0.010000  99.000000 000
#   0.048194  19.749544 001
#   0.050227  18.909730 002
#   0.052301  18.119968 003
#   0.054418  17.376352 004
#   0.056576  16.675347 005
#   0.058776  16.013753 006
#   0.061018  15.388664 007
#   0.063301  14.797437 008
#   0.065627  14.237665 009
#   0.067994  13.707152 010
#   0.070403  13.203888 011
#   0.072854  12.726037 012
#   0.075347  12.271910 013
#   0.077882  11.839961 014
#   0.080459  11.428765 015
#   0.083077  11.037011 016
#   0.085738  10.663489 017
#   0.088440  10.307082 018
#   0.091185   9.966759 019
#   0.093971   9.641562 020
#   0.096800   9.330607 021
#   0.099670   9.033071 022
#   0.102583   8.748191 023
#   0.105538   8.475257 024
#   0.108535   8.213608 025
#   0.111574   7.962628 026
#   0.114656   7.721742 027
#   0.117780   7.490415 028
#   0.120946   7.268145 029
#   0.124155   7.054463 030
#   0.127406   6.848931 031
#   0.130700   6.651136 032
#   0.134036   6.460693 033
#   0.137415   6.277241 034
#   0.140836   6.100439 035
#   0.144301   5.929968 036
#   0.147808   5.765526 037
#   0.151359   5.606830 038
#   0.154952   5.453615 039
#   0.158589   5.305627 040


snap_z = {'006':16.013753, '008':14.797437,  '011': 13.203888,  '014':11.839961, '016': 11.037011,  '019': 9.966759, '022': 9.033071, '026': 7.962628, '030': 7.054463, '036': 5.929968, '037': 5.765526, '038': 5.606830, '039':5.453615, '040':5.305627, '041':5.162630, '042': 5.024400, '043':4.890725, '044':4.761405, '045':4.636251, '046': 4.515082}

tables = []

volume = (100)**3

for snap in ['014', '016', '019', '022', '026', '030', '036', '042']:

    print(snap)
    z = snap_z[snap]

    d = h5py.File(f'/Users/stephenwilkins/Dropbox/Research/data/simulations/simba/m100n1024_{snap}.hdf5','r')

    log10Mstar = np.log10(np.array(d[F'galaxy_data/dicts'].get('masses.stellar_30kpc'), dtype = np.float64))
    # ok = np.where(mstar>1e8)[0]

    print(np.min(log10Mstar), np.median(log10Mstar), np.max(log10Mstar))

    Ms = np.array(d[F'galaxy_data/dicts'].get(f'absmag.i1500'), dtype = np.float64)

    print(np.min(Ms), np.median(Ms), np.max(Ms))

    bin_w = 0.2

    bin_edges = np.arange(-22., -18., bin_w)
    M = 0.5*(bin_edges[:-1]+bin_edges[1:])

    N, _ = np.histogram(Ms, bins = bin_edges)
    print(N)


    phi = N/volume

    s = N>0

    t = Table()
    t.add_column(Column(data = z*np.ones(np.sum(s)), name = 'z'))
    t.add_column(Column(data = M[s], name = 'M', unit = 'mag'))
    t.add_column(Column(data = np.log10(phi[s]/bin_w), name = 'log10phi', unit = 'dex(Mpc^-3 mag^-1)'))
    tables.append(t)

table = vstack(tables)

table.meta['x'] = 'M'
table.meta['y'] = 'log10phi'
table.meta['name'] = 'Simba'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = []

out_name = 'simba'

table.write(f'models/binned/{out_name}.ecsv', format = 'ascii.ecsv', overwrite=True)
