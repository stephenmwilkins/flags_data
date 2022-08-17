# Generates ecsv files for FLAGS dataset from the Caesar files generated from the SIMBA simulations.
# Usage:
#    % python flags_scalings.py WIND TYPE MODELS
# WIND = s50 for the Simba main run
# TYPE = GSMF for stellar mass fcn, SFR for SFR fcn, LUV for UV1500 lum fcn, LUV_nodust for raw UVLF
# MODELS = m100n1024 m50n1024  to combine the 100 Mpc/h and 50 Mpc/h 2x1024^3 runs.
# Caesar files for m100n1024 are available at simba.roe.ac.uk.  m50n1024 available upon request.

import pylab as plt
import numpy as np
import sys
import os
import caesar
import plotmedian as pm
from yt.units import km, s, Mpc

WIND = sys.argv[1]
TYPE = sys.argv[2]
MODEL = sys.argv[3:]

SNAPNUM = [22,26,30,36,42]

TYPES = [TYPE]
LOGLIMITS = {'HI':-3, 'H2':-4, 'dust':-6, 'BH':-5, 'SFR':-11, 'Zgas':-2, 'Zstar':-2}
colors = ('b', 'crimson', 'g', 'c', 'm', 'y', 'k')

def massScaling(sims):
    mstarall = np.zeros(0)
    fracall = np.zeros(0)
    for sim in sims:
        mstar = np.array([g.masses['stellar'].d for g in sim.galaxies])
        if TYPE == 'HI': mass = np.array([g.masses['HI'].d for g in sim.galaxies]) 
        elif TYPE == 'H2': mass = np.array([g.masses['H2'].d for g in sim.galaxies]) 
        elif TYPE == 'dust': mass = np.array([g.masses['dust'].d for g in sim.galaxies]) 
        elif TYPE == 'BH': mass = np.array([g.masses['bh'].d for g in sim.galaxies]) 
        elif TYPE == 'SFR': mass = np.array([g.sfr.d for g in sim.galaxies]) 
        elif TYPE == 'Zgas': mass = np.array([g.metallicities['sfr_weighted'].d for g in sim.galaxies]) 
        elif TYPE == 'Zstar': mass = np.array([g.metallicities['stellar'].d for g in sim.galaxies]) 
        elif TYPE == 'LUV': mass = np.array([10**g.absmag['i1500'] for g in sim.galaxies]) 
        elif TYPE == 'LUV_nodust': mass = np.array([10**g.absmag['i1500_nodust'] for g in sim.galaxies]) 
        elif TYPE == 'Halo': mass = np.array([h.masses['virial'].d for h in sim.halos]) 
        else: sys.exit('Type %s not recognized'%TYPE)

        frac = mass/mstar
        if TYPE == 'Zgas' or TYPE == 'Zstar': frac = mass/0.0134
        frac = np.clip(frac,10**LOGLIMITS[TYPE],None)

        H0 = (100*sim.simulation.hubble_constant) * km / s / Mpc
        rhocrit = 3.*H0.to('1/s')**2 / (8*np.pi*sim.simulation.G)
        mlim = 32*rhocrit.to('Msun/kpc**3')*sim.simulation.boxsize.to('kpc')**3*sim.simulation.omega_baryon/sim.simulation.effective_resolution**3/sim.simulation.scale_factor**3 # galaxy mass resolution limit
        mstarall = np.concatenate((mstarall,mstar[mstar>mlim]))
        fracall = np.concatenate((fracall,frac[mstar>mlim]))

    nbins = min(12,int(len(mstarall)/4))
    massbin,fracbin,ebinlo,ebinhi = pm.runningmedian(np.log10(mstarall),np.log10(fracall),erange=[15.8,84],bins=nbins)
    #print(nbins, massbin, fracbin, ebinlo, ebinhi)
    return massbin, fracbin, ebinlo, ebinhi

def write_ecsv(z, xval, yval, errlo, errhi):
    if TYPE == 'SFR': 
        ylabel = 'log10sSFR'
        yfields = 'log10sSFR_P15.8 log10sSFR_P50.0 log10sSFR_P84.2'
    elif TYPE == 'HI': 
        ylabel = 'log10fHI'
        yfields = 'log10fHI_P15.8 log10fHI_P50.0 log10fHI_P84.2'
    elif TYPE == 'H2': 
        ylabel = 'log10fH2'
        yfields = 'log10fH2_P15.8 log10fH2_P50.0 log10fH2_P84.2'
    elif TYPE == 'dust': 
        ylabel = 'log10fdust'
        yfields = 'log10fdust_P15.8 log10fdust_P50.0 log10fdust_P84.2'
    elif TYPE == 'BH': ylabel = 'log10fBH'
    elif TYPE == 'Zgas': 
        ylabel = 'log10Zgas'
        yfields = 'log10Zgas_P15.8 log10Zgas_P50.0 log10Zgas_P84.2'
    elif TYPE == 'Zstar': 
        ylabel = 'log10Zstar'
        yfields = 'log10Zstar_P15.8 log10Zstar_P50.0 log10Zstar_P84.2'
    elif 'LUV' in TYPE: 
        ylabel = 'log10LUV'
        yfields = 'log10LUV_P15.8 log10LUV_P50.0 log10LUV_P84.2'
    else: sys.exit('Type %s not recognised'%TYPE)

    # if file doesn't exist, write the header info
    if not os.path.isfile(fname): 
        with open(fname, "w") as f:
            f.write('# %s 1.0\
                     \n# ---\
                     \n# datatype:\
                     \n# - {name: z, datatype: float64}\
                     \n# - {name: log10Mstar, unit: dex(solMass), datatype: float64}\
                     \n# - {name: log10sSFR_P15.8, unit: dex(solMass / yr), datatype: float64}\
                     \n# - {name: log10sSFR_P50.0, unit: dex(solMass / yr), datatype: float64}\
                     \n# - {name: log10sSFR_P84.2, unit: dex(solMass / yr), datatype: float64}\
                     \n# meta: !!omap\
                     \n# - {name: Simba}\
                     \n# - {dataset_type: binned}\
                     \n# - redshifts: [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]\
                     \n# - references: [2021MNRAS.500.2127L]\
                     \n# - {x: log10Mstar}\
                     \n# - {y: %s}\
                     \n# schema: astropy-2.0\
                     \nz log10Mstar %s\n' % ('%ECSV', ylabel, yfields) )

    # write out the data
    with open(fname, "a") as f:
        for i in range(len(xval)):
            if not np.isnan(yval[i]): f.write('%.2f %.2f %.4f %.4f %.4f\n'%(z, xval[i], yval[i]-errlo[i], yval[i], yval[i]+errlo[i]))


###########################################################################

if __name__ == '__main__':

    fig,ax = plt.subplots()

    fname = "simba_scaling_%s.ecsv"%TYPE
    if os.path.isfile(fname): os.remove(fname)  # remove the target file

    minmass = 1.e20
    for k in range(0,len(SNAPNUM)):
        # loop over models to find max phi at each bin location
        sims = []
        for j in range(0,len(MODEL)):
            # load caesar file
            caesarfile = '/home/rad/data/%s/%s/Groups/%s_%03d.hdf5' % (MODEL[j],WIND,MODEL[j],SNAPNUM[k])
            if os.path.isfile(caesarfile):
                sims.append(caesar.load(caesarfile))
            else:
                print('Could not find caesar file %s'%caesarfile)
                continue
        massbin,fracbin,ebinlo,ebinhi = massScaling(sims)
        if massbin[0] < minmass: minmass = massbin[0]
        redshift = np.round(sims[0].simulation.redshift,0)
        if k==len(SNAPNUM)-1: ax.errorbar(massbin, fracbin, yerr=[ebinlo,ebinhi], color=colors[k], label='z=%g'%redshift)
        else: ax.plot(massbin, fracbin, color=colors[k], label='z=%g'%redshift)
        write_ecsv(redshift, massbin, fracbin, ebinlo, ebinhi)

    ax.legend(loc='best',fontsize=12)
    plt.xlim(minmass-0.2,)
    plt.xlabel(r'$\log M_* [M_\odot]$',fontsize=16)
    plt.ylabel(r'$\log$ Quantity/$M_*$',fontsize=16)
    plt.minorticks_on()
    plt.subplots_adjust(hspace=.0)

    figname = 'simba_scaling_%s.pdf' % (TYPE)
    plt.savefig(figname,bbox_inches='tight')
    plt.show()

