# Generates ecsv files for FLAGS dataset from the Caesar files generated from the SIMBA simulations.
# Usage:
#    % python flags.py WIND TYPE MODELS
# WIND = s50 for the Simba main run
# TYPE = GSMF for stellar mass fcn, SFR for SFR fcn, LUV for UV1500 lum fcn, LUV_nodust for raw UVLF
# MODELS = m100n1024 m50n1024  to combine the 100 Mpc/h and 50 Mpc/h 2x1024^3 runs.
# Combination is done by computing LFs of each and then taking highest value among models in each bin.
# Caesar files for m100n1024 are available at simba.roe.ac.uk.  m50n1024 available upon request.

import pylab as plt
import numpy as np
import sys
import os
import caesar
import function as fu
from yt.units import km, s, Mpc

WIND = sys.argv[1]
TYPE = sys.argv[2]
MODEL = sys.argv[3:]

SNAPNUM = [22,26,30,36,42]

TYPES = [TYPE]
LOGZERO = -10
colors = ('b', 'crimson', 'g', 'c', 'm', 'y', 'k')

def massFunc(sim,jwind):
    for curType in TYPES:
        galpos = np.array([g.pos.d for g in sim.galaxies])
        if curType == 'GSMF': 
            mass = np.array([g.masses['stellar'].d for g in sim.galaxies])
        elif curType == 'HI': mass = np.array([g.masses['HI'].d for g in sim.galaxies]) 
        elif curType == 'H2': mass = np.array([g.masses['H2'].d for g in sim.galaxies]) 
        elif curType == 'dust': mass = np.array([g.masses['dust'].d for g in sim.galaxies]) 
        elif curType == 'SFR': mass = np.array([g.sfr.d for g in sim.galaxies]) 
        elif curType == 'LUV': mass = np.array([g.absmag['i1500'] for g in sim.galaxies]) 
        elif curType == 'LUV_nodust': mass = np.array([g.absmag_nodust['i1500'] for g in sim.galaxies]) 
        elif curType == 'Halo': 
            mass = np.array([h.masses['virial'].d for h in sim.halos]) 
            galpos = np.array([h.pos.d for h in sim.halos]) 
        volume = sim.simulation.boxsize.to('Mpc').d**3 * (1+sim.simulation.redshift)**3
        if curType == 'SFR': 
            galpos = galpos[mass>1.e-2]
            mass = mass[mass>1.e-2]
            mlim = np.log10(min(mass))
        elif 'LUV' in curType: 
            mass = 10**(-mass)  # NOTE: switch sign for mags to make histogramming work
            mlim = 15
        else:
            H0 = (100*sim.simulation.hubble_constant) * km / s / Mpc
            rhocrit = 3.*H0.to('1/s')**2 / (8*np.pi*sim.simulation.G)
            mlim = np.log10(16*rhocrit.to('Msun/kpc**3')*sim.simulation.boxsize.to('kpc')**3*sim.simulation.omega_baryon/sim.simulation.effective_resolution**3/sim.simulation.scale_factor**3) # galaxy mass resolution limit: 32 gas particle masses

        #print(mlim,sim.simulation.redshift,min(np.log10(mass)),max(np.log10(mass)))
        nbin = int( (max(np.log10(mass))-mlim)/BINSIZE ) + 1
        x,y,sig = fu.cosmic_variance(mass, galpos, sim.simulation.boxsize, volume, nbin=nbin, minmass=mlim)
        ngbin, edges = np.histogram(np.log10(mass),bins=nbin,range=(mlim,max(np.log10(mass))))
        print(mlim,np.log10(max(mass)),len(mass),nbin,ngbin)
        return np.log10(x),np.log10(np.clip(y,10**LOGZERO,None)),sig

def write_ecsv(z, xval, yval, errval):
    if TYPE == 'GSMF': xlabel = 'log10Mstar'
    elif TYPE == 'SFR': xlabel = 'log10SFR'
    elif TYPE == 'HI': xlabel = 'log10MHI'
    elif TYPE == 'H2': xlabel = 'log10MH2'
    elif 'LUV' in TYPE: xlabel = 'MUV'
    else: sys.exit('Type %s not recognised'%TYPE)

    # if file doesn't exist, write the header info
    if not os.path.isfile(fname): 
        with open(fname, "w") as f:
            f.write('# %s 1.0\n# ---\n# datatype:\n# - {name: z, datatype: float64, description: redshift}\n# - {name: log10L, datatype: float64, unit: dex(erg s^-1 Hz^-1)}\n# - {name: log10phi, datatype: float64, unit: dex(Mpc^-3 dex^-1)}\n# - {name: err_log10phi, datatype: float64, unit: dex(Mpc^-3 dex^-1)}\n# meta: !!omap\n# - {name: Simba}\n# - {x: %s}\n# - {y: log10phi}\n# - references: [2019MNRAS.486.2827D]\n# - redshifts: [5.0, 6.0, 7.0, 8.0, 9.0]\n# - {model type: hydro}\n# - {type: binned}\n# schema: astropy-2.0\nz %s log10phi err_log10phi\n'%(('%ECSV',xlabel,xlabel)))

    # write out the data
    with open(fname, "a") as f:
        for i in range(len(xval)):
            if 'LUV' in TYPE: xval[i] = -xval[i]   # switch sign back for mags
            f.write('%.2f %.2f %.4f %.4f\n'%(z, xval[i], yval[i], errval[i]))


###########################################################################

if __name__ == '__main__':

    fig,ax = plt.subplots()

    if TYPE == 'SFR':
        BINSIZE = 0.2
        mbin = np.arange(-2,3,BINSIZE)  # these are our final desired bins
    elif 'LUV' in TYPE:
        BINSIZE = 0.2
        mbin = np.arange(15,24,BINSIZE)  # these are our final desired bins
    else:
        BINSIZE = 0.2
        mbin = np.arange(7.8,12,BINSIZE)  # these are our final desired bins
    fname = "simba_%s.ecsv"%TYPE
    if os.path.isfile(fname): os.remove(fname)  # remove the target file

    minmass = 1.e20
    for k in range(0,len(SNAPNUM)):
        phi = np.zeros(len(mbin))+LOGZERO  # initialize log phi to small number
        ephi = np.zeros(len(mbin))+LOGZERO  # initialize dlogphi to small number
        # loop over models to find max phi at each bin location
        for j in range(0,len(MODEL)):
            # load caesar file
            caesarfile = '/home/rad/data/%s/%s/Groups/%s_%03d.hdf5' % (MODEL[j],WIND,MODEL[j],SNAPNUM[k])
            if os.path.isfile(caesarfile):
                sim = caesar.load(caesarfile)
            else:
                print('Could not find caesar file %s'%caesarfile)
                continue
            # compute mass function
            mb,ph,eph = massFunc(sim,j)
            if mb[0] < minmass: minmass = mb[0] # set lower limit
            # add "zero" values to end of mass function
            mb = np.insert(mb, 0, mb[0]-0.001, axis=0)
            ph = np.insert(ph, 0, LOGZERO, axis=0)
            eph = np.insert(eph, 0, LOGZERO, axis=0)
            mb = np.append(mb, mb[-1]+0.001)
            ph = np.append(ph, LOGZERO)
            eph = np.append(eph, LOGZERO)
            # (log)linearly interpolate mass function to desired mass bins
            phinterp = np.interp(mbin, mb, ph)
            ephinterp = np.interp(mbin, mb, eph)
            # store max value among models for each bin
            phi = np.where( phi>phinterp, phi, phinterp)  
            ephi = np.where( ephi>ephinterp, ephi, ephinterp)  
        # ignore low bins for plotting
        mbin_plot = mbin[phi>-6]
        ephi_plot = ephi[phi>-6]
        phi_plot = phi[phi>-6]
        redshift = np.round(sim.simulation.redshift,0)
        #ax.plot(mbin_plot, phi_plot, color=colors[k], label='z=%g'%redshift)
        ax.errorbar(mbin_plot, phi_plot, yerr=ephi_plot, color=colors[k], label='z=%g'%redshift)
        write_ecsv(redshift, mbin_plot, phi_plot, ephi_plot)
        #print(mbin)
        #print(phi)


    ax.legend(loc='best',fontsize=12)
    plt.xlim(minmass,)
    plt.xlabel(r'$\log M_{%s} [M_\odot]$'%TYPE,fontsize=16)
    plt.ylabel(r'$\log \Phi [Mpc^{-3}$]',fontsize=16)
    plt.minorticks_on()
    plt.subplots_adjust(hspace=.0)

    figname = 'df_%s_%s.pdf' % (TYPE,WIND)
    plt.savefig(figname,bbox_inches='tight')
    plt.show()

