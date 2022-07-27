

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as units

geo = 4.*np.pi*(100.*10.*3.0867*10**16)**2 # factor relating the L to M in cm^2
log10geo = np.log10(geo)



def mlabel(l):
    l_ = l.replace(' ', '\ ')
    return rf'$\rm {l}$'


def fnu_to_m(fnu):

    return -2.5*np.log10(fnu/1E9) + 8.9 # -- assumes flux in nJy

def m_to_fnu(m):

    return 1E9 * 10**(-0.4*(m - 8.9)) # -- flux returned nJy

def fnu_to_Lnu(fnu, cosmo, z):

    """ convert flux to luminosity including the band stretching effect """

    return fnu*(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)/(1E9 * 1E23 * (1.+z))

def Lnu_to_fnu(Lnu, cosmo, z):

    """ convert luminosity to flux including the band stretching effect """

    return 1E9 * 1E23 * Lnu * (1.+ z)/(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)


def Lnu_to_M(Lnu):

    return -2.5*np.log10(Lnu/geo)-48.6

def M_to_Lnu(M):

    return 10**(-0.4*(M+48.6)) * geo


def log10Lnu_to_M(log10Lnu):

    return -2.5*log10Lnu-log10geo-48.6

def M_to_log10Lnu(M):

    return -0.4*(M+48.6) + log10geo


# def DM(cosmo, z):
#     luminosity_distance = cosmo.luminosity_distance(z).to('pc').value
#     return 5*np.log10(luminosity_distance/(np.sqrt(1.+z)*10.))

# def M_to_m(M, cosmo, z):
#     return M + DM(z, cosmo = cosmo)

# def m_to_M(m, cosmo, z):
#     return m - DM(z, cosmo = cosmo)



def bin_centres(bin_edges):
    return 0.5*(bin_edges[:-1]+bin_edges[1:])


def simple_fig(fig_size = 3.5):

    if type(fig_size) == float or type(fig_size) == int:
        fig_size = (fig_size, fig_size)

    fig = plt.figure(figsize = fig_size)

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax


labels = {}
labels['Mstar'] = r'M_{\star}'
labels['LUV'] = r'L_{FUV}'
labels['SFR'] = 'SFR'

tunits = {}
tunits['Mstar'] = r'M_{\odot}'
tunits['LUV'] = r'erg\ s^{-1}\ Hz^{-1}'
tunits['SFR'] = r'M_{\odot}\ yr^{-1}'

def label(x, x_unit):
    if isinstance(x_unit, units.DexUnit):
        u = rf'{x_unit.physical_unit:latex}'
        x_ = x[5:] # remove log10
        x_ = x
        if x_ in labels.keys(): x_ = labels[x_]

        l = rf'$\rm log_{{10}}({x_}/{u.replace("$", "")})$'

    else:
        if x_unit:
            u = rf'/{x_unit:latex}'
        else:
            u = ''
        if x in labels.keys(): x = labels[x]
        l = rf'$\rm {x}{u.replace("$", "")})$'

    print(x, x_unit, l)
    return l



def label_(x):
    return rf'$\rm log_{{10}}({labels[x]}/{tunits[x]})$'
