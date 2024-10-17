import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.colors import LogNorm
import pandas as pd
import glob
import os

from castor_etc.background import Background
from castor_etc.photometry import Photometry
from castor_etc.sources import ExtendedSource, GalaxySource, PointSource
from castor_etc.telescope import Telescope

import spectral_interpolation_castor as spec_interp

######################################################################################################################################################################################
#### #### #### #### ####                         Functions to set up the telescope, background, and source objects.                                                #### #### #### #### 
######################################################################################################################################################################################


# Adopts default CASTOR telescope parameters
def config_telescope():
    MyTelescope = Telescope()
    return MyTelescope


# Adopts default background parameters defined in FORECASTOR
def config_background():
    MyBackground = Background()
    return MyBackground


# Configure source at particular redshift from model SED at particular phase
def config_source(type, model, phase, z):

    # Create point source
    MySource = PointSource()

    # Import SED, checking first if file exist, otherwise creating spectrum using spectral interpolation functions
    if os.path.isfile('/users/deckersm/CASTOR/Templates/individual_spectral_templates/{}/SED_{}_{}_{}d.dat'.format(type, type, model, phase)) == True:
        filename = '/users/deckersm/CASTOR/Templates/individual_spectral_templates/{}/SED_{}_{}_{}d.dat'.format(type, type, model, phase)
    elif os.path.isfile('/users/deckersm/CASTOR/Templates/individual_spectral_templates/interpolated_spectra/SED_{}_{}_{}d.dat'.format(type, type, model, phase)) == True:
        filename = '/users/deckersm/CASTOR/Templates/individual_spectral_templates/interpolated_spectra/SED_{}_{}_{}d.dat'.format(type, type, model, phase)
    else:
        spectrum = spec_interp.create_spec_at_phase(type, model, phase)
        spectrum.to_csv('/users/deckersm/CASTOR/Templates/individual_spectral_templates/interpolated_spectra/SED_{}_{}_{}d.dat'.format(type, type, model, phase), index = False, encoding='utf-8', sep = ' ')
        filename = '/users/deckersm/CASTOR/Templates/individual_spectral_templates/interpolated_spectra/SED_{}_{}_{}d.dat'.format(type, type, model, phase)
    
    MySource.use_custom_spectrum(filename)

    # Putting the spectrum at distance d (correcting flux and wavelength)
    c = 3e5 # km/s
    H0 = 70 # km/sec/Mpc
    d = c*z / H0
    d_cm = d * 3.086e+24 # converting Mpc to cm
    d_pc = d * 1e6 # converting Mpc to pc

    #if type == 'snia' and model == 'hsiao':
    #   MySource.spectrum = MySource.spectrum * 5e47 # Scaling the Hsiao spectrum to have a reasonable peak absolute magnitue
    
    # Correcting the flux for redshift 
    #MySource.spectrum = MySource.spectrum / (4 * math.pi * (d_pc**2)) / (1. + z)

    MySource.spectrum = MySource.spectrum * ((10**2)/(d_pc**2)) / (1. + z)

    # Correcting wavelength for redshift
    MySource.wavelengths = MySource.wavelengths * (1. + z)

    return MySource


######################################################################################################################################################################################
#### #### #### #### ####                         Functions to compute photometry with ETC and produce light curve.                                                 #### #### #### #### 
######################################################################################################################################################################################


# Perform photometry on source in CASTOR filters and return apparent mag and SNR from given exposure time
# Vectorized get_photometry
def get_photometry(MyTelescope, MyBackground, MySource, ebv=0.01, time=100):
    
    MyPhot = Photometry(MyTelescope, MySource, MyBackground)
    MyPhot.use_optimal_aperture(quiet=True)
    
    mags = np.array([MySource.get_AB_mag(MyTelescope)[band] for band in MyTelescope.passbands])
    snrs = np.array([MyPhot.calc_snr_or_t(t=time, reddening=ebv, quiet=True)[band] for band in MyTelescope.passbands])
    
    return MyTelescope.passbands, mags, snrs

# Function to plot the output light curve in apparent mag
def visualise_lc(results):
    fig = plt.figure(dpi=150, figsize = (8,6), facecolor = 'w')
    ax = fig.add_subplot(111)
    
    colours = ['plum', 'slateblue', 'mediumaquamarine']
    c=0
    for band in MyTelescope.passbands:
        ax.errorbar(results.loc[results['filter']==band, 'phase'], results.loc[results['filter']==band, 'mag'], yerr = results.loc[results['filter']==band, 'mag_err'], fmt = '.', linestyle = None, label = band, color = colours[c])
        c+=1

        
    ax.invert_yaxis()
    ax.set_xlabel('Time since explosion [d]')
    ax.set_ylabel('Apparent magnitude')
    ax.legend()
    plt.show()
    plt.close()


# Loops over all available SEDs for a particular model to produce light curve as seen through CASTOR at a given redshift
def create_lc(type, model, z, ra, dec, ebv, number, MyTelescope, MyBackground, cadence = 1.0, start_time = 0, max_phase = 50):
    results = pd.DataFrame(columns = ['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'time', 'phase', 'filter', 'mag', 'mag_err', 'snr'])
    #MyBackground = config_background()
    #MyTelescope = config_telescope()

    files = glob.glob('/users/deckersm/CASTOR/Templates/individual_spectral_templates/{}/SED_{}_{}_*d.dat'.format(type, type, model))
    phases_ = []
    for f in files:
        #phases_.append(float(f.split('/')[-1].split('_')[3].replace('d.dat', '')))
        phases_.append(float(f.split('/')[-1].split('_')[3].replace('d.dat', '')))

    phases = np.arange(float(np.nanmin(phases_)), float(np.nanmax(phases_)), float(cadence))

    for phase in phases:
        if phase < max_phase:
            MySource = config_source(type, model, phase, z)
            filters, mags, snrs = get_photometry(MyTelescope, MyBackground, MySource, ebv)
            
            for i in range(len(filters)):
                if snrs[i]!=0:
                    results.loc[len(results), ['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'time', 'phase', 'filter', 'mag', 'mag_err', 'snr']] = number, type, model, z, ra, dec, ebv, phase + start_time, phase, filters[i], mags[i], 1./snrs[i], snrs[i]
                else:
                    pass 
    return results





