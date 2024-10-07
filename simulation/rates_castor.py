from astropy.cosmology import Planck18 as cosmo  # Using Planck18 cosmology
from astropy.cosmology import FlatLambdaCDM 
from astropy import units as u
from astropy.coordinates import SkyCoord
from dustmaps.planck import PlanckQuery
from dustmaps.sfd import SFDQuery
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.colors import LogNorm
import pandas as pd
import glob
import sys
import os
import extinction
import time
import psutil
from astropy.coordinates import SkyCoord
import dustmaps.sfd
import gc


from castor_etc.background import Background
from castor_etc.photometry import Photometry
from castor_etc.sources import ExtendedSource, GalaxySource, PointSource
from castor_etc.telescope import Telescope



######################################################################################################################################################################################
#### #### #### #### ####                         Functions to inject tranients at particular redshifts depending on rates                                          #### #### #### #### 
######################################################################################################################################################################################



# Star formation rate density function (Madau & Dickinson 2014) -- check for more up to date versions (2018>) check Graur (2017)
def sfrd_madau_dickinson(z):
    sfrd_0 = 0.015  # M_sun yr^-1 Mpc^-3 at z=0
    return sfrd_0 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)


# Volumetric rate scaling function to ensure our rate is as a function of redshift
def md14_rate_evolution(z, r_sn_0, z0):
    sfrd_z = sfrd_madau_dickinson(z)
    sfrd_0 = sfrd_madau_dickinson(z0)  # SFRD at z=0
    return r_sn_0 * (sfrd_z / sfrd_0)

# SN Ia rate split at z=1.0, parametrisation taken from Hounsell (2018). However, ELASTICC implementation has beta = -0.5 whereas Hounsell (2018) has beta = 0.5 for z > 1.0?? -> adopted from plasticc code
def snia_rate(z):
    #if z <= 1.0:
    #    return 2.5e-5 * (1. + z)**1.5
    #else:
    #    return 9.7e-5 * (1. + z)**(0.5)
    return np.piecewise(z, [z <= 1.0, z > 1.0], [lambda x:2.5e-5 * (1. + x)**1.5, lambda x:9.7e-5 * (1. + x)**(0.5)])
   

"""
# Function returning the CC volumetric rate at redshift intervals taken from Strolger (2015) -> adopted from plasticc code
def cc_rate(z):

    fx = np.piecewise(z, [(z <= 0.08) * (z > 0), 
                            (z <= 0.29) * (z > 0.08), 
                            (z <= 0.47) * (z > 0.29), 
                            (z <= 0.72) * (z > 0.29), 
                            (z <= 1.56) * (z > 0.72), 
                            (z > 1.56)], 
                            [0.72e-4, 1.33e-4, 1.81e-4, 3.91e-4, 3.22e-4, 3.76e-4])

    return fx
"""

def cc_rate(z):
    strolger = pd.DataFrame(pd.read_csv('../Rates/strolger2015_cc_rate.csv'))
    rates = []
    for i in range(len(z)):
        try:
            index = np.where(strolger['z'] == np.round(z[i], 2))[0][0]
            rates.append(strolger['rate'][index])
            return rates
        except IndexError:
            print(np.round(z[i], 2))


# TDE rate at z = 0 taken from Van Velzen (2017) and evolution with redshift assumed to follow SFRD -> adopted from plasticc code
def tde_rate(z):
    return md14_rate_evolution(z, 1.0e-6, 0)

# SN Iax rate at z = 0 taken from ELASTICC and evolution with redshift assumed to follow SFRD -> adopted from plasticc code
def iax_rate(z):
    return md14_rate_evolution(z, 6.0e-6, 0)

# SN 91bg rate at z = 0 taken from ELASTICC and evolution with redshift following powerlaw -> adopted from plasticc code
def sn91bg_rate(z):
    return 3e-6 * (1+z)**1.5

# SLSN rate at z = 0 taken from ELASTICC and evolution with redshift following SFRD -> adopted from plasticc code
def slsn_rate(z):
    return md14_rate_evolution(z, 2.0e-8, 0)

# SN Ib/c rate assumed to be 30% of CC rate -- following Shivvers (2017) -> adopted from plasticc code
def snibc_rate(z):
    return cc_rate(z) * 0.3

def transient_rate(type, z):
    if type == 'snia':
        return snia_rate(z)
    elif type == 'snii':
        return cc_rate(z)
    elif type == 'snibc':
        return snibc_rate(z)
    elif type == 'slsn':
        return slsn_rate(z)
    elif type == 'tde':
        return tde_rate(z)
    elif type == '91bg':
        return sn91bg_rate(z)
    elif type == 'iax':
        return iax_rate(z)


