
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
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor
import time
import psutil
from astropy.coordinates import SkyCoord
import dustmaps.sfd
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import freeze_support
import multiprocessing as mp
from functools import partial
import gc

# importing functions from other files
import utils_castor as utils
import rates_castor as rates
import simulation_functions_castor as simul
import statistics_castor as stats


from castor_etc.background import Background
from castor_etc.photometry import Photometry
from castor_etc.sources import ExtendedSource, GalaxySource, PointSource
from castor_etc.telescope import Telescope


# To run this file, run the following command from the terminal:
# python redshift_statistics.py <transient type> <minimum redshift> <maximum redshift>

# Executing this file will run through all available models for the transient type
# And put them at the range of redshifts specified
# Run them through FORECASTOR to obtain the SNR and mag for each data point in the light curves
# Then run detection statistics to determine which light curves would classify as detected and at what phase/mag they are detected

# Requires template filenames to be of the format:
# path_to_template_folder/{type}/SED_{type}_{model}_{phase}d.dat

# Also requires the following folders to be set up within the ETC folder:
# redshift results/  -> which will store the output light curves for each model at each redshift as seen by CASTOR
# statistics/   -> which will store the detection statistics for each model at the range of redshifts specified by the user

######################################################################################################################################################################################
#### #### #### #### ####                         Lines below will simulate survey using the functions defined above.                                               #### #### #### #### 
######################################################################################################################################################################################




if __name__ == '__main__':

    # Count how long simulation takes to run
    start = time.time()

    type = sys.argv[1]
    max_z = float(sys.argv[2])
    band = sys.argv[3]

    if len(sys.argv) > 4:
        cadence = sys_argv[4]
    else:
        cadence = 1

    # Globally defining background and telescope as these won't change between transients
    MyTelescope = utils.config_telescope();
    MyBackground = utils.config_background();

    models = []
    files = glob.glob('/Users/maximedeckers/Documents/RSE/CASTOR/Templates/individual_spectral_templates/{}/SED_{}_*d.dat'.format(type, type))
    for f in files:
        models.append(f.split('/')[-1].split('_')[2])

    models = list(set(models))
    
    freeze_support()  # Optional: Needed only if creating a frozen executable
    
    all_results = simul.populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, cadence = cadence)


    print('Finished simulating light curves, now running statistics \n')
    overview = stats.statistics(all_results, max_z, type, band = band)
    overview.to_csv('results/statistics_{}_{}_{}.csv'.format(type, max_z, band), index = False)

    end = time.time()
    print('Total run time for simulation = {} seconds'.format(end - start))


