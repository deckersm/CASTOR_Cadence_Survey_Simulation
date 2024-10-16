import pandas as pd
import glob
import sys
import os
import time

# importing functions from other files
import utils_castor as utils
import rates_castor as rates
import simulation_functions_castor as simul
import statistics_castor as stats


# To run this file, run the following command from the terminal:
# python main.py <transient type> <maximum redshift> <survey cadence>

# This command will extract the volumetric rate based on the transient shape and inject transients out to the maximum redshift, with redshifts sampled according to the volumetric rate
# It will next choose a random ra and dec from within a specified field and apply Milky Way extinction

# For each redshift in the array it chooses one of the models of the transient type, where the range of models is already representative of the luminosity function

# Then it runs all the spectral templates of that particular model through FORECASTOR to obtain the SNR and mag for each data point in the light curve
# Then it runs detection statistics to determine which light curves would classify as detected and at what phase/mag they are detected

# Requires template filenames to be of the format:
# path_to_template_folder/{type}/SED_{type}_{model}_{phase}d.dat

# Also requires the following folders to be set up within the ETC folder:
# results/  -> which will store the output light curves for each model at each redshift as seen by CASTOR (filename: results....) 
# and a statistics file summarising detection statistics on each model (filename: statistics......)

######################################################################################################################################################################################

if __name__ == '__main__':

    # Count how long simulation takes to run
    start = time.time()

    # Extracts the transient type, maximum redshift, and survey cadence from the terminal input
    type = sys.argv[1]
    max_z = float(sys.argv[2])

    # If no survey cadence is specified, it is assumed to be 1.0 d
    if len(sys.argv) > 3:
        cadence = sys.argv[3]
    else:
        cadence = 1.0

    # Globally defining background and telescope as these won't change between transients
    MyTelescope = utils.config_telescope();
    MyBackground = utils.config_background();

    # Extracting all available models from the spectral templates folder
    models = []
    files = glob.glob('/users/deckersm/CASTOR/Templates/individual_spectral_templates/{}/SED_{}_*d.dat'.format(type, type))
    for f in files:
        models.append(f.split('/')[-1].split('_')[2])
    models = list(set(models))
    
    # Populating the redshift range with the transients, extracting the light curves as they would be seen by CASTOR
    all_results = simul.populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, cadence = cadence)

    # Running detection statistics in the three CASTOR filters
    print('Finished simulating light curves, now running statistics \n')
    for band in ['uv', 'u', 'g']:
        overview = stats.statistics(all_results, max_z, type, band = band, cadence = cadence)
        overview.to_csv('results/statistics_{}_{}_{}_{}d.csv'.format(type, max_z, band, cadence), index = False)

    end = time.time()
    print('Total run time for simulation = {} seconds'.format(end - start))

