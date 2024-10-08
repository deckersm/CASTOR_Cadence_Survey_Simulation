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
# python main.py <transient type> <maximum redshift> <filter for statistics> <survey cadence>

# This command will extract the volumetric rate based on the transient shape and inject transients out to the maximum redshift, with redshifts sampled according to the volumetric rate
# It will next choose a random ra and dec from within a specified field and apply Milky Way extinction

# For each redshift in the array it chooses one of the models of the transient type, where the range of models is already representative of the luminosity function

# Then it runs al the spectral templates of that particular model through FORECASTOR to obtain the SNR and mag for each data point in the light curve
# Then run detection statistics to determine which light curves would classify as detected and at what phase/mag they are detected

# Requires template filenames to be of the format:
# path_to_template_folder/{type}/SED_{type}_{model}_{phase}d.dat

# Also requires the following folders to be set up within the ETC folder:
# results/  -> which will store the output light curves for each model at each redshift as seen by CASTOR (filename: results....) 
# and a statistics file summarising detection statistics on each model (filename: statistics......)

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
    
    all_results = simul.populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, cadence = cadence)

    print('Finished simulating light curves, now running statistics \n')
    overview = stats.statistics(all_results, max_z, type, band = band)
    overview.to_csv('results/statistics_{}_{}_{}.csv'.format(type, max_z, band), index = False)

    end = time.time()
    print('Total run time for simulation = {} seconds'.format(end - start))

