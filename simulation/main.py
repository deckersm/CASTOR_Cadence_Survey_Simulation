import pandas as pd
import glob
import sys
import os
import time
import argparse
import numpy as np

# importing functions from other files
import utils_castor as utils
import rates_castor as rates
import simulation_functions_castor as simul
import statistics_castor as stats


#filepath_to_castor_folder = '/users/deckersm/CASTOR/'
filepath_to_castor_folder = '/Users/maximedeckers/Documents/RSE/CASTOR/CASTOR_Cadence_Survey_Simulation/'


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

    def collect_as(coll_type):
        class Collect_as(argparse.Action):
            def __call__(self, parser, namespace, values, options_string=None):
                setattr(namespace, self.dest, coll_type(values))
        return Collect_as

    # Count how long simulation takes to run
    start = time.time()

    parser = argparse.ArgumentParser(description='Code to simulate CASTOR transient survey')

    parser.add_argument('--type', '-t', help='Type of transient')
    parser.add_argument('--max_redshift', '--max_z', '-z', help='Maximum redshift to simulate transients out to', type=np.float64)
    parser.add_argument('--min_redshift', '--min_z', help='Minimum redshift to simulate transients at', type=np.float64)
    parser.add_argument('--cadence', '-c', help='Cadence of survey', type=np.float64)      
    parser.add_argument('--exposure', '-e', help='Exposure time of each visit of survey', type=np.float64)
    parser.add_argument('--length_survey', '-l', help='Duration of survey in days', type=np.float64)
    parser.add_argument('--field_ra', '-r', help='ras of the centers of the fields for the simulation, separated by spaces', type=np.float64, action=collect_as(np.array))
    parser.add_argument('--field_dec', '-d', help='decs of the centers of the fields for the simulation, separated by spaces', type=np.float64, action=collect_as(np.array))
    parser.add_argument('--simul_type', '-s', help='Simulation type (full = simulate transients based on volumetric rates, test = inject transients are regular redshift intervals to check when detection efficiency drops off')
    parser.add_argument('--number_redshifts', help='Number of redshifts between minz and maxz to inject transients to check detection efficiency', type=int)

    args = parser.parse_args()
    
    if args.type != None:
        type = args.type
    else:
        print("Error: No transient type passed")
    
    if args.max_redshift != None:
        max_z = args.max_redshift
    else:
        print("Error: No maximum redshift passed")
    
    if args.min_redshift != None:
        min_z = args.min_redshift
    else:
        min_z = 0.01
    
    if args.cadence != None:
        cadence = args.cadence
    else:
        cadence = 1.0

    if args.exposure != None:
        exposure = args.exposure
    else:
        exposure = 100.0
    
    if args.length_survey != None:
        length_survey = args.length_survey
    else:
        length_survey = 182.0

    if args.number_redshifts != None:
        number_redshifts = args.number_redshifts
    else:
        number_redshifts = 10

    ra_array = args.field_ra
    dec_array = args.field_dec

    # Globally defining background and telescope as these won't change between transients
    MyTelescope = utils.config_telescope();
    MyBackground = utils.config_background();

    # Extracting all available models from the spectral templates folder
    models = []
    files = glob.glob(filepath_to_castor_folder + 'Templates/{}/SED_{}_*d.dat'.format(type, type))
    for f in files:
        models.append(f.split('/')[-1].split('_')[2])
    models = list(set(models))

    if args.simul_type == 'full':

        # If no coordinates for the fields are provided, then run simulation on the four LSST deep drilling fields
        if ra_array == None:
            centers = {
                "ELAIS_S1": (9.45, -44.0),
                "XMM_LSS": (35.708, -4.75),
                "Extended_Chandra": (53.125, -28.1),
                "COSMOS": (150.1, 2.182)
            }

            for field_name, (ra_center, dec_center) in centers.items():        
                # Populating the redshift range with the transients, extracting the light curves as they would be seen by CASTOR
                all_results = simul.populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, cadence = cadence, exposure = exposure, survey_time = 182, c_ra = ra_center, c_dec = dec_center)

        else:
            for ra_center, dec_center in zip(ra_array, dec_array):
                all_results = simul.populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, cadence = cadence, exposure = exposure, survey_time = 182, c_ra = ra_center, c_dec = dec_center)

        # Running detection statistics in the three CASTOR filters
        print('Finished simulating light curves, now running statistics \n')
        for band in ['uv', 'u', 'g']:
            overview = stats.statistics(all_results, max_z, type, band = band, cadence = cadence, exposure = exposure)
            overview.to_csv('results/statistics_{}_{}_{}_{}d_{}s.csv'.format(type, max_z, band, cadence, exposure, c_ra = ra_center, c_dec = dec_center), index = False)

    elif args.simul_type == 'test':

        print(f'Running quick simulation to test detection efficiencies for {type} transients between {min_z} < z < {max_z}')

        all_results = simul.populate_redshift_range_test(type, models, max_z, MyTelescope, MyBackground, min_z = min_z, number_redshifts = number_redshifts)
        
        for band in ['uv', 'u', 'g']:
            overview = stats.statistics(all_results, max_z, type, band = band, test = True, number_redshifts = number_redshifts)
            overview.to_csv(f'results/statistics_{type}_{max_z}_{band}_{cadence}d_{exposure}s.csv', index = False)



    end = time.time()
    print('Total run time for simulation = {} seconds'.format(end - start))

