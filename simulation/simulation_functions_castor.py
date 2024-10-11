import numpy as np
import math
import pandas as pd
import glob
import os
import extinction
from astropy import units as u
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
from astropy.cosmology import FlatLambdaCDM 


# importing functions from other files
import utils_castor as utils
import rates_castor as rates



def redshift_samples(type = 'snia', z_min = 0.001, z_max = 0.5, survey_time = 1, survey_area = 20):
    
    # Adding units to survey time
    survey_time = survey_time * u.yr

    fraction_sky = survey_area / 41253
    
    # Number of redshift bins
    z_bins = np.linspace(z_min, z_max, 1000)

    # Extracting rate depending on type of transient and adding units
    rate = rates.transient_rate(type, z_bins)
    rate = rate * u.Mpc**(-3) * u.yr**(-1)

    # Calculate the comoving volume element dV/dz (Mpc^3/sr) for each redshift bin
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dVdz = cosmo.differential_comoving_volume(z_bins)  # Mpc^3/sr
    
    # Multiply by 4π (the whole sky in steradians) and the sky fraction to get the comoving volume over the survey area
    dVdz_survey = dVdz * 4 * np.pi * fraction_sky  # Mpc^3 for the fraction of the sky
    
    # Integrate to find the total comoving volume surveyed up to z_max, weighted by the rate
    total_volume = np.trapz(dVdz_survey * rate, z_bins)  # Mpc^3 * yr^-1
    
    # Calculate total number of transients over the survey period (N_total is now a scalar)
    N_total = (total_volume * survey_time).decompose()  # Dimensionless (just a number of transients)
    N_total = N_total.value  # Get rid of any remaining units, making N_total a scalar

    # Probability distribution proportional to dV/dz
    pdf = (dVdz * rate).value
    pdf /= np.sum(pdf)  # Normalize to get a proper probability distribution
    
    # Sample redshifts according to the distribution
    N_samples = int(N_total)  # Total number of transients to sample
    sampled_redshifts = np.random.choice(z_bins, size=N_samples, p=pdf)

    return sampled_redshifts



# Function to pull an array of random ra and dec combination within the field to match the array of redshifts
def ra_dec_samples(c_ra, c_dec, radius, number):

    ra = np.random.uniform(c_ra - radius, c_ra + radius, number)
    dec = np.random.uniform(c_dec - radius, c_dec + radius, number)

    return ra, dec


def time_samples(duration_survey, number, time_before_survey = 30):
    time_array = np.random.uniform(-time_before_survey, duration_survey, number)
    return time_array



# Preload the SFD map at the start and reuse it for all calculations
#dustmaps.sfd.fetch()  # Ensure dust maps are fetched locally
sfd = SFDQuery()  # This loads the SFD map

def mw_ebv(ra, dec):
    coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
    ebv = sfd(coords)  # Reuse the preloaded sfd object
    return ebv
    


# Parallelized function for simulating light curves at a given redshift
def process_single_redshift(i, type, models, redshift_array, ra_array, dec_array, time_array, cadence, MyTelescope, MyBackground):
    z = redshift_array[i]
    ra, dec = ra_array[i], dec_array[i]
    time = time_array[i]

    ebv = mw_ebv(ra, dec)  # Calculate extinction at this RA and Dec
    model = models[int(np.random.uniform(0, len(models), 1))]
    results = utils.create_lc(type, model, z, ra, dec, ebv, i, MyTelescope, MyBackground, start_time = time, cadence = cadence)
    return results


def populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, c_ra=9.45, c_dec=-44.0, radius=1.75, num_cores=8, cadence = 1):
    results_filename = f'results/results_{type}_{max_z}.csv'
    redshift_filename = f'results/redshift_array_{type}_{max_z}.npy'

    # Check if redshift array file exists, else create it
    if os.path.exists(redshift_filename):
        redshift_array = np.load(redshift_filename)
        print(f"Loaded redshift array from {redshift_filename}\n")
        print('Simulating a total of {} {} transients between z = 0.001 and z = {} \n'.format(len(redshift_array), type, max_z))

    else:
        # Create a new redshift array and save it for consistency
        redshift_array = redshift_samples(type = type, z_max = max_z)
        np.save(redshift_filename, redshift_array)
        print(f"Generated and saved new redshift array to {redshift_filename}\n")
        print('Simulating a total of {} {} transients between z = 0.001 and z = {} \n'.format(len(redshift_array), type, max_z))

    processed_count = 0
    num_transients = len(redshift_array)
    ra_array, dec_array = ra_dec_samples(c_ra, c_dec, radius, num_transients)
    time_array = time_samples(365, num_transients)

    # Check if the results file already exists
    if os.path.exists(results_filename):
        # Load the existing results
        all_results = pd.read_csv(results_filename)
        print(f"Loaded existing results from {results_filename}\n")

        # Check which numbered transients are already processed
        processed_numbers = list(set(all_results['number']))
        remaining_numbers = list(set(np.arange(0, num_transients, 1)) - set(processed_numbers))
        print(f"{len(processed_numbers)} transients already processed, {len(remaining_numbers)} left to process.\n")

        # If all desired transients have been processed, exit
        if len(remaining_numbers) <= 0:

            print("All desired transients have already been processed.\n")
            return all_results  # Return the existing results
    else:
        # If the file doesn't exist, create an empty DataFrame
        all_results = pd.DataFrame(columns=['number', 'type', 'model', 'z', 'phase', 'filter', 'mag', 'mag_err', 'snr'])
        remaining_numbers = np.arange(0, num_transients, 1)
        print(f"No existing results found for {type} with maxz = {max_z}, starting fresh.\n")


    for number in remaining_numbers:

        try:
            result = process_single_redshift(number, type, models, redshift_array, ra_array, dec_array, time_array, cadence, MyTelescope, MyBackground)
            
            if len(all_results) != 0:
                all_results = pd.concat([all_results, result])
            else:
                all_results = result

            # Save the results incrementally to avoid losing progress
            all_results.to_csv(results_filename, index=False)
            print(f'Simulated transient {number+1}/{num_transients}\n')

        except Exception as e:
            print(f"Error processing transient {number+1}: {e}\n")
            print()
    # Reset index and save the final results to a CSV
    all_results = all_results.reset_index(drop=True)
    all_results.to_csv(results_filename, index=False)
    print(f"All results saved to {results_filename}\n")

    return all_results

