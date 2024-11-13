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


# Function which extracts redshifts based on the volumetric rate of a particular transient type
# Defaults assume the survey runs for 1 year and covers 20 deg^2 of the sky
def redshift_samples(type = 'snia', z_min = 0.00, z_max = 0.5, survey_time = 365.25, survey_radius = 1.75):
    
    # Adding units to survey time
    survey_time = survey_time * u.yr

    # Convert sky area to a fraction of the total sky (41253 deg² in total)
    survey_area = math.pi * survey_radius**2
    fraction_sky = survey_area / 41253

    # Number of redshift bins
    z_bins = np.linspace(z_min, z_max, 1000)

    # Extract the transient rate for the type, converting units to Mpc^-3 yr^-1
    rate = rates.transient_rate(type, z_bins)
    rate = rate * u.Mpc**(-3) * u.yr**(-1)

    # Calculate the comoving volume element dV/dz (Mpc³/sr) for each redshift bin
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dVdz = cosmo.differential_comoving_volume(z_bins)  # Mpc³/sr

    # Multiply by 4π steradians and the sky fraction to get comoving volume over the survey area
    dVdz_survey = dVdz * 4 * np.pi * fraction_sky  # Mpc³

    # Apply time dilation correction to the survey time in each redshift bin
    # Effective survey time in each bin is scaled by (1 + z)
    survey_time_years = survey_time / 365.25
    time_dilation = survey_time_years / (1 + z_bins)

    # Integrate to find the total comoving volume surveyed up to z_max, weighted by the rate
    total_volume = np.trapz(dVdz_survey * rate * time_dilation, z_bins)  # Mpc³ * yr^-1

    # Calculate the total number of transients over the survey period
    N_total = total_volume.decompose().value  # Dimensionless (total number of transients)

    # Probability distribution proportional to dV/dz * rate, normalized to 1
    pdf = (dVdz * rate * time_dilation).value
    pdf /= np.sum(pdf)  # Normalize to get a proper probability distribution

    # Sample redshifts according to the probability distribution
    N_samples = int(N_total)  # Total number of transients to sample
    sampled_redshifts = np.random.choice(z_bins, size=N_samples, p=pdf)

    return sampled_redshifts


# Function to pull an array of random ra and dec combination within a specified radius of the center of a field to match the array of redshifts
# Currently assumes field to be square, not circular but can easily be adapted 
def ra_dec_samples(c_ra, c_dec, radius, number):

    ra = np.random.uniform(c_ra - radius, c_ra + radius, number)
    dec = np.random.uniform(c_dec - radius, c_dec + radius, number)

    return ra, dec

# Function to pull random start times for each transients from before the start of the survey (default = 30 d) until the end of the survey
def time_samples(survey_time, number, time_before_survey = 30, start_time = 0):
    time_array = np.random.uniform(-time_before_survey + start_time, survey_time + start_time, number)
    return time_array



# Preload the SFD map at the start and reuse it for all calculations
#dustmaps.sfd.fetch()  # Ensure dust maps are fetched locally
sfd = SFDQuery()  # This loads the SFD map

# Function to extract Milky Way extinction along the line of sight of the transient 
def mw_ebv(ra, dec):
    coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
    ebv = sfd(coords)  # Reuse the preloaded sfd object
    return ebv
    


# Function the absorbs redshift, ra, dec, time, ebv, and outputs the model light curve
def process_single_redshift(i, type, models, redshift_array, ra_array, dec_array, time_array, cadence, exposure, MyTelescope, MyBackground):
    z = redshift_array[i]
    ra, dec = ra_array[i], dec_array[i]
    time = time_array[i]

    ebv = mw_ebv(ra, dec)  # Calculate extinction at this RA and Dec
    model = models[int(np.random.uniform(0, len(models), 1))]
    results = utils.create_lc(type, model, z, ra, dec, ebv, i, MyTelescope, MyBackground, start_time = time, cadence = cadence, exposure = exposure)
    return results


# Function that iterates of the redshift array created and appends a light curve for each transient to the results file
# This will check whether a results file already exists and continue where it last left off if some results were already produced
# Also checks for existing redshift array file to ensure we don't start from scratch every time the function is ran
def populate_redshift_range(type, models, max_z, MyTelescope, MyBackground, min_z = 0, c_ra=9.45, c_dec=-44.0, radius=1.75, cadence = 1.0, exposure = 100, survey_time = 365.25, start_time = 0, test = False):
    results_filename = f'results/results_{type}_{max_z}_{cadence}d_{exposure}s_{c_ra}_{c_dec}.csv'
    redshift_filename = f'results/redshift_array_{type}_{max_z}_{c_ra}_{c_dec}.npy'

    # Check if redshift array file exists, else create it
    if os.path.exists(redshift_filename):
        redshift_array = np.load(redshift_filename)
        print(f'Loaded redshift array from {redshift_filename}\n')
        print(f'Simulating a total of {len(redshift_array)} {type} transients between z = 0.001 and z = {max_z} \n')

    else:
        # Create a new redshift array and save it for consistency
        redshift_array = redshift_samples(type = type, z_min = min_z, z_max = max_z, survey_time = survey_time, survey_radius = radius)
        np.save(redshift_filename, redshift_array)
        print(f'Generated and saved new redshift array to {redshift_filename}\n')
        print(f'Simulating a total of {len(redshift_array)} {type} transients between z = 0.001 and z = {max_z} \n')

    # Initialise processed_count and checks how many transients need to be generated
    processed_count = 0
    num_transients = len(redshift_array)
    
    # Producing the ra, dec, and time arrays
    ra_array, dec_array = ra_dec_samples(c_ra, c_dec, radius, num_transients)
    time_array = time_samples(survey_time, num_transients, start_time = start_time)

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
        #numbers = np.arange(0, num_transients, 1)

        print(f"No existing results found for {type} with maxz = {max_z}, starting fresh.\n")

    # Iterating over the remaining_numbers (i.e. those that have not been previously processed)
    for number in remaining_numbers:
        try:
            # Produce a light curve based on the input
            result = process_single_redshift(number, type, models, redshift_array, ra_array, dec_array, time_array, cadence, exposure, MyTelescope, MyBackground)
            
            # Append the light curve to the results file
            if len(all_results) != 0:
                all_results = pd.concat([all_results, result])
            else:
                all_results = result

            # Save the results incrementally to avoid losing progress
            all_results.to_csv(results_filename, index=False)
            print(f'Simulated transient {number+1}/{num_transients}\n')

        except Exception as e:
            print(f"Error processing transient {number+1}: {e}\n")

    # Reset index and save the final results to a CSV
    all_results = all_results.reset_index(drop=True)
    all_results.to_csv(results_filename, index=False)
    print(f"All results saved to {results_filename}\n")

    return all_results

# Function the absorbs redshift, ra, dec, time, ebv, and outputs the model light curve
def process_single_redshift_test(i, type, model, redshift, ra_array, dec_array, time_array, MyTelescope, MyBackground):
    ra, dec = ra_array[i], dec_array[i]
    time = time_array[i]
    ebv = mw_ebv(ra, dec)  # Calculate extinction at this RA and Dec
    results = utils.create_lc(type, model, redshift, ra, dec, ebv, i, MyTelescope, MyBackground, start_time = time)
    return results

# Function that iterates of the redshift array created and appends a light curve for each transient to the results file
# This will check whether a results file already exists and continue where it last left off if some results were already produced
# Also checks for existing redshift array file to ensure we don't start from scratch every time the function is ran
def populate_redshift_range_test(type, models, max_z, MyTelescope, MyBackground, min_z = 0.01, c_ra=9.45, c_dec=-44.0, radius=1.75, cadence = 1.0, exposure = 100, survey_time = 365.25, number_redshifts = 10):

    results_filename = f'results/results_{type}_{max_z}_{cadence}d_{exposure}s_{number_redshifts}_test.csv'
    redshift_filename = f'results/redshift_array_{type}_{max_z}_{c_ra}_{c_dec}_{number_redshifts}_test.npy'

    # Check if redshift array file exists, else create it
    if os.path.exists(redshift_filename):
        redshift_array = np.load(redshift_filename)
        print(f'Loaded redshift array from {redshift_filename}\n')
        print(f'Simulating a total of {len(redshift_array)} {type} transients between z = 0.001 and z = {max_z} \n')

    else:
        # Create a new redshift array and save it for consistency
        redshift_array = np.linspace(min_z, max_z, number_redshifts)
        np.save(redshift_filename, redshift_array)
        print(f'Generated and saved new redshift array to {redshift_filename}\n')
        print(f'Simulating a total of {len(redshift_array)} {type} transients between z = 0.001 and z = {max_z} \n')

    # Initialise processed_count and checks how many transients need to be generated
    processed_count = 0
    num_transients = len(redshift_array) * len(models)
    # Producing the ra, dec, and time arrays
    ra_array, dec_array = ra_dec_samples(c_ra, c_dec, radius, num_transients)
    time_array = time_samples(survey_time, num_transients)

    # Check if the results file already exists
    if os.path.exists(results_filename):
        # Load the existing results
        all_results = pd.read_csv(results_filename)
        print(f"Loaded existing results from {results_filename}\n")

        # Check which numbered transients are already processed
        processed_numbers = list(set(all_results['number']))
        remaining_numbers = (len(models) * number_redshifts) - len(processed_numbers)
        print(f"{len(processed_numbers)} transients already processed, {remaining_numbers} left to process.\n")

        # If all desired transients have been processed, exit
        if remaining_numbers <= 0:

            print("All desired transients have already been processed.\n")
            return all_results  # Return the existing results
            
    else:
        # If the file doesn't exist, create an empty DataFrame
        all_results = pd.DataFrame(columns=['number', 'type', 'model', 'z', 'phase', 'filter', 'mag', 'mag_err', 'snr'])
        print(f"No existing results found for {type} with maxz = {max_z}, starting fresh.\n")

    count = 0
    try:
        for z in redshift_array:
            for m in models:

                if len(all_results.loc[(all_results['z']==z) & (all_results['model']==m)]) == 0:
                    result = process_single_redshift_test(count, type, m, z, ra_array, dec_array, time_array, MyTelescope, MyBackground)
                    # Append the light curve to the results file
                    if len(all_results) != 0:
                        all_results = pd.concat([all_results, result])
                    else:
                        all_results = result
                    all_results.to_csv(results_filename, index=False)

                print(f'Simulated transient {count+1}/{num_transients}\n')
                count += 1

    except Exception as e:
        print(f"Error processing transient {count+1}: {e}\n")

    # Reset index and save the final results to a CSV
    all_results = all_results.reset_index(drop=True)
    all_results.to_csv(results_filename, index=False)
    print(f"All results saved to {results_filename}\n")

    return all_results

