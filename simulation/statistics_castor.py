import numpy as np
import math
import pandas as pd
import glob
import os
from astropy.cosmology import FlatLambdaCDM 
from astropy import units as u
import extinction

# importing functions from other files
import utils_castor as utils
import rates_castor as rates
import simulation_functions_castor as simul


# Basic light curve detection condition, default conditions check whether there are a minimum of 2 detections with a SNR > 5
def lc_detected(lc, snr_lim = 5, n_det_above_snr = 2):
    if len(lc.loc[lc['snr'] >= snr_lim]) >= n_det_above_snr:
        return True
    else:
        return False
        
def lc_detected_plateau(lc, snr_lim = 2, n_det_above_snr = 1):
    if len(lc.loc[lc['snr_plt'] >= snr_lim]) >= n_det_above_snr:
        return True
    else:
        return False

def plateau_phase_detection(lc, snr_lim = 2):
    index = np.where(lc['phase']==lc['phase'][len(lc)-1])[0][0]
    if lc.loc[index, 'snr_plt']>= snr_lim:
       return lc.loc[index, 'phase'] 
        
def mag_plateau_detection(lc, snr_lim = 2):
    index = np.where(lc['phase']==lc['phase'][len(lc)-1])[0][0]
    if lc.loc[index, 'snr_plt']>= snr_lim:
       return lc.loc[index, 'mag_plt']   

# Extracts the phase of first detection assuming the SNR limit from phase 0 documents, although this can be changed
def first_detection(lc, snr_lim = 5):
    return np.nanmin(lc.loc[lc['snr'] >= snr_lim, 'phase'])

# Extract magnitude at first detection
def mag_first_detection(lc, snr_lim = 5):
    index = np.where(lc['phase']==np.nanmin(lc.loc[lc['snr'] >= snr_lim, 'phase']))[0][0]
    return lc.loc[index, 'mag']

# Function that fits polynomial around peak to extract the magnitude and time of peak
# Relies on initial estimate of peak from just finding the brighest point in the light curve and fits polynomial 7 days either side of this
def mag_at_peak(lc):
    est_t0 = lc['phase'][np.where(lc['mag']==np.nanmin(lc['mag']))[0][0]]
    lc_around_peak = lc.loc[(lc['phase'] < est_t0 + 7) & (lc['phase'] > est_t0 - 7)]

    lc_around_peak = lc_around_peak.sort_values(by = ['phase'])
    lc_around_peak = lc_around_peak.reset_index(drop=True)

    # Fitting polynomial around peak to extract peak mag and time
    z = np.polyfit(lc_around_peak['phase'].astype(float), lc_around_peak['mag'].astype(float), 2)
    p = np.poly1d(z)
    x = np.arange(est_t0 - 7, est_t0 + 7, 0.01)

    peak_mag = np.nanmin(p(x))
    t0 = x[np.where(p(x)==peak_mag)[0][0]]

    return peak_mag, t0


# More detailed light curve detection condition, adapt to individual needs
def lc_detected_useful(lc, snr_lim = 5, n_det_above_snr = 3):
    peak_time = mag_at_peak(lc)[1]
    pre_peak_lc = lc.loc[lc['phase'] < peak_time]

    # Checks how many detections pass SNR limit pre-peak, default requires 3 to pass 
    if len(pre_peak_lc.loc[pre_peak_lc['snr'] >= snr_lim]) >= n_det_above_snr:
        return True
    else:
        return False

# Function to calculate absolute peak magnitude (default in g band) by taking into account redshift and Milky Way extinction
def abs_peak_mag(mb, ebv, z, band = 'g', r_v = 3.1):

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dismod=cosmo.distmod(z).value

    a_v = r_v * ebv

    if band == 'g':
        wave = np.array([4708.0]) # Effective wavelengths in Angstrom taken from phase 0
    elif band == 'u':
        wave = np.array([3426.0])
    elif band == 'uv':
        wave = np.array([2209.0])

    filter_extinction = extinction.ccm89(wave, a_v, r_v)[0] # extinction in particular filter in units of magnitudes

    return mb - filter_extinction - dismod



# Function to process each light curve and pass it to all the above statistics functions
def process_light_curve(i, df, band, snr_lim=5, n_det_above_snr=2, plateau=False):
    if plateau==False:
        overview = pd.DataFrame(columns=['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'detected', 'detected_useful', 'phase_detected', 't0', 'mag_peak', 'abs_mag_peak', 'mag_detect'])
    elif plateau==True:
        overview = pd.DataFrame(columns=['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'detected', 'detected_useful', 'detected_plateau', 'phase_detected', 't0', 'mag_peak', 
                                         'abs_mag_peak', 'mag_detect', 'phase_plateau', 'mag_plateau'
                                         ])       
    
    # Locating a single light curve within the larger results file
    lc_all_filters = df.loc[(df['number'] == i)]
    lc_all_filters = lc_all_filters.reset_index(drop=True)

    lc = lc_all_filters.loc[(lc_all_filters['filter'] == band)]
    lc = lc.reset_index(drop=True)
    
    if len(lc) > 0:
        print(f'Processing statistics for light curve number {i}')
        try:
            # Extracting key parameters for transient
            ra, dec, ebv, z, t, m = lc['ra'][0], lc['dec'][0], lc['ebv'][0], lc['z'][0], lc['type'][0], lc['model'][0]
            
            # Passing light curve to basic and more complex detection checks (checks detections in any filter at this point)
            det = lc_detected(lc_all_filters)
            det_useful = lc_detected_useful(lc_all_filters)

            # If light curve has passed basic detection criteria, check all other statistics in specified filter
            if det:
                try:
                    first_det = first_detection(lc)
                    mag_first_det = mag_first_detection(lc)
                    mag_peak, t0 = mag_at_peak(lc)
                    abs_mag_peak = abs_peak_mag(mag_peak, ebv, z)
                    if plateau==True:
                        det_plateau = lc_detected_plateau(lc_all_filters)
                        if det_plateau:
                            mag_plateau_det = mag_plateau_detection(lc)
                            plt_phase_det = plateau_phase_detection(lc)
                        else:
                            mag_plateau_det = np.nan
                            plt_phase_det = np.nan
                    elif plateau==False: 
                        pass
                except ValueError:
                    first_det = np.nan
                    mag_first_det = np.nan
                    mag_peak, t0 = np.nan, np.nan
                    abs_mag_peak = np.nan
            else:
                first_det = np.nan
                mag_first_det = np.nan
                mag_peak, t0 = np.nan, np.nan
                abs_mag_peak = np.nan

            # Saving all detection statistics to overview file
            if plateau==False:
                overview.loc[0] = [i, t, m, z, ra, dec, ebv, det, det_useful, first_det, t0, mag_peak, abs_mag_peak, mag_first_det]
            elif plateau==True:
                overview.loc[0] = [i, t, m, z, ra, dec, ebv, det, det_useful, det_plateau, first_det, t0, mag_peak, abs_mag_peak, mag_first_det, plt_phase_det, mag_plateau_det]

        except Exception as e:
            print(f"Error processing light curve number {i}: {e}")
            return pd.DataFrame()

    return overview



# Function that iterates over the results file and passes each light curve to the process_light_curve function
# Checks if statistics have previously been run for any of the light curves and if so, picks up where it left off
def statistics(df, max_z, type, snr_lim=5, n_det_above_snr=2, checkpoint_interval=10, band='g', cadence = 1.0, exposure = 100, c_ra=9.45, c_dec=-44.0, test = False, number_redshifts = 10, plateau=False):

    # Check if redshift array file exists, else create it
    # redshift_filename = f'results/redshift_array_{type}_{max_z}_{c_ra}_{c_dec}.npy'
    
    # Load existing statistics if they exist
    overview_file = f'results/statistics_{type}_{max_z}_{band}_{cadence}d_{exposure}s_{c_ra}_{c_dec}.csv'

    if test == True:
        #redshift_filename = f'results/redshift_array_{type}_{max_z}_{c_ra}_{c_dec}_{number_redshifts}_test.npy'
        overview_file = f'results/statistics_{type}_{max_z}_{band}_{cadence}d_{exposure}s_{number_redshifts}_test.csv'

    #if os.path.exists(redshift_filename):
    #    redshift_array = np.load(redshift_filename)
    #    print(f"Loaded redshift array from {redshift_filename}\n")

    #else:
    #    print(f"Cannot find {redshift_filename}, likely light curves have not been simulated with these survey parameters. \n")

    if os.path.isfile(overview_file):
        overview = pd.DataFrame(pd.read_csv(overview_file))

        # Checks how many transients we need to run, and how many have already been processed
        num_transients = len(df)
        numbers_total = np.arange(0, num_transients, 1)
        numbers_completed = list(set(overview['number']))
        numbers = list(set(numbers_total) - set(numbers_completed))
        print(f'{len(numbers_completed)} already processed, processing remaining {len(numbers)} \n')

    else:
        overview = pd.DataFrame()
        numbers = list(set(df['number']))

    if test == False:


        # Iterates over remaining numbers to be processed and passes light curves to process_light_curve function
        for number in numbers:
            result = process_light_curve(number, df, band, snr_lim=snr_lim, n_det_above_snr=n_det_above_snr, plateau=plateau)
            
            # Appending results to overview file
            if len(overview) != 0:
                overview = pd.concat([overview, result], ignore_index = True)
            else:
                overview = result
    else:
         # Iterates over remaining numbers to be processed and passes light curves to process_light_curve function
        for number in list(set(df['number'])):
            result = process_light_curve(number, df, band, snr_lim=snr_lim, n_det_above_snr=n_det_above_snr)
            
            # Appending results to overview file
            if len(overview) != 0:
                overview = pd.concat([overview, result], ignore_index = True)
            else:
                overview = result

    # Save the final dataframe to CSV
    overview.to_csv(overview_file, index=False)
    return overview




