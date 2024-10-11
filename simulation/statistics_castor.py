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


# Basic light curve detection condition
def lc_detected_basic(lc, snr_lim = 5, n_det_above_snr = 2):
    if len(lc.loc[lc['snr'] >= snr_lim]) >= n_det_above_snr:
        return True
    else:
        return False

# More detailed light curve detection condition 
def lc_detected(lc, snr_lim = 5, n_det_above_snr = 2):
    if len(lc.loc[lc['snr'] >= snr_lim]) >= n_det_above_snr:
        return True
    else:
        return False

# Extracts the phase of first detection assuming the SNR limit from phase 0 documents, although this can be changed
def first_detection(lc, snr_lim = 5):
    return np.nanmin(lc.loc[lc['snr'] >= snr_lim, 'phase'])

# Extract magnitude at first detection
def mag_first_detection(lc, snr_lim = 5):
    index = np.where(lc['phase']==np.nanmin(lc.loc[lc['snr'] >= snr_lim, 'phase']))[0][0]
    return lc.loc[index, 'mag']

# Currently just provides an estimate of peak mag but not reliable, need to fit models before FORECASTOR and extract
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


# More detailed light curve detection condition 
def lc_detected_useful(lc, snr_lim = 5, n_det_above_snr = 3):
    peak_time = mag_at_peak(lc)[1]
    pre_peak_lc = lc.loc[lc['phase'] < peak_time]

    if len(pre_peak_lc.loc[pre_peak_lc['snr'] >= snr_lim]) >= n_det_above_snr:
        return True
    else:
        return False


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



# Modified statistics function for parallel processing
def process_light_curve(i, df, band, snr_lim=5, n_det_above_snr=2):
    overview = pd.DataFrame(columns=['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'detected', 'detected_useful', 'phase_detected', 't0', 'mag_peak', 'abs_mag_peak', 'mag_detect'])

    lc_all_filters = df.loc[(df['number'] == i)]
    lc_all_filters = lc_all_filters.reset_index(drop=True)

    lc = lc_all_filters.loc[(lc_all_filters['filter'] == band)]
    lc = lc.reset_index(drop=True)

    if len(lc) > 0:
        print(f'Processing statistics for light curve number {i}')

        try:
            ra, dec, ebv, z, t, m = lc['ra'][0], lc['dec'][0], lc['ebv'][0], lc['z'][0], lc['type'][0], lc['model'][0]
            det = lc_detected(lc_all_filters)
            det_useful = lc_detected_useful(lc_all_filters)

            if det:
                try:
                    first_det = first_detection(lc)
                    mag_first_det = mag_first_detection(lc)
                    mag_peak, t0 = mag_at_peak(lc)
                    abs_mag_peak = abs_peak_mag(mag_peak, ebv, z)
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

            overview.loc[0] = [i, t, m, z, ra, dec, ebv, det, det_useful, first_det, t0, mag_peak, abs_mag_peak, mag_first_det]

        except Exception as e:
            print(f"Error processing light curve number {i}: {e}")
            return pd.DataFrame()

    return overview




def statistics(df, max_z, type, snr_lim=5, n_det_above_snr=2, checkpoint_interval=10, band='g', n_cores=8):

    # Check if redshift array file exists, else create it
    redshift_filename = f'results/redshift_array_{type}_{max_z}.npy'

    if os.path.exists(redshift_filename):
        redshift_array = np.load(redshift_filename)
        print(f"Loaded redshift array from {redshift_filename}\n")

    else:
        # Create a new redshift array and save it for consistency
        redshift_array = simul.redshift_samples(type = type, z_max = max_z)
        np.save(redshift_filename, redshift_array)
        print(f"Generated and saved new redshift array to {redshift_filename}\n")

    # Load existing statistics if they exist
    overview_file = 'results/statistics_{}_{}_{}.csv'.format(type, max_z, band)
    if os.path.isfile(overview_file):
        overview = pd.read_csv(overview_file, names=['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'detected', 'detected_useful', 'phase_detected', 't0', 'mag_peak', 'abs_mag_peak', 'mag_detect'])
    else:
        overview = pd.DataFrame(columns=['number', 'type', 'model', 'z', 'ra', 'dec', 'ebv', 'detected', 'detected_useful', 'phase_detected', 't0', 'mag_peak', 'abs_mag_peak', 'mag_detect'])
    
    num_transients = len(redshift_array)
    numbers_total = np.arange(0, num_transients, 1)
    numbers_completed = list(set(overview['number']))

    numbers = list(set(numbers_total) - set(numbers_completed))

    print(f'{len(numbers_completed)} already processed, processing remaining {len(numbers)} \n')


    for number in numbers:
        result = process_light_curve(number, df, band, snr_lim=snr_lim, n_det_above_snr=n_det_above_snr)
        if len(overview) != 0:
            overview = pd.concat([overview, result])
        else:
            overview = result
    # Save the final dataframe to CSV
    overview.to_csv(overview_file, index=False)
    return overview




