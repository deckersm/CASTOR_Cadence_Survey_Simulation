# CASTOR_Cadence_Survey_Simulation

## Description
This software simulates transients for the CASTOR Cadence Survey: A Wide-Field UV Time Domain Survey (https://www.castormission.org/). This software relies on the usage of FORECASTOR ETC (https://github.com/CASTOR-telescope/ETC) to determine the SNR of each spectrum at a given redshift. 



## Features
- Simulates light curves for transient events at different redshifts by feeding spectral time series into FORECASTOR ETC.
- Incorporates volumetric rate estimates for all the different transient types to inject a realistic number into the simulation at realistic redshift distributions.
- Outputs detection statistics based on signal-to-noise ratios of detections.
- Includes the option to run a simplified version of the simulation to obtain redshift-dependent detection efficiencies. 

## Installation
Clone this repository:
git clone https://github.com/yourusername/yourproject.git

pip install -r requirements.txt

## Usage
To run a simulation:

python castor_simulation.py <transient_type> <max_redshift> <filter_for_statistics> <survey_cadence>

## 6. Configuration
- The software is currently set up to run a simulation based on the CASTOR survey parameters, these can be changed to simulate other surveys. To change the telescope parameters, we point the user to the FORECASTOR ETC (https://github.com/CASTOR-telescope/ETC) for further details. These can then be directly implemented in the `config_telescope()` function in `utils.py` file where all the FORECASTOR ETC functions are called.
- The current version includes templates for: SNe Ia, II, Ibc, Iax, 91bg, SLSNe, and TDEs, but other classes of transients can be added by the user given that the spectral templates have sufficient wavelength coverage (1000 \AA -- 10 000 \AA). Transient templates should be stored in the `/Templates/` folder and follow the naming convention `SED_{type}_{model}_{phase}.dat`. The template library should be representative of the intrinsic transient population (e.g. in terms of the luminosity distribution for most transient classes). 
- The code runs the simulation directly from the spectral templates, rather than first producing models based on these templates. As a result, it is limited to running the cadence of the spectral templates (or integer multiples of this cadence). In future version we will include an interpolation function that enables the used to interpolate the spectral time series so that the survey simulations can be run at any desired cadence.
  



## Contact
If you have any questions, feel free to open an issue or contact me at [maxime.deckers@port.ac.uk](mailto:maxime.deckers@port.ac.uk).
