# CASTOR_Cadence_Survey_Simulation

## Description
This package simulates transients for the CASTOR Cadence Survey: A Wide-Field UV Time Domain Survey (https://www.castormission.org/) and relies on the usage of FORECASTOR ETC (https://github.com/CASTOR-telescope/ETC) to determine the SNR of each spectrum at a given redshift. 



## Features
- Simulates light curves for transient events at different redshifts by feeding spectral time series into FORECASTOR ETC.
- Incorporates volumetric rate estimates for all the different transient types to inject a realistic number into the simulation at realistic redshift distributions.
- Outputs detection statistics based on signal-to-noise ratios of detections.
- Includes the option to run a simplified version of the simulation to obtain redshift-dependent detection efficiencies. 

## Installation
Clone this repository:
git clone [https://github.com/deckersm/CASTOR_Cadence_Survey_Simulation.git](https://github.com/deckersm/CASTOR_Cadence_Survey_Simulation.git)

pip install -r requirements.txt

## Usage
To run a simulation, make sure a `results` folder is set up inside the `simulation` folder and run the following command:

python main.py <<arguments>>

Required arguments:

`--type`, `-t` = transient type

`--max_redshift`, `-z` = maximum redshift out to which transients are simulated


Optional arguments:

`--cadence`, `-c` = cadence of transients survey (default = 1.0 d).

`--exposure`, `-e` = exposure time of each visit (default = 100 s).

`--length_survey`, `-l` = duration of the survey in days (default = 182.625 d).

`--simul_type`, `-s` = 'full' runs a complete simulation with transients injected according to their volumetric rate, 'test' runs a small test survey just to extract detection efficiencies as a function of redshift (default = 'full').

`--number_redshifts` = the number of redshifts samples to populate between min_z and max_z if a 'test' simulation is ran (default = 10).

`--min_redshift`, `--min_z` = minimum redshift to injecting transients at - only used for 'test' simulations (default = 0.01).

`--field_ra`, `-r` (`--field_dec`, `-d`) = ra and dec (in deg) of the center of the fields covered by the survey, if multiple, separate with spaces.

'--plateau' = Caculate detection for TDEs at the plateau. It only works when type = TDE (default = False).


This will produce several output files in the `results` folder:

- A `redshift_array_<transient_type>__<max_redshift>.npy` file which contains the redshifts of the transients simulated.
- A `results_<transient_type>_<max_redshift>_<cadence>d_<exposure_time>s_<ra>_<dec>.csv` file which contains the light curves of all the transients simulated.
- a `statistics_<transient_type>_<max_redshift>_<filter>_<cadence>d_<exposure_time>s_<ra>_<dec>.csv` file which contains the detection statistics for all the generated transients in a particular filter.

- if a 'test' version of the simulation is ran, these filenames will have an additional 'test' suffix

Detailed descriptions of the contents of the latter two files are available in `docs/file_descriptors.txt`. 


## Configuration
- The software is currently set up to run a simulation based on the CASTOR survey parameters, these can be changed to simulate other surveys. To change the telescope parameters, we point the user to the FORECASTOR ETC (https://github.com/CASTOR-telescope/ETC) for further details. These can then be directly implemented in the `config_telescope()` function in `utils.py` file where all the FORECASTOR ETC functions are called.
- The current version includes templates for: SNe Ia, II, Ibc, Iax, 91bg, SLSNe, and TDEs, but other classes of transients can be added by the user given that the spectral templates have sufficient wavelength coverage (1000 \AA -- 10 000 \AA). Transient templates should be stored in the `/Templates/` folder and follow the naming convention `SED_{type}_{model}_{phase}.dat`. The template library should be representative of the intrinsic transient population (e.g. in terms of the luminosity distribution for most transient classes). 
- The code runs the simulation directly from the spectral templates, rather than first producing models based on these templates. As a result, it is limited to running the cadence of the spectral templates (or integer multiples of this cadence). In future version we will include an interpolation function that enables the used to interpolate the spectral time series so that the survey simulations can be run at any desired cadence.
  



## Contact
If you have any questions, feel free to open an issue or contact me at [maxime.deckers@port.ac.uk](mailto:maxime.deckers@port.ac.uk).
