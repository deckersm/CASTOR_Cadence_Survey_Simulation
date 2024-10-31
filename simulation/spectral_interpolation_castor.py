import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
from mpl_toolkits.mplot3d import Axes3D

filepath_to_castor_folder = '/users/deckersm/CASTOR/'
#filepath_to_castor_folder = '/Users/maximedeckers/Documents/RSE/CASTOR/CASTOR_Cadence_Survey_Simulation/'


# Function that produces a 2D surface in wavelength and phase space from existing spectral templates
def create_2d_surface(type, model, cadence = 0.1, plot = False):

    # Finding all the spectral templates for this particular transient type and model
    files = glob.glob(filepath_to_castor_folder + 'Templates/{}/SED_{}_{}_*d.dat'.format(type, type, model))

    # Saving all the fluxes and phases for this particular transient type and model
    phases, fluxes = [],[]
    for f in files:
        df = pd.DataFrame(pd.read_csv(f, names = ['lamb', 'flux'], comment = '#', sep = '\s+'))
        fluxes.append(df['flux'])
        phases.append(float(f.split('/')[-1].split('_')[3].replace('d.dat', '')))

    wavelengths = df['lamb'].values
    phases = np.asarray(phases)
    fluxes = np.asarray(fluxes)
    
    # Sort the phases and fluxes so that the phase array is always increasing (requirement for interpolation function)
    sorted_indices = np.argsort(phases)
    phases = phases[sorted_indices]
    fluxes = fluxes[sorted_indices, :]

    # Interpolate across both phases and wavelengths    
    interp_func = RectBivariateSpline(phases, wavelengths, fluxes)

    # Interpolating to new phases and the original wavelength grid (default assumes new cadence = 1.0 d)
    fine_phases = np.arange(np.round(np.nanmin(phases), 0), np.round(np.nanmax(phases), 0), cadence)
    fine_phases = np.round(fine_phases, 1)
    fine_wavelengths = wavelengths  # Can keep the wavelength grid or interpolate here too e.g. fine_wavelengths = np.arange(0, 10000, 100)

    # Generate interpolated fluxes
    interpolated_fluxes = interp_func(fine_phases, fine_wavelengths)

    # Plotting the 2D surface if plot = True
    if plot == True:
        
        # Create a meshgrid for plotting (phases and wavelengths)
        P, W = np.meshgrid(fine_phases, fine_wavelengths)
        
        # Transpose the flux data to match the meshgrid shape
        Z = interpolated_fluxes.T
        
        # Plotting the 2D surface using 3D plot
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot the surface
        surf = ax.plot_surface(W, P, Z, cmap='viridis')
        
        # Adding labels and a color bar
        ax.set_ylabel('Phase (Days)')
        ax.set_xlabel('Wavelength (Angstrom)')
        ax.set_zlabel('Flux')
        ax.set_title('Interpolated Spectral Surface')
        ax.set_xlim(0, 10000)
        # Add a color bar for flux values
        fig.colorbar(surf)
        
        plt.show()
        

    return fine_phases, fine_wavelengths, interpolated_fluxes


# Function that will output a spectrum at any given phase (integer) by pulling from the 2D interpolated spectral surface
def create_spec_at_phase(type, model, phase):
    
    # Creating the 2D surface based on existing templates
    phases, wavelengths, fluxes = create_2d_surface(type, model)

    # Finding the index of the required phase and extracting the spectrum at this phase
    index = np.where(list(phases) == np.round(phase, 1))[0][0]
    spectrum = pd.DataFrame(data = {'#lamb':wavelengths, 'flux':fluxes[index]})

    return spectrum

