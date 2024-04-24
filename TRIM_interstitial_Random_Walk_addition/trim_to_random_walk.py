import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
from scipy.interpolate import interp1d
import load_trim_insert_interstitials as lti
from copy import deepcopy
####TODO: the first interstitial is at the same place where there was a vacancy
def calculate_interstitials(track_data, energy_data, energy_distance_file, initial_position):
    # Load the energy to distance mapping
    energy_distance_df = pd.read_csv(energy_distance_file)
    energy_distance_df.rename(columns={'Column_A': 'Distance (A)'}, inplace=True)

    # Setup interpolation function for energy to distance
    energy_interpolator = interp1d(energy_distance_df['Energy (eV)'], energy_distance_df['Distance (A)'],
                                   kind='linear', fill_value='extrapolate')

    # Calculate displacement vectors from the initial position
    initial_position = initial_position
    vectors = track_data - initial_position

    # Convert from micrometers to angstroms (1 um = 10000 A)
    vectors *= 10000  # for energy-to-distance scaling in angstroms

    # Calculate distances for each energy value
    distances = energy_interpolator(energy_data)

    # Normalize vectors to unit vectors
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    unit_vectors = np.where(norms != 0, vectors / norms, 0)

    # Calculate interstitial positions (in angstroms)
    interstitial_positions_A = track_data * 10000 + unit_vectors * distances[:, None]

    # Convert positions from angstroms back to micrometers
    interstitial_positions_um = interstitial_positions_A / 10000

    return interstitial_positions_um

particle = "Si"
# expects something like python trim_to_random_walk.py ../Silicon/0.001
parser = argparse.ArgumentParser(
    description="Analyze TRIM output. \
    Requires 1 argument: folder containing the TRIM output files."
)
folder = parser.add_argument(
    "folder", metavar="str", type=str, nargs="+", help="folder to load."
)
args = parser.parse_args()
folder = args.folder[0]
PKA, no_inc = lti.load_data(folder, particle)

for pka in PKA:
    track = pka.recoil_positions
    track.append(pka.ion_position)
    track = np.asarray(track)  
    energy_i = np.append(pka.recoil_energies, pka.last_ion_energy)
    interstitial_track = np.copy(track)
    interstitial_track[-1] = pka.last_ion_position
    interstitial_positions = calculate_interstitials(interstitial_track, energy_i, "Energy_to_distance.csv", initial_position=pka.ion_position)
    x, y, z = zip(*track)
    x_i, y_i, z_i = zip(*interstitial_positions)
    for array in [x_i, y_i, z_i, x, y, z]:
        array = np.asarray(array)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z,"o", label='Vacancies', color='blue')
    ax.plot(x_i, y_i, z_i,"x", label='Interstitials', color='red')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.show()
