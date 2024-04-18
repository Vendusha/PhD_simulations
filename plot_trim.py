import matplotlib.pyplot as plt
import numpy as np
import load_trim
import seaborn as sns
import csv
import argparse
import _pickle as cPickle
import optics
import os
import math
from collections import OrderedDict
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
import json
import trim_helper as th
import trim_plotting_helper as tph
particle = "Si"
#particle = "Ian"
parser = argparse.ArgumentParser(
    description="Analyze TRIM output. \
    Requires 4 arguments: folder, n, eps/xi, method=dbscan/xi  \
    Folder must containing exactly 1.txt and 1.dat file"
)
folder = parser.add_argument(
    "folder", metavar="str", type=str, nargs="+", help="folder to load."
)
n = parser.add_argument("n", metavar="int", type=int, nargs="+", help="min_samples")
xi = parser.add_argument("xi", metavar="float", type=float, nargs="+", help="eps/xi")
method = parser.add_argument(
    "method", metavar="str", type=str, nargs="+", help="method=dbscan/xi"
)
args = parser.parse_args()
folder = args.folder[0]
n = args.n[0]
xi = args.xi[0]
method = args.method[0]
PKA, no_inc = load_trim.load_data(folder, particle)
# generate 15 differen colors
colors_proton_number = {
    "He": "blue",
    "Li": "green",
    "Be": "green",
    "B": "green",
    "C": "green",
    "Ni": "green",
    "O": "green",
    "F": "green",
    "Ne": "green",
    "Na": "green",
    "Mg": "green",
    "Mg": "yellow",
    "Al": "orange",
    "Si": "red",
    "P": "pink",
}
colors_proton_number = dict(colors_proton_number)
colors = list(colors_proton_number.values())
x_3D, y_3D, z_3D = [], [], []
color_array = []
color_array_cluster = []
total_count = 0
population_counts = [
    0,
    0,
    0,
    0,
    0,
]  # vacancy, divacancy, trivacancy, tetravacancy, >4vacancy
# setting the depth projection graphs
#depth_samples =[0,9,49,80]
depth_samples = []
#depth_samples = [9]
depth_samples_len = len(depth_samples)
x_section = [[] for i in range(depth_samples_len)]
y_section = [[] for i in range(depth_samples_len)]
z_section = [[] for i in range(depth_samples_len)]
x_centers_section = [[] for i in range(depth_samples_len)]
y_centers_section = [[] for i in range(depth_samples_len)]
z_centers_section = [[] for i in range(depth_samples_len)]
population_section = [[] for i in range(depth_samples_len)]
energy_section = [[] for i in range(depth_samples_len)]
density_ball_section = [[] for i in range(depth_samples_len)]
density_cuboid_section = [[] for i in range(depth_samples_len)]
color_density_section = [[] for i in range(depth_samples_len)]
divacancy_section = [[] for i in range(depth_samples_len)]
trivacancy_section = [[] for i in range(depth_samples_len)]
tetravacancy_section = [[] for i in range(depth_samples_len)]
pentavacancy_section = [[] for i in range(depth_samples_len)]
counts_section = [Counter() for i in range(depth_samples_len)]
isolated_count_section = [0 for i in range(depth_samples_len)]
section_count = [0 for i in range(depth_samples_len)]
color_origins_section = [[] for i in range(depth_samples_len)]
color_cluster_section = [[] for i in range(depth_samples_len)]
cluster_center_section = [[] for i in range(depth_samples_len)]
parent_id = 0
process_count = 0
count_inside_beam = 0
count_outside_beam = 0
count_outside_beam_upper = 0
count_outside_beam_lower = 0
########
Cluster_populations_counter = Counter()
Defects_isolated = 0
Defets_clustered = 0
No_clusters = 0
Defects_total = 0
Cluster_populations = []
Cluster_volumes_cuboid = []
Cluster_lengths_cuboid = []
Cluster_radius = []
Clusters_distance = []
Clusters_distance_x = []
Isolated_defects_distance_x = []
Isolated_defects_distance = []
Energy_isolated = []
Energy_clustered = []
Cluster_energies = []
Cluster_densities_cuboid = []
Cluster_densities_sphere = []
Energy_primary_recoil = []
Phonon_energy_ion = []
Phonon_energy_recoil = []
#Energy_all = []
for pka in PKA:
    # pka.flatten()
    track = pka.recoil_positions
    # if np.size(track) == 0:
    #     print("trouble")
    #     continue
    track.append(pka.ion_position)
    track = np.asarray(track)
    #Energy_all.extend(pka.recoil_energies)
    # if some item from pka.recoil_energies is bigger then 20 000, print the energy out
    pka.recoil_energies = np.append(pka.recoil_energies, pka.ion_energy)
    pka.add_max_min() 
    # uncomment after debugging
    pka.cluster_analysis(track, n, xi, method, showPlot=False, cutoff_1=0, cutoff_2=1)
    if np.size(pka.cluster_populations) > 0:
        Cluster_populations_counter += Counter(pka.cluster_populations)
    Defects_isolated += pka.defects_isolated
    Defets_clustered += pka.defects_clustered
    No_clusters += pka.no_clusters
    Defects_total += pka.defects_isolated + pka.defects_clustered
    Cluster_populations.extend(pka.cluster_populations)
    Cluster_volumes_cuboid.extend(pka.cluster_volumes_cuboid)
    Cluster_lengths_cuboid.extend(pka.cluster_lengths_cuboid)
    Cluster_radius.extend(pka.cluster_radius)
    #print("Printing 1 event")
    #print(pka.isolated_defect_distance)
    Isolated_defects_distance.extend(pka.isolated_defect_distance)
    Isolated_defects_distance_x.extend(pka.isolated_defect_distance_x)
    Clusters_distance.extend(pka.clustered_defect_distance)
    Clusters_distance_x.extend(pka.clustered_defect_distance_x)
    Energy_isolated.append(pka.energy_isolated)
    Energy_clustered.append(pka.energy_clustered)
    Cluster_energies.extend(pka.cluster_energies)
    Cluster_densities_cuboid.extend(pka.cluster_densities_cuboid)
    Cluster_densities_sphere.extend(pka.cluster_densities_sphere)
    Energy_primary_recoil.append(pka.energy_primary_recoils)
    x, y, z = zip(*track)
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    x_3D.append(x)
    y_3D.append(y)
    z_3D.append(z)
    color_array.append([colors[pka.ion_origin - 2]] * np.size(x))
    color_array_cluster.append(pka.cluster_colors)
    total_count += np.size(x)
    # plot 3D track, title is the total count, colors are color_array_cluster colors
    #tph.plot_single_event(x, y, z, pka.cluster_colors, total_count)
    indeces_outside_beam_upper = np.argwhere((0 < y) & (y < 1) & (1 < z) & (z < 2))
    indeces_outside_beam_lower = np.argwhere((0 < y) & (y < 1) & (-1 < z) & (z < 0))
    if indeces_outside_beam_upper.size != 0:
        count_outside_beam_upper += indeces_outside_beam_upper.size
    if indeces_outside_beam_lower.size != 0:
        count_outside_beam_lower += indeces_outside_beam_lower.size
    indeces_inside_beam = np.argwhere((0 < y) & (y < 1) & (0 < z) & (z < 1))
    if indeces_inside_beam.size != 0:
        count_inside_beam += indeces_inside_beam.size
        count_outside_beam += np.size(x) - indeces_inside_beam.size
    else:
        count_outside_beam += np.size(x)
    for i in range(0, depth_samples_len):
        indeces = np.argwhere(
            (depth_samples[i] < x)
            & (x < depth_samples[i] + 1)
            & (0 < y)
            & (y < 1)
            & (0 < z)
            & (z < 1)
        )
        if indeces.size != 0:
            x_section[i].extend(x[indeces].ravel())
            y_section[i].extend(y[indeces].ravel())
            z_section[i].extend(z[indeces].ravel())
            color_origins_section[i].extend([colors[pka.ion_origin - 2]] * indeces.size)
            track_section = np.asarray([x[indeces], y[indeces], z[indeces]]).T
            if len(pka.recoil_energies) != len(pka.recoil_energies):
                print(
                    "Error: Number of recoil energies does not match number of recoil positions"
                )
            pka_section = load_trim.RecoilsFrom1PKA()
            pka_section.add_origins(
                pka.primary_event, pka.ion_origin, pka.ion_energy, pka.ion_position
            )
            indeces = np.asarray(indeces).ravel()
            pka_section.recoil_energies = np.take(pka.recoil_energies, indeces)
            pka_section.recoil_positions = track_section
            pka_section.cluster_analysis(
                track_section[0], n, xi, method, showPlot=False, cutoff_1=0, cutoff_2=1
            )
            if len(pka_section.divacancy) > 0:
                divacancy_section[i].extend(pka_section.divacancy)
            if len(pka_section.trivacancy) > 0:
                trivacancy_section[i].extend(pka_section.trivacancy)
            if len(pka_section.tetravacancy) > 0:
                tetravacancy_section[i].extend(pka_section.tetravacancy)
            if len(pka_section.pentavacancy) > 0:
                pentavacancy_section[i].extend(pka_section.pentavacancy)
            isolated_count_section[i] += pka_section.defects_isolated
            # color_cluster_section[i].append(pka.cluster_colors[i] for i in indeces) good for later when doing the whole analysis
            color_cluster_section[i].extend(pka_section.cluster_colors)
            section_count[i] += indeces.size
            counts_section[i] += Counter(pka_section.cluster_populations)
            if len(pka_section.cluster_populations) > 0:
                x_centers, y_centers, z_centers = zip(*(pka_section.cluster_centers))               
                x_centers = np.asarray(x_centers)
                y_centers = np.asarray(y_centers)
                z_centers = np.asarray(z_centers)
                x_centers_section[i].extend(x_centers.ravel())
                y_centers_section[i].extend(y_centers.ravel())
                z_centers_section[i].extend(z_centers.ravel())
                population_section[i].extend(pka_section.cluster_populations)
                energy_section[i].extend(np.asarray(pka_section.cluster_energies)/1000)
                density_ball_section[i].extend(pka_section.cluster_densities_sphere)
                density_cuboid_section[i].extend(pka_section.cluster_densities_cuboid)
                #marker_population_section[i].extend(pka.cluster_populations[indeces_cluster_centers].ravel())
                #color_energy_section[i].extend(pka.cluster_energies[indeces_cluster_centers].ravel())
                #marker_energy_section[i].extend(pka.cluster_energies[indeces_cluster_centers].ravel())
                #energy_section[i].extend(pka.cluster_energies[indeces_cluster_centers].ravel())
                #density_section[i]
    process_count += 1
    if folder.split("_")[-1] == "electrons" or folder.split("_")[-1] == "electrons/":
        if process_count % 2000 == 0:
            print(str(process_count) + "/" + str(len(PKA)) + " events processed.")
    else:
        if process_count % 20 == 0:
            print(str(process_count) + "/" + str(len(PKA)) + " events processed.")
x_3D = np.concatenate(x_3D).ravel()
y_3D = np.concatenate(y_3D).ravel()
z_3D = np.concatenate(z_3D).ravel()
color_array = np.concatenate(color_array).ravel()
color_array_cluster = np.concatenate(color_array_cluster).ravel()
print(Defects_total/1000)
print(Defects_isolated/Defects_total)
# print(x_section[i])
####################starting to add additional tracks####################
#tph.plot(
#    depth_samples, depth_samples_len, x_section, y_section, z_section, 
#    color_origins_section, counts_section, isolated_count_section, section_count, 
#    x_3D, y_3D, z_3D, color_array, total_count, count_inside_beam, count_outside_beam, 
#    count_outside_beam_upper, count_outside_beam_lower, divacancy_section, trivacancy_section, folder, x_centers_section,
#    y_centers_section, z_centers_section, population_section, energy_section, density_ball_section, density_cuboid_section, "_before_correction")
# plt.figure()
# #plot histogram of the energies
# plt.hist(Energy_all, bins=10000, alpha=0.5, color='b')
# np.savetxt("Energy_recoil_distribution.txt", Energy_all)
# plt.xlabel("Energy (eV)")
# plt.ylabel("Counts")
# plt.yscale('log')
# plt.show()
dict = {
    "Cluster_populations_counter": Cluster_populations_counter,
    "Defects_isolated": Defects_isolated,
    "Defects_clustered": Defets_clustered,
    "No_clusters": No_clusters,
    "Defects_total": Defects_total,
    "Cluster_populations": Cluster_populations,
    "Cluster_volumes_cuboid": Cluster_volumes_cuboid,
    "Cluster_lengths_cuboid": Cluster_lengths_cuboid,
    "Cluster_radius": Cluster_radius,
    "Clusters_distance": Clusters_distance,
    "Clusters_distance_x": Clusters_distance_x,
    "Isolated_defects_distance": Isolated_defects_distance,
    "Isolated_defects_distance_x": Isolated_defects_distance_x,
    "Energy_isolated": Energy_isolated,
    "Energy_clustered": Energy_clustered,
    "no_inc": no_inc,
    "Cluster_energies": Cluster_energies,
    "Cluster_densities_sphere": Cluster_densities_sphere,
    "Cluster_densities_cuboid": Cluster_densities_cuboid,
    "Energy_primary_recoil": Energy_primary_recoil
}
json_data = json.dumps(dict)
with open(str(folder) + "/total_count_ion_included"+str(xi)+".json", "w") as f:
    f.write(json_data)
    
# print("Starting to add additional tracks")
# total_norm_no_tracks = 100000000
# #total_norm_no_tracks = 100000
# if total_norm_no_tracks < no_inc:
#     print("Error: total_norm_no_tracks<no_inc")
#     exit()
# generate_additional_tracks = False
# if generate_additional_tracks:
#     (
#         x_section,
#         y_section,
#         z_section,
#         divacancy_section,
#         trivacancy_section,
#         tetravacancy_section,
#         pentavacancy_section,
#         color_origins_section,
#         color_cluster_section,
#         section_count,
#         counts_section,
#         x_centers_section,
#         y_centers_section,
#         z_centers_section,
#         population_section,
#         energy_section,
#         density_ball_section,
#         density_cuboid_section
#     ) = th.generate_random_tracks(
#         int(total_norm_no_tracks - no_inc),
#         PKA,
#         depth_samples,
#         0,
#         1,
#         0,
#         1,
#         x_section,
#         y_section,
#         z_section,
#         color_origins_section,
#         divacancy_section,
#         trivacancy_section,
#         tetravacancy_section,
#         pentavacancy_section,
#         isolated_count_section,
#         color_cluster_section,
#         section_count,
#         counts_section,
#         colors,
#         n,
#         xi,
#         method,
#         x_centers_section,
#         y_centers_section,
#         z_centers_section,
#         population_section,
#         energy_section,
#         density_ball_section,
#         density_cuboid_section,
#         no_inc
#     )
#     ####################ending to add additional tracks####################
#     print("Data prepared for plotting")
#     tph.plot(
#         depth_samples, depth_samples_len, x_section, y_section, z_section, 
#         color_origins_section, counts_section, isolated_count_section, section_count, 
#         x_3D, y_3D, z_3D, color_array, total_count, count_inside_beam, count_outside_beam, 
#         count_outside_beam_upper, count_outside_beam_lower, divacancy_section, trivacancy_section, folder, x_centers_section,
#         y_centers_section, z_centers_section, population_section, energy_section, density_ball_section, density_cuboid_section,"_after_correction")

# # the incident and total_count is the total number of vacancies
# #plt.show()
