import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import matplotlib.lines as mlines
import csv

# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.spatial import KDTree
from itertools import islice
import matplotlib.animation as animation
from scipy import stats
import matplotlib.patches as patches
import argparse
import sklearn.cluster
import _pickle as cPickle

DBSCAN = sklearn.cluster.DBSCAN
import math
from matplotlib import gridspec
import optics
import os


def centers_from_borders(borders):
    """From histogram, calculates the centers of everything"""
    return borders[:-1] + np.diff(borders) / 2


def analyze_event(filename, neighbours=[9], r=[4]):
    X = []
    Y = []
    Z = []
    X_end = []
    Y_end = []
    Z_end = []
    Energy = []
    track = []
    track_vacancies = []
    with open(filename) as csvfile:
        for row in csv.reader(csvfile, delimiter=" "):
            # Z_number, x, y, z, energy, x_end, y_end, z_end = int(row[0]), float(row[1]),float(row[2]),float(row[3]),float(row[4]), float(row[6]), float(row[7]), float(row[8])
            Z_number, x, y, z, energy, x_end, y_end, z_end = (
                int(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
                float(row[4]),
                float(row[5]),
                float(row[6]),
                float(row[7]),
            )
            X.append(x * 1000)
            Y.append(y * 1000)
            Z.append(z * 1000)
            X_end.append(x_end * 1000)
            Y_end.append(y_end * 1000)
            Z_end.append(z_end * 1000)
            Energy.append(energy * 1000)
            track.append([x_end * 1000, y_end * 1000, z_end * 1000])
            track_vacancies.append([x * 1000, y * 1000, z * 1000])
    track = np.asarray(track)
    track_vacancies = np.asarray(track_vacancies)
    # r = [0.001, 0.002, 0.003, 0.005]
    # neighbours = [3, 4, 5]
    columns = 1
    no_subfig = np.size(r) * np.size(neighbours)
    # dbsc_helper.dbsc_vendula(no_subfig, track,r,neighbours, columns)
    # threshold_array= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    threshold_array = [4]
    (
        cluster_populations,
        cluster_sizes,
        single_defects,
        defects_in_clusters,
        cluster_populations_optics,
        cluster_sizes_optics,
        no_clusters,
        no_clusters_optics,
        no_clusters_merged,
        cluster_populations_merged,
        cluster_sizes_merged,
        energy_isolated,
        energy_clustered,

    ) = optics.optics_vendula(
        track,
        track_vacancies,
        neighbours,
        threshold_array,
        Energy
    )
    return (
        cluster_populations,
        cluster_sizes,
        single_defects,
        defects_in_clusters,
        cluster_populations_optics,
        cluster_sizes_optics,
        no_clusters,
        no_clusters_optics,
        no_clusters_merged,
        cluster_populations_merged,
        cluster_sizes_merged,
        energy_isolated,
        energy_clustered

    )


n = 20
palette = sns.color_palette(None, n)
parser = argparse.ArgumentParser(description="Visualize events.")
parser.add_argument(
    "folder", metavar="str", type=str, nargs="+", help="folder to load."
)
filen = parser.parse_args()
folder = filen.folder[0]
events = 0
Cluster_populations = np.array([])
Cluster_sizes = np.array([])
Cluster_populations_merged = np.array([])
Cluster_sizes_merged = np.array([])
Single_defects = np.array([])
Ratio = np.array([])
Ratio_no_cl = np.array([])
Defects_in_clusters = np.array([])
Cluster_population_optics = np.array([])
Cluster_sizes_optics = np.array([])
No_clusters_merged = np.array([])
No_clusters = np.array([])
No_clusters_optics = np.array([])
neighbours = [9]
r = [4]
for file in os.listdir(folder):
    if file.endswith(".txt"):
        events += 1
        print(file)
        (
            cluster_populations,
            cluster_sizes,
            single_defects,
            defects_in_clusters,
            cluster_populations_optics,
            cluster_sizes_optics,
            no_clusters,
            no_clusters_optics,
            no_clusters_merged,
            cluster_populations_merged,
            cluster_sizes_merged,
            energy_isolated,
            energy_clustered

        ) = analyze_event((os.path.join(folder, file)), r=r, neighbours=neighbours)
        Cluster_populations = np.append(Cluster_populations, cluster_populations)
        Cluster_populations_merged = np.append(Cluster_populations_merged, cluster_populations_merged)
        Cluster_sizes = np.append(Cluster_sizes, cluster_sizes)
        Cluster_sizes_merged = np.append(Cluster_sizes_merged, cluster_sizes_merged)
        Single_defects = np.append(Single_defects, single_defects)
        if (single_defects>0 or defects_in_clusters>0):
            Ratio = np.append(Ratio, single_defects/(single_defects+defects_in_clusters))
        if (no_clusters>0 or single_defects>0):
            Ratio_no_cl = np.append(Ratio_no_cl, no_clusters/(no_clusters_optics+single_defects))
        Single_defects = np.append(Single_defects, single_defects)
        Defects_in_clusters = np.append(Defects_in_clusters, defects_in_clusters)
        Cluster_population_optics = np.append(
            Cluster_population_optics, cluster_populations_optics
        )
        Cluster_sizes_optics = np.append(Cluster_sizes_optics, cluster_sizes_optics)
        No_clusters = np.append(No_clusters, no_clusters)
        No_clusters_merged = np.append(No_clusters_merged, no_clusters_merged)
        No_clusters_optics = np.append(No_clusters_optics, no_clusters_optics)
# print(Cluster_populations)
energy = 100
fig, axs = plt.subplots(3, 1, figsize=(12, 12))
bins = np.linspace(0, 50, 51)

# Number of clusters optics
y, bin_edges = np.histogram((No_clusters_optics), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_no_clusters_optics = y
x_no_clusters_optics = x
y = np.asarray(y / events)
y_error = np.sqrt(y_no_clusters_optics) / events
axs[0].bar(x, y, yerr=y_error)
# axs[1, 0].set_title('Axis [1, 0]')
axs[0].set_xlabel("# of clusters [-]")
axs[0].set_ylabel("Probability [-]")
axs[0].set_title("Probability of the # of clusters created by " + str(energy) + " keV Si OPTICS")


# Cluster size optics
bins = np.linspace(0, 50, 250)
y, bin_edges = np.histogram((Cluster_sizes_optics), bins=bins)
x = centers_from_borders(bin_edges)
y_cluster_sizes_optics = y
y_error = np.sqrt(y_cluster_sizes_optics)
x_cluster_sizes_optics = x
c = np.sum(y)
if c>0:
    y = np.asarray(y) / c
    y_error = np.sqrt(np.asarray(y_cluster_sizes_optics)) / c
axs[1].step(x, y, color="red")
axs[1].fill_between(x, y - y_error, y + y_error, step="pre", alpha=0.5, color="red")
# axs[1].fill_between(x, y-y_error, y+y_error, step="pre", alpha=0.5, color="red")
axs[1].set_xlabel("Cluster size [nm]")
axs[1].set_ylabel("Probability [-]")
axs[1].set_title("Probability of a size of cluster created by " + str(energy) + " keV Si OPTICS")
#Cluster population optics
bins = np.linspace(0, 50, 51)
y, bin_edges = np.histogram((Cluster_population_optics), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_cluster_populations_optics = y
y_error = np.sqrt(y)
x_cluster_populations_optics = x
c = np.sum(y)
if c>0:
    y = np.asarray(y) / c
    y_error = np.sqrt(y_cluster_populations_optics) / c
axs[2].bar(x, y, yerr=y_error)
axs[2].set_xlabel("Number of Si-recoils in a cluster [-]")
axs[2].set_ylabel("Probability [-]")
axs[2].set_title(
    "Probability of a population of cluster created by "
    + str(energy)
    + " keV Si OPTICS"
)
fig.tight_layout()
label = str(folder) + "_" + str(neighbours[0]) + "_neighbours_" + str(r[0]) + "_eps"
plt.savefig(label + ".png")

bins = np.linspace(0, 100, 101)
fig1, axs1 = plt.subplots(3, 1, figsize=(12, 12))
y, bin_edges = np.histogram((Single_defects), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_single_defects = y
x_single_defects = x
y = np.asarray(y/events)
y_error = np.sqrt(y_single_defects) / events
axs1[0].bar(x, y, yerr=y_error)
# axs[1, 0].set_title('Axis [1, 0]')
axs1[0].set_xlabel("# of clusters [-]")
axs1[0].set_ylabel("Probability [-]")
axs1[0].set_title(
    "Probability of the # of single defects created by " + str(energy) + " keV Si OPTICS"
)

bins = np.linspace(0, 300, 301)
y, bin_edges = np.histogram((Defects_in_clusters), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_cluster_defects = y
x_cluster_defects = x
y = np.asarray(y / events)
y_error = np.sqrt(y_cluster_defects) / events
axs1[1].bar(x, y, yerr=y_error)
# axs[1, 0].set_title('Axis [1, 0]')
axs1[1].set_xlabel("# of clusters [-]")
axs1[1].set_ylabel("Probability [-]")
axs1[1].set_title("Probability of the # of defects in cluster of" + str(energy) + " keV Si OPTICS")

bins = np.linspace(0, 100, 101)
y, bin_edges = np.histogram((Ratio*100), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_ratio = y
x_ratio = x
y = np.asarray(y / events)
# y_error = np.sqrt(y_cluster_defects) / events
axs1[2].bar(x,y)
# axs[1, 0].set_title('Axis [1, 0]')
axs1[2].set_xlabel("Ratio [%]")
axs1[2].set_ylabel("Probability [-]")
axs1[2].set_title("Ratio" + str(energy) + " keV Si OPTICS")

bins = np.linspace(0, 100, 101)
y, bin_edges = np.histogram((Ratio_no_cl*100), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_ratio_no_cl = y
x_ratio_no_cl = x


fig2, axs2 = plt.subplots(3, 1, figsize=(12, 12))
bins = np.linspace(0, 50, 51)

# Number of clusters merged
y, bin_edges = np.histogram((No_clusters_merged), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_no_clusters_merged = y
x_no_clusters_merged = x
y = np.asarray(y / events)
y_error = np.sqrt(y_no_clusters_merged) / events
axs2[0].bar(x, y, yerr=y_error)
# axs[1, 0].set_title('Axis [1, 0]')
axs2[0].set_xlabel("# of clusters [-]")
axs2[0].set_ylabel("Probability [-]")
axs2[0].set_title("Probability of the # of clusters created by " + str(energy) + " keV Si merged")


# Cluster size optics
bins = np.linspace(0, 50, 250)
y, bin_edges = np.histogram((Cluster_sizes_merged), bins=bins)
x = centers_from_borders(bin_edges)
y_cluster_sizes_merged = y
y_error = np.sqrt(y_cluster_sizes_merged)
x_cluster_sizes_merged = x
c = np.sum(y)
if c>0:
    y = np.asarray(y) / c
    y_error = np.sqrt(np.asarray(y_cluster_sizes_merged)) / c
axs2[1].step(x, y, color="red")
axs2[1].fill_between(x, y - y_error, y + y_error, step="pre", alpha=0.5, color="red")
# axs[1].fill_between(x, y-y_error, y+y_error, step="pre", alpha=0.5, color="red")
axs2[1].set_xlabel("Cluster size [nm]")
axs2[1].set_ylabel("Probability [-]")
axs2[1].set_title("Probability of a size of cluster created by " + str(energy) + " keV Si merged")
#Cluster population optics
bins = np.linspace(0, 50, 51)
y, bin_edges = np.histogram((Cluster_populations_merged), bins=bins)
x = centers_from_borders(bin_edges) - 0.5
y_cluster_populations_merged = y
y_error = np.sqrt(y)
x_cluster_populations_merged = x
c = np.sum(y)
if c>0:
    y = np.asarray(y) / c
    y_error = np.sqrt(y_cluster_populations_merged) / c
axs2[2].bar(x, y, yerr=y_error)
axs2[2].set_xlabel("Number of Si-recoils in a cluster [-]")
axs2[2].set_ylabel("Probability [-]")
axs2[2].set_title(
    "Probability of a population of cluster created by "
    + str(energy)
    + " keV Si OPTICS"
)
fig2.tight_layout()
label = str(folder) + "_" + str(neighbours[0]) + "_neighbours_" + str(r[0]) + "_eps"
plt.savefig(label + ".png")
# bins = np.linspace(0, 100, 200)
# fraction = np.asarray((Single_defects*100/ (Single_defects+Defects_in_clusters)))
# print(fraction)
# y, bin_edges = np.histogram((fraction), bins=bins)
# x = centers_from_borders(bin_edges) - 0.5
# y_fraction = y/events
# x_fraction = x
# # y = np.asarray((Single_defects/ (Single_defects+Defects_in_clusters)))
# # y_error = np.sqrt(Single_defects/Defects_in_clusters) / events
# axs1[2].bar(x_fraction, y_fraction)
# # axs[1, 0].set_title('Axis [1, 0]')
# axs1[2].set_xlabel("# of clusters [-]")
# axs1[2].set_ylabel("Probability [-]")
# axs1[2].set_title("Percent of single si-recoils" + str(energy) + " keV Si OPTICS")

# plt.show()


all = [
    x_no_clusters_optics,
    y_no_clusters_optics,
    x_cluster_sizes_optics,
    y_cluster_sizes_optics,
    x_cluster_populations_optics,
    y_cluster_populations_optics,
    x_no_clusters_merged,
    y_no_clusters_merged,
    x_cluster_sizes_merged,
    y_cluster_sizes_merged,
    x_cluster_populations_merged,
    y_cluster_populations_merged,
    x_ratio,
    y_ratio,
    x_ratio_no_cl,
    y_ratio_no_cl,
    Single_defects,
    Defects_in_clusters,
]
cPickle.dump(all, open(label, "wb"))

# y, bin_edges = np.histogram((Cluster_population), bins=bins)
# x = centers_from_borders(bin_edges)-0.5
# axs[2, 1].bar(x, y)
# # axs[0, 0].set_title('Axis [0, 0]')
# axs[2,1].set_xlabel("Number of Si-recoils in a cluster, "+str(events)+" events")
# axs[2,1].set_ylabel("Counts")
