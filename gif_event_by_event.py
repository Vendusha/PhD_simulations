import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import matplotlib.lines as mlines
import csv
import plot_h as h
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

def analyze_event(filename,e, neighbours=9, r=4,par=0.1, method="dbscan"):
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
            Z_number, x, y, z, energy, x_end, y_end, z_end = (
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
                float(row[4]),
                float(row[5]),
                float(row[6]),
                float(row[7]),
            )
            X.append(x)
            Y.append(y)
            Z.append(z)
            X_end.append(x_end)
            Y_end.append(y_end)
            Z_end.append(z_end)
            Energy.append(energy)
            track.append([x_end, y_end, z_end])
            track_vacancies.append([x, y, z])
    track = np.asarray(track)
    track_vacancies = np.asarray(track_vacancies)
    Energy = np.asarray(Energy)
    event = optics.Event(filename, e)
    event.clustering(track, track_vacancies, Energy, neighbours, r, par=par, method=method, showPlot=False)
    #You have changed that to track_vacancies!!!!
    # event.clustering(track_vacancies, track, Energy, neighbours, r, showPlot=False)
    return(event)
n = 20
palette = sns.color_palette(None, n)

parser = argparse.ArgumentParser(description="Visualize events.")
parser.add_argument(
    "folder", metavar="str", type=str, nargs="+", help="folder to load."
)

filen = parser.parse_args()
folder = filen.folder[0]
neighbours = int(filen.folder[1])
par = float(filen.folder[2])
method = str(filen.folder[3])
events = 0
e = float(folder)
# neighbours = 9
r = 4
Cluster_populations = np.array([])
Cluster_sizes = np.array([])
Cluster_sizes_ball = np.array([])
Defects_isolated = np.array([])
Defects_clustered = np.array([])
Defects_total = np.array([])
Energy_isolated = np.array([])
Energy_clustered = np.array([])
No_clusters = np.array([])
Energy_init = np.array([])


for file in os.listdir(folder):
    if file.endswith(".txt") and file[0]=="S":
        print(file)
        events += 1
        event = analyze_event((os.path.join(folder, file)),e, r=r, neighbours=neighbours, par=par, method=method)
        Cluster_populations = np.append(Cluster_populations, event.cluster_populations)
        Cluster_sizes = np.append(Cluster_sizes, event.cluster_sizes)
        Cluster_sizes_ball = np.append(Cluster_sizes_ball, event.cluster_sizes_ball)
        No_clusters = np.append(No_clusters, event.no_clusters)
        Energy_isolated = np.append(Energy_isolated, event.energy_isolated)
        Energy_clustered = np.append(Energy_clustered, event.energy_clustered)
        Defects_isolated = np.append(Defects_isolated, event.defects_isolated)
        Defects_clustered = np.append(Defects_clustered, event.defects_clustered)
        Defects_total = np.append(Defects_total, event.defects_total)
        Energy_init = np.append(Energy_init, event.energy)


Events = [events]
all = [Cluster_populations, Cluster_sizes, Cluster_sizes_ball,
      No_clusters, Energy_isolated, Energy_clustered,Defects_clustered, Defects_total, Defects_isolated, Events]
label = str(folder) + "_" + str(neighbours) + "_neighbours_" + str(r) + "_eps_"+str(par)+"_par"
cPickle.dump(all, open(label, "wb"))
