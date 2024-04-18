import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import matplotlib.lines as mlines
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
from matplotlib import rcParams

rcParams["text.usetex"] = True

def centers_from_borders(borders):
    """From histogram, calculates the centers of everything"""
    return borders[:-1] + np.diff(borders) / 2

#toPlot
# Number of clusters distribution, cluster size distribution, cluster population distribution
# as above but merged
# Energy_isolated Energy_clustered Energy_ratio
# Number of clusters vs isolated events
# Everything above per keV

def plot_isolated(defects_isolated, events, energy, showPlot = False, ax="None"):
    bins = np.linspace(0, 600, 601)
    y_orig, bin_edges = np.histogram((defects_isolated), bins=bins)
    x = centers_from_borders(bin_edges)
    y_error_orig = np.sqrt(y_orig)
    y_error = np.sqrt(y_orig) / events
    y = np.asarray(y_orig/events)
    if showPlot:
        ax.bar(x, y, yerr=y_error, label=str(energy)+" keV")
        ax.set_xlabel("\# of Defects$_{isolated}$ [-]")
        ax.set_ylabel("Probability [-]")
        ax.set_title("Probability of Si event creating a particular number of isolated defects.")
    return(x, y_orig, y_error_orig)

def plot_clustered(defects_clustered, events, energy, showPlot = False, ax="None"):
    bins = np.linspace(0, 600, 601)
    y_orig, bin_edges = np.histogram((defects_clustered), bins=bins)
    x = centers_from_borders(bin_edges)
    y_error_orig = np.sqrt(y_orig)
    y_error = np.sqrt(y_orig) / events
    y = np.asarray(y_orig/events)
    if showPlot:
        ax.bar(x, y, yerr=y_error, label=str(energy)+" keV")
        ax.set_xlabel("\# of Defects$_{clustered}$ [-]")
        ax.set_ylabel("Probability [-]")
        ax.set_title("Probability of Si event creating a particular number of clustered defects.")
    return(x, y_orig, y_error_orig)


def plot_ratio_single_clustered (defects_isolated, defects_clustered, events, energy, showPlot = False, ax="None"):
    bins = np.linspace(0, 100, 101)
    y, bin_edges = np.histogram((defects_isolated/(defects_isolated+defects_clustered)*100), bins=bins)
    x = centers_from_borders(bin_edges)
    if showPlot:
        ax.bar(x, y, label=str(energy)+" keV")
        ax.set_xlabel("\# of Defects$_{clustered}$ [-]")
        ax.set_ylabel("Probability [-]")
        ax.set_title("Probability of Si event creating a particular number of clustered defects.")
    return(x, y)

# Energy

def plot_energy_isolated(energy_isolated, events, energy, showPlot = False, ax="None"):
    bins = np.linspace(0, 600, 601)
    y_orig, bin_edges = np.histogram((defects_isolated), bins=bins)
    x = centers_from_borders(bin_edges)
    y_error_orig = np.sqrt(y_orig)
    y_error = np.sqrt(y_orig) / events
    y = np.asarray(y_orig/events)
    if showPlot:
        ax.bar(x, y, yerr=y_error, label=str(energy)+" keV")
        ax.set_xlabel("\# of Defects$_{isolated}$ [-]")
        ax.set_ylabel("Probability [-]")
        ax.set_title("Probability of Si event creating a particular number of isolated defects.")
    return(x, y_orig, y_error_orig)

def plot_clustered(defects_clustered, events, energy, showPlot = False, ax="None"):
    bins = np.linspace(0, 600, 601)
    y_orig, bin_edges = np.histogram((defects_clustered), bins=bins)
    x = centers_from_borders(bin_edges)
    y_error_orig = np.sqrt(y_orig)
    y_error = np.sqrt(y_orig) / events
    y = np.asarray(y_orig/events)
    if showPlot:
        ax.bar(x, y, yerr=y_error, label=str(energy)+" keV")
        ax.set_xlabel("\# of Defects$_{clustered}$ [-]")
        ax.set_ylabel("Probability [-]")
        ax.set_title("Probability of Si event creating a particular number of clustered defects.")
    return(x, y_orig, y_error_orig)


def plot_ratio_single_clustered (defects_isolated, defects_clustered, events, energy, showPlot = False, ax="None"):
    bins = np.linspace(0, 100, 101)
    y, bin_edges = np.histogram((defects_isolated/(defects_isolated+defects_clustered)*100), bins=bins)
    x = centers_from_borders(bin_edges)
    if showPlot:
        ax.bar(x, y, label=str(energy)+" keV")
        ax.set_xlabel("\# of Defects$_{clustered}$ [-]")
        ax.set_ylabel("Probability [-]")
        ax.set_title("Probability of Si event creating a particular number of clustered defects.")
    return(x, y)
# show = True
# fig0, ax0 = plt.subplots(3, 1, figsize=(12, 12))
# e=Energy_init[0]
# x_isolated, y_isolated, y_err_isolated = h.plot_isolated(Defects_isolated, events, e, showPlot= show, ax=ax0[0])
# x_clustered, y_clustered, y_err_clustered = h.plot_clustered(Defects_clustered, events, e, showPlot = show, ax=ax0[1])
# x_clustered, y_clustered = h.plot_ratio_single_clustered(Defects_isolated, Defects_clustered, events, e, showPlot = show, ax=ax0[2])



# x_isolated, y_isolated, y_err_isolated = h.plot_isolated(Defects_isolated, events, e)
# x_clustered, y_clustered, y_err_clustered = h.plot_clustered(Defects_clustered, events, e)



# analyze_event("./Silicon/trial/Si_recoil_0.100000_MeV_cutoff_0.021_keV_1.txt")
plt.show()
