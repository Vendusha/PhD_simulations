"""The file to analyze results from the simulation
"""
from __future__ import print_function
from cycler import cycler
import matplotlib.pyplot as plt
import SimpleHists as sh
import numpy as np
from SimpleHists.plotutils import *


def centers_from_borders(borders):
    """From histogram, calculates the centers of everything"""
    return borders[:-1] + np.diff(borders) / 2

particle_type="proton"
hist = sh.HistCollection("PKA_log.shist")
color_array = ["black", "green", "red", "blue"]
# extra_label = [" protons 200 MeV", " protons 20 MeV", " neutrons 200 MeV", " neutrons 20 MeV"]
# plt.figure(0)
# plt.figure(1)
# plt.figure(2)
# for i in range(0, 4):
e_label = ["Elastic", "Inelastic", "Coulomb"]
###########Coulomb##############
if particle_type == "proton":
    no_of_incident = hist.hist("Primary_"+particle_type).integral
Elastic = hist.hist("Recoil_spectra_PKA_Elastic")
Inelastic = hist.hist("Recoil_spectra_PKA_Inelastic")
Coulomb = hist.hist("Recoil_spectra_Coulomb")
plt.figure(0)
histograms = [Elastic, Inelastic, Coulomb]
for i, histogram in enumerate([Elastic, Inelastic, Coulomb]):
    contents, edges = histogram.histogram()
    error_G4 = np.sqrt(contents)/no_of_incident
    x = np.logspace(0, np.max(edges), np.size(contents))
    y = contents/no_of_incident
    plt.step(x/10**9, y, color=color_array[i], label=e_label[i])
    plt.fill_between(x/10**9, y-error_G4, y+error_G4, step="pre", color=color_array[i], alpha = 0.45)
plt.legend()
plt.xlabel("Recoil Energy [MeV]")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Normalized counts of the PKA")
plt.title("Type of interaction")
plt.show()
