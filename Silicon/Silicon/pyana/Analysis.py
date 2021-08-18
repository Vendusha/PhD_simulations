from __future__ import print_function
import pylab as pl
import logging
import numpy as np
#from SimpleHists.plotutils import overlay
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from scipy import stats
from datetime import datetime
from PyAna import *
import PyAnaUtils.bins as bins
import PyAnaUtils.divide as divide
from pathlib import Path
import SimpleHists as sh
from SimpleHists.plotutils import *
from multiprocessing import Process
import Core.FindData
#small routine to readout the dependency of the ratio: percentage of the neutrons that are scattered in the Vanadium foil and detected by the Helium tubes/ number of all the events

x=range(1,20)
y=range(1,20)
Energy_aray = [0_2,0_7,3,6,10,14,18]
Histogram_array = [Recoil_spectra_all_PKA, Recoil_spectra_high_Z, Recoil_spectra_low_Z, Recoil_spetra_SKA, Recoil_spectra_Cascade]
xlabels=["PKA Energy [keV]", "PKA Energy [keV]", "PKA Energy [keV]", "SKA Energy [keV]", "All Cascade energy [keV]"]
ylabels=["Frequency [PKA-1 keV-1]", "Frequency [PKA-1 keV-1]", "Frequency [PKA-1 keV-1]", "SKA Energy [SKA-1 keV-1]", "All Cascade energy [KA-1 keV-1]"]
plt.figure()
for index, key in enumerate Histogram_array:
    for energy in Energy_array:
        IntensityHist = energy+"./"energy+".shist"
        Intensity=sh.HistCollection(IntensityHist)
        Intensity=Intensity.hist(Recoil_spectra)
        Intensity.norm()
        x,y = Intensity.curve()
        plt.plot(x,y,label=str(energy))
    plt.legend()
    plt.xlabel("PKA Energy [keV]")
    plt.ylabel("Frequency [PKA-1 keV-1]")
    plt.yscale("log")
    plt.xscale("log")
    plt.show()
