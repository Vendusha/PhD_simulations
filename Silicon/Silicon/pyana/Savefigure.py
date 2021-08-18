from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import pylab as pl
import logging
import numpy as np
from SimpleHists.plotutils import overlay
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
Detector=[]
for i in range(36):
    DetectorString="Detector"+str(i+1)
    DetectorStringHist="Detector"+str(i+1)+".shist"
    Detector.append(sh.HistCollection(DetectorStringHist))
    h_time_of_flight=Detector[i].hist("Time_of_flight_micros")
    x,y = h_time_of_flight.curve()   
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure()
    plt.xlabel(r'Time of flight [\mu s]',fontsize=16)
    plt.ylabel(r'Counts [-]',fontsize=16)
    plt.plot(x,y)
    plt.savefig(DetectorString+'TOF.png')
    detector_spectra_neutrons=Detector[i].hist("Detector_spectra_neutrons")
    x,y =detector_spectra_neutrons.curve()
    plt.figure()
    plt.plot(x,y)
    plt.xlabel(r'Detector neutrons [meV]',fontsize=16)
    plt.ylabel(r'Counts [-]',fontsize=16)
    plt.savefig(DetectorString+'neutrons.png')
    detector_spectra_gammas=Detector[i].hist("Detector_spectra_gammas")
    x,y =detector_spectra_gammas.curve()
    plt.figure()
    plt.plot(x,y)
    plt.yscale('log')
    plt.xlabel(r'Detector gammas [MeV]',fontsize=16)
    plt.ylabel(r'Counts [-]',fontsize=16)
    plt.savefig(DetectorString+'gammas.png')
