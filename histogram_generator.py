#!/usr/bin/env python3


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import Core.FindData
import SimpleHists as sh
from scipy.signal import savgol_filter
#########################
# datafile = Core.FindData("SFGenerators","PuBe_data.txt")
# print(datafile)
E2,I2 = np.loadtxt("TriangularSpectra.txt", delimiter=",", unpack=True)
arg_sort = np.argsort(E2)
I2 = I2[arg_sort]
E2 = E2[arg_sort]
x = np.linspace(np.min(E2), np.max(E2), 1500)
y = np.zeros(1500)
y[0] = I2[0]
y[-1] = I2[-1]
for i in range(1, len(x)-1):
    if y[i] < 0:
        print(y[i])
    idx = np.searchsorted(E2, x[i], side='left')
    # if idx == len(E2)-1:
        # idx = idx-1
    # high_idx = np.searchsorted(E2, x[i], side='right')
    y[i] = (I2[idx]-I2[idx-1])/(E2[idx]-E2[idx-1])*(x[i]-E2[idx-1])+I2[idx-1]
    # y[i] = (I2[idx+1]-I2[idx])
# plt.plot(E2, I2, '+')
plt.plot(x, y, 'x')
plt.plot(E2, I2, '+')
plt.show()
E = np.ascontiguousarray(x, dtype=np.float64)
I = np.ascontiguousarray(y, dtype=np.float64)
hc = sh.HistCollection()
h = hc.book1D("Ljubljana_reactor", len(x) ,0 ,max(x),"Ljubljana_reactor")
h.setComment('PuBe energy spectrum created from convolution of multipleGauss functions.')
h.setXLabel('Energy (MeV)')
h.setYLabel('Normalized intensity (a.u.)')
h.fill(E, I)
h.norm()
h.plot()
hc.saveToFile("Ljubljana_reactor",True)
