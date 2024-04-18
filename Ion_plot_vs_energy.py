import os
import json
import numpy as np
import re
import lindhardt as lndh
import plot_niel_trim
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"xtick.labelsize": 15})
plt.rcParams.update({"ytick.labelsize": 15})
plt.rcParams.update({"axes.labelsize": 15})
from scipy.interpolate import CubicSpline
Elements = [
    "Hydrogen",
    "Helium",
    "Lithium",
    "Beryllium",
    "Boron",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Fluor",
    "Neon",
    "Sodium",
    "Magnesium",
    "Aluminium",
    "Silicon",
    "Phosphorus",
]
#define two figures
fig, ax1 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()
ax2 = ax1.twinx()
plt.tight_layout()
for idx, element in enumerate(Elements):
    #if idx is odd
    if idx % 2 == 0:
        continue
    Energy, NIEL_energy, NIEL_energy_vac, Range, Vacancy_per_ion, Single_vacancies, Divacancies, Trivacancies, Tetra_vacancies, Penta_vacancies, vacancies_6_10, vacancies_11_50, vacancies_51_100, vacancies_101_200, vacancies_201_500, vacancies_above_500, Total_vacancies= np.loadtxt(
        "NIEL_OPTICS_" + str(element) + ".txt", unpack=True
    )
    # Energy, NIEL_energy, NIEL_energy_vac, Range, Vacancy_per_ion = np.loadtxt(
    #     "NIEL_" + str(element) + ".txt", unpack=True
    # )
    Energy = Energy/1000
    if element == "Sodium" or element =="Neon":
        x = np.linspace(Energy[:-1].min(), Energy[:-1].max(), 10000)
    else:
        x = np.linspace(Energy.min(), Energy.max(), 10000)  
    #plt.figure(fig1.number)
    line=ax1.plot(Energy, NIEL_energy, "x",label=str(element))
    color = line[-1].get_color()
    if element == "Sodium" or element=="Neon":
        f = interp1d(Energy[:-1], NIEL_energy[:-1], kind="slinear")
        y = f(x)
        ax1.plot(x, y, "--", color=color)
        plt.tight_layout()
    else:
        f = interp1d(Energy, NIEL_energy, kind="slinear")
        y = f(x)
        ax1.plot(x, y, "--", color=color)
    #plt.figure(fig2.number)
    line=ax2.plot(Energy, NIEL_energy_vac, "o",label=str(element)+"")
    color = line[-1].get_color()
    if element == "Sodium" or element=="Neon":
        f = interp1d(Energy[:-1], NIEL_energy_vac[:-1], kind="slinear")
        y = f(x)
        ax2.plot(x, y, "--", color=color)
    else:
        f = interp1d(Energy, NIEL_energy_vac, kind="slinear")
        y = f(x)
        ax2.plot(x, y, "--", color=color)
    line = ax3.plot(Energy, NIEL_energy_vac/NIEL_energy*100, "o",label=str(element))
    color = line[-1].get_color()
    y_smooth = savgol_filter(NIEL_energy_vac/NIEL_energy*100, window_length=11, polyorder=2)
    ax3.plot(Energy, y_smooth, "--", color=color)
    print(idx)
    if idx in {1,5,9,13}: #TOFIX
        print(idx)
        line=ax4.plot(Energy, Single_vacancies/Total_vacancies, "x")
        color = line[-1].get_color()
        ax4.plot(Energy, (Total_vacancies-Single_vacancies)/Total_vacancies, "o",color=color, label=str(element))
        y_smooth = savgol_filter(Single_vacancies/Total_vacancies, window_length=11, polyorder=2)
        ax4.plot(Energy, y_smooth, "--", color=color)
        y_smooth = savgol_filter((Total_vacancies-Single_vacancies)/Total_vacancies, window_length=11, polyorder=2)
        ax4.plot(Energy, y_smooth, "--", color=color)
        plt.tight_layout()
ax3.set_xscale("log")
ax3.set_ylabel("NIEL vac/NIEL [%]")
ax3.set_xlabel("Energy of the PKA [MeV]")
ax3.legend(loc="best")
plt.tight_layout()
ax1.set_ylabel("Energy deposited in NIEL [keV]")
ax2.set_ylabel("Energy deposited in NIEL vac [keV]")
ax1.set_xlabel("Energy of the PKA [MeV]")
ax1.set_xscale("log")
ax1.set_yscale("log")
plt.tight_layout()
ax2.set_yscale("log")
plt.tight_layout()
lines1, labels1 = ax2.get_legend_handles_labels()
ax1.legend(lines1, labels1, loc="best")
handles = [plt.Line2D([], [], marker='x', linestyle='None', color='black'),
           plt.Line2D([], [], marker='o', linestyle='None', color='black')]
labels = ['NIEL', 'NIEL vac']
ax2.legend(handles, labels, loc="lower right")
plt.tight_layout()
plt.savefig("NIEL_vs_energy.png",dpi=900)
ax4.set_xlabel("Energy of the PKA [MeV]")
ax4.set_ylabel("Fraction of defects")
ax4.set_xscale("log")
handles_4 = [plt.Line2D([], [], marker='x', linestyle='None', color='black'),
           plt.Line2D([], [], marker='o', linestyle='None', color='black')]
labels_4 = ['Isolated defects', 'Clustered defects']
legend=ax4.legend(loc="best")
legend_4 = ax4.legend(handles_4, labels_4, loc = 'lower right')
ax4.add_artist(legend_4)
ax4.add_artist(legend)
# plt.figure(fig1.number)
# plt.xlabel("Energy of the PKA [MeV]")
# plt.ylabel("Energy deposited in NIEL [keV]")
# plt.legend()
# plt.figure(fig2.number)
# plt.xlabel("Energy of the PKA [MeV]")
# plt.ylabel("Energy deposited in NIEL vac [keV]")
# plt.legend()
#plt.xscale("log")
plt.show()
