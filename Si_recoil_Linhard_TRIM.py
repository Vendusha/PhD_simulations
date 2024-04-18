import os
import json
import numpy as np
import re
import lindhardt as lndh
import plot_niel_trim
import matplotlib.pyplot as plt
import argparse

def log(x, a, b):
    return a + b * np.log(x)


Na = 6.022e23
mSi = 28.086
rhoSi = 2.33
Ntarg = (Na * rhoSi) / mSi  # target density (number per cm3)
rhoSi = 2.33
Nthick = 0.01


def Srim_output_ana(item_path):
    vacancy_path = os.path.join(item_path, "VACANCY.txt")
    ioniz_path = os.path.join(item_path, "IONIZ.txt")
    phonon_path = os.path.join(item_path, "PHONON.txt")
    e2rec_path = os.path.join(item_path, "E2RECOIL.txt")
    range_path = os.path.join(item_path, "RANGE.txt")
    tdata_path = os.path.join(item_path, "TDATA.txt")
    with open(range_path, "r") as f:
        for line in f:
            if "Ion Average Range =" in line:
                range = float(line.split()[4])
                break
    with open(tdata_path, "r") as f:
        for line in f:
            if "Average Vacancy/Ion   =" in line:
                vacancy_per_ion = float(line.split()[3])
                break
    e_binding = 2
    e_init = 500
    depth, vac_ion, vac_recoil = np.loadtxt(
        open(vacancy_path, "rt").readlines()[:-1], skiprows=34, unpack=True
    )
    depth, ioniz_ion, ioniz_recoil = np.loadtxt(ioniz_path, skiprows=25, unpack=True)
    depth, phonon_ion, phonon_recoil = np.loadtxt(phonon_path, skiprows=23, unpack=True)
    depth, e2rec_ion, e2rec_recoil = np.loadtxt(e2rec_path, skiprows=24, unpack=True)
    depth_width = depth[1] - depth[0]
    vac_ion_tot = np.sum(vac_ion) * depth_width * 2
    vac_recoil_tot = np.sum(vac_recoil) * depth_width * 2
    ioniz_ion_tot = np.sum(ioniz_ion) * depth_width
    ioniz_recoil_tot = np.sum(ioniz_recoil) * depth_width
    phonon_ion_tot = np.sum(phonon_ion) * depth_width
    phonon_recoil_tot = np.sum(phonon_recoil) * depth_width
    e2rec_ion_tot = np.sum(e2rec_ion) * depth_width
    e2rec_recoil_tot = np.sum(e2rec_recoil) * depth_width
    tot_loss = (
        vac_ion_tot
        + ioniz_ion_tot
        + phonon_ion_tot
        + vac_recoil_tot
        + ioniz_recoil_tot
        + phonon_recoil_tot
    )
    vac_recoil_tot -= vac_ion_tot
    I1 = ioniz_ion_tot * 100 / tot_loss
    I2 = ioniz_recoil_tot * 100 / tot_loss
    V1 = vac_ion_tot * 100 / tot_loss
    V2 = vac_recoil_tot * 100 / tot_loss
    P1 = phonon_ion_tot * 100 / tot_loss
    P2 = phonon_recoil_tot * 100 / tot_loss
    #print(range)
    return (I1, I2, V1, V2, P1, P2, range, vacancy_per_ion)


def load_data(root_folder):
    Energy = []
    No_inc = []
    Total_vacancies = []
    Single_vacancies = []
    Divacancies = []
    Trivacancies = []
    Tetra_vacancies = []
    Penta_vacancies = []
    vacancies_6_10 = []
    vacancies_11_50 = []
    vacancies_51_100 = []
    vacancies_101_200 = []
    vacancies_201_500 = []
    vacancies_above_500 = []
    Cluster_energies = []
    Energy_TRIM = []
    Energy_primary_recoils = []
    norm = False
    Energy_phonons_ion = []
    Energy_phonons_recoil = []
    Ion_ionization = []
    Ion_vacancy = []
    Ion_phonon = []
    Recoil_ionization = []
    Recoil_vacancy = []
    Recoil_phonon = []
    Single_vacancies_no_ion = []
    Total_vacancies_no_ion = []
    Energy_no_ion = []
    No_inc_no_ion = []
    Energy_Ian = []
    Single_vacancies_Ian = []
    Total_vacancies_Ian = []
    No_inc_Ian = []
    Range = []
    Vacancy_per_ion = []
    Energy_vacancies = []
    ###### example of what all is in the dict
    # dict = {
    #     "Cluster_populations_counter": Cluster_populations_counter,
    #     "Defects_isolated": Defects_isolated,
    #     "Defects_clustered": Defets_clustered,
    #     "No_clusters": No_clusters,
    #     "Defects_total": Defects_total,
    #     "Cluster_populations": Cluster_populations,
    #     "Cluster_volumes_cuboid": Cluster_volumes_cuboid,
    #     "Cluster_lengths_cuboid": Cluster_lengths_cuboid,
    #     "Cluster_radius": Cluster_radius,
    #     "Energy_isolated": Energy_isolated,
    #     "Energy_clustered": Energy_clustered,
    #     "no_inc": no_inc,
    #     "Cluster_energies": Cluster_energies,
    #     "Cluster_densities_sphere": Cluster_densities_sphere,
    #     "Cluster_densities_cuboid": Cluster_densities_cuboid,
    # }
    for item in sorted(
        os.listdir(root_folder),
        key=lambda x: (
            float(re.search(r"\d+\.\d+|\d+", x).group())
            if re.search(r"\d+\.\d+|\d+", x)
            else 0,
            x,
        ),
    ):
        item_path = os.path.join(root_folder, item)
        if item == ".vscode":
            continue
        elif item == "GUI_summary.txt":
            continue
        energy = float(item.split("_")[0])
        print(energy)
        # if energy > 0.005:
        #     break

        if os.path.isdir(item_path):
            # Iterate over all the files in the subfolder
            for file_name in os.listdir(item_path):
                if file_name == "total_count_no_ion.json":
                    file_path = os.path.join(item_path, file_name)
                    with open(file_path, "r") as json_file:
                        data = json.load(json_file)
                        Energy_no_ion.append(energy)
                        Single_vacancies_no_ion.append(data["Defects_isolated"])
                        Total_vacancies_no_ion.append(data["Defects_total"])
                        No_inc_no_ion.append(data["no_inc"])
                if file_name == "total_count_Ians_counting.json":
                    file_path = os.path.join(item_path, file_name)
                    with open(file_path, "r") as json_file:
                        data = json.load(json_file)
                        Energy_Ian.append(energy)
                        Single_vacancies_Ian.append(data["Defects_isolated"])
                        Total_vacancies_Ian.append(data["Defects_total"])
                        No_inc_Ian.append(data["no_inc"])
                if file_name =="COLLISON.txt":
                    Energy.append(energy)
                    i1, i2, v1, v2, p1, p2, range, vacancy_per_ion = Srim_output_ana(item_path)
                    Vacancy_per_ion.append(vacancy_per_ion)
                    Range.append(range)
                    Ion_ionization.append(i1)
                    Ion_vacancy.append(v1)
                    Ion_phonon.append(p1)
                    Recoil_ionization.append(i2)
                    Recoil_vacancy.append(v2)
                    Recoil_phonon.append(p2)
                if file_name == "total_count_ion_included"+str(xi)+".json":
                    file_path = os.path.join(item_path, file_name)
                    print(file_path)
                    # Open the JSON file and load the data
                    with open(file_path, "r") as json_file:
                        data = json.load(json_file)
                        Cluster_population_counter = data["Cluster_populations_counter"]
                        Total_vacancies.append(data["Defects_total"])
                        Single_vacancies.append(data["Defects_isolated"])
                        Energy_vacancies.append(energy)
                        # total_vacncies_check = sum(np.asarray(data['Defects_clustered'])+np.asarray(data['Defects_isolated']))
                        Divacancies.append(
                            np.nan_to_num(Cluster_population_counter.get("2"))
                        )
                        Trivacancies.append(
                            np.nan_to_num(Cluster_population_counter.get("3"))
                        )
                        Tetra_vacancies.append(
                            np.nan_to_num(Cluster_population_counter.get("4"))
                        )
                        Penta_vacancies.append(
                            np.nan_to_num(Cluster_population_counter.get("5"))
                        )
                        vacancies_6_10.append(
                            np.nan_to_num(
                                sum(
                                    value
                                    for key, value in Cluster_population_counter.items()
                                    if int(key) > 6 and int(key) <= 10
                                )
                            )
                        )
                        vacancies_11_50.append(
                            np.nan_to_num(
                                sum(
                                    value
                                    for key, value in Cluster_population_counter.items()
                                    if int(key) > 11 and int(key) <= 20
                                )
                            )
                        )
                        vacancies_51_100.append(
                            np.nan_to_num(
                                sum(
                                    value
                                    for key, value in Cluster_population_counter.items()
                                    if int(key) > 51 and int(key) <= 100
                                )
                            )
                        )
                        vacancies_101_200.append(
                            np.nan_to_num(
                                sum(
                                    value
                                    for key, value in Cluster_population_counter.items()
                                    if int(key) > 101 and int(key) <= 200
                                )
                            )
                        )
                        vacancies_201_500.append(
                            np.nan_to_num(
                                sum(
                                    value
                                    for key, value in Cluster_population_counter.items()
                                    if int(key) > 201 and int(key) <= 500
                                )
                            )
                        )
                        vacancies_above_500.append(
                            np.nan_to_num(
                                sum(
                                    value
                                    for key, value in Cluster_population_counter.items()
                                    if int(key) > 500
                                )
                            )
                        )
                        Energy_primary_recoil = np.concatenate(
                            data["Energy_primary_recoil"]
                        ).ravel()
                        Energy_TRIM.append(
                            (
                                np.sum(data["Energy_clustered"])
                                + np.sum(data["Energy_isolated"])
                            )
                            / 10 ** 6
                        )  # in MeV
                        Energy_primary_recoils.append(
                            np.sum(Energy_primary_recoil) / 10 ** 6
                        )
                        No_inc.append(data["no_inc"])
                        if energy == 1:
                            norm = 1 / (data["Defects_total"] / data["no_inc"])
    Divacancies = [0 if x is None else x for x in Divacancies]
    Trivacancies = [0 if x is None else x for x in Trivacancies]
    Tetra_vacancies = [0 if x is None else x for x in Tetra_vacancies]
    Penta_vacancies = [0 if x is None else x for x in Penta_vacancies]
    vacancies_6_10 = [0 if x is None else x for x in vacancies_6_10]
    vacancies_11_50 = [0 if x is None else x for x in vacancies_11_50]
    vacancies_51_100 = [0 if x is None else x for x in vacancies_51_100]
    vacancies_101_200 = [0 if x is None else x for x in vacancies_101_200]
    vacancies_201_500 = [0 if x is None else x for x in vacancies_201_500]
    vacancies_above_500 = [0 if x is None else x for x in vacancies_above_500]
    return (
        np.asarray(Ion_ionization),
        np.asarray(Ion_vacancy),
        np.asarray(Ion_phonon),
        np.asarray(Recoil_ionization),
        np.asarray(Recoil_vacancy),
        np.asarray(Recoil_phonon),
        np.asarray(Energy),
        np.asarray(Total_vacancies) / np.asarray(No_inc),
        np.asarray(Single_vacancies) / np.asarray(No_inc),
        np.asarray(Total_vacancies_no_ion) / np.asarray(No_inc_no_ion),
        np.asarray(Single_vacancies_no_ion) / np.asarray(No_inc_no_ion),
        np.asarray(Energy_no_ion),
        np.asarray(Energy_Ian),
        np.asarray(Total_vacancies_Ian) / np.asarray(No_inc_Ian),
        np.asarray(Single_vacancies_Ian) / np.asarray(No_inc_Ian),
        np.asarray(Divacancies) / np.asarray(No_inc),
        np.asarray(Trivacancies) / np.asarray(No_inc),
        np.asarray(Tetra_vacancies) / np.asarray(No_inc),
        np.asarray(Penta_vacancies) / np.asarray(No_inc),
        np.asarray(vacancies_6_10) / np.asarray(No_inc),
        np.asarray(vacancies_11_50) / np.asarray(No_inc),
        np.asarray(vacancies_51_100) / np.asarray(No_inc),
        np.asarray(vacancies_101_200) / np.asarray(No_inc),
        np.asarray(vacancies_201_500) / np.asarray(No_inc),
        np.asarray(vacancies_above_500) / np.asarray(No_inc),
        Range,
        Vacancy_per_ion,
        np.asarray(Energy_vacancies),
        np.asarray(No_inc)
    )

parser = argparse.ArgumentParser(
    description="Load all files of particular xi. \
    Requires 1 argument: xi"
)
parser.add_argument("xi", type=str, help="xi value")
args = parser.parse_args()
xi = args.xi
# create a dict for Ion_ionization, Ion_vacancy, Ion_phonon, .....and load the data with
# keys 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 for the different ions (Hydrogen, Helium,... Phosphorus)
# load user input
Ion_ionization = {}
Ion_vacancy = {}
Ion_phonon = {}
Recoil_ionization = {}
Recoil_vacancy = {}
Recoil_phonon = {}
Energy = {}
Total_vacancies = {}
Single_vacancies = {}
Total_vacancies_no_ion = {}
Single_vacancies_no_ion = {}
Energy_no_ion = {}
Energy_Ian = {}
Total_vacancies_Ian = {}
Single_vacancies_Ian = {}
Divacancies = {}
Trivacancies = {}
Tetra_vacancies = {}
Penta_vacancies = {}
vacancies_6_10 = {}
vacancies_11_50 = {}
vacancies_51_100 = {}
vacancies_101_200 = {}
vacancies_201_500 = {}
vacancies_above_500 = {}
Range = {}
Vacancy_per_ion = {}
Energy_vacancies = {}
No_inc = {}
Elements = [
    "Hydrogen",
    "Helium",
    "Lithium",
    "Beryllium",
    "Boron",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Fluorine",
    "Neon",
    "Sodium",
    "Magnesium",
    "Aluminium",
    "Silicon",
    "Phosphorus",
]
# load_data for all the elements
for z in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
    (
        Ion_ionization[z],
        Ion_vacancy[z],
        Ion_phonon[z],
        Recoil_ionization[z],
        Recoil_vacancy[z],
        Recoil_phonon[z],
        Energy[z],
        Total_vacancies[z],
        Single_vacancies[z],
        Total_vacancies_no_ion[z],
        Single_vacancies_no_ion[z],
        Energy_no_ion[z],
        Energy_Ian[z],
        Total_vacancies_Ian[z],
        Single_vacancies_Ian[z],
        Divacancies[z],
        Trivacancies[z],
        Tetra_vacancies[z],
        Penta_vacancies[z],
        vacancies_6_10[z],
        vacancies_11_50[z],
        vacancies_51_100[z],
        vacancies_101_200[z],
        vacancies_201_500[z],
        vacancies_above_500[z],
        Range[z],
        Vacancy_per_ion[z],
        Energy_vacancies[z],
        No_inc[z],
    ) = load_data(Elements[z - 1])
# (
#     Energy_GUI,
#     Ion_ionization_GUI,
#     Recoil_ionization_GUI,
#     Ion_vacancy_GUI,
#     Recoil_vacancy_GUI,
#     Ion_phonon_GUI,
#     Recoil_phonon_GUI,
#     Total_vacancies_GUI,
# ) = np.loadtxt("Silicon/GUI_summary.txt", skiprows=1, unpack=True)
#plt.figure()
# plt.plot(Energy[14]*1000, Single_vacancies[14]/Total_vacancies[14]*100, "x", color="b", label="Vendula Optics TRIM")
# plt.plot(Energy_no_ion[14]*1000, Single_vacancies_no_ion[14]/Total_vacancies_no_ion[14]*100, "x", color="g", label="Vendula TRIM")
# plt.plot(
#     [
#         0.03,
#         0.04,
#         0.05,
#         0.06,
#         0.07,
#         0.08,
#         0.09,
#         0.1,
#         0.2,
#         0.3,
#         0.5,
#         0.7,
#         1.0,
#         2.0,
#         5.0,
#         10.0,
#         20.0,
#         100.0,
#     ],
#     [
#         77.3,
#         63.3,
#         55.2,
#         51.2,
#         45.8,
#         42.4,
#         40.7,
#         37.9,
#         29.4,
#         27.5,
#         24.8,
#         26.9,
#         24.5,
#         26.6,
#         25.7,
#         26.6,
#         27.4,
#         28.4,
#     ],
#     "-x",
#     color="r",
#     label="Ian TRIM",
# )
# plt.plot(
#     Energy_Ian[14] * 1000,
#     Single_vacancies_Ian[14] / Total_vacancies_Ian[14] * 100,
#     "x",
#     color="magenta",
#     label="Vendula TRIM monolayer + Optics",
# )

index = 51
# NIEL_GUI_lattice_6 = Ion_vacancy_GUI_lattice_6 + Recoil_vacancy_GUI_lattice_6 + Ion_phonon_GUI_lattice_6 + Recoil_phonon_GUI_lattice_6
NIEL = {}
NIEL_energy_vac = {}
NIEL_energy = {}
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
    NIEL[i] = (
        np.asarray(Ion_vacancy[i])
        + np.asarray(Recoil_vacancy[i])
        + np.asarray(Ion_phonon[i])
        + np.asarray(Recoil_phonon[i])
    )
    NIEL_energy_vac[i] = np.asarray(Ion_vacancy[i]) + np.asarray(Recoil_vacancy[i])
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
    NIEL_energy[i]= np.multiply(np.asarray(NIEL[i][0:]) / 100, np.asarray(Energy[i]))
    NIEL_energy_vac[i] = np.multiply(np.asarray(NIEL_energy_vac[i][0:]) / 100, np.asarray(Energy[i]))
    np.savetxt(
        "NIEL_" + str(Elements[i-1]) + ".txt",
        np.column_stack(
            (
                np.asarray(Energy[i]) * 1000, # in keV
                np.asarray(NIEL_energy[i]) * 1000, # in keV
                np.asarray(NIEL_energy_vac[i]) * 1000, # in keV
                np.asarray(Range[i])/10000, # in um
                np.asarray(Vacancy_per_ion[i])
            )
        ),
        delimiter="\t",
    )
    matching_indexes = np.where(np.isin(Energy_vacancies[i], Energy[i]))[0]
    #principal debugging
    print("Energy_vacancies")
    print(Energy_vacancies)
    print("Energy")
    print(Energy)
    print("matching_indexes")
    print(matching_indexes)
    np.savetxt(
        "NIEL_OPTICS_" + str(Elements[i-1]) + ".txt",
        np.column_stack(
            (
                np.asarray(Energy[i][matching_indexes]) * 1000, # in keV
                np.asarray(NIEL_energy[i][matching_indexes]) * 1000, # in keV
                np.asarray(NIEL_energy_vac[i][matching_indexes]) * 1000, # in keV
                np.take(Range[i],matching_indexes)/10000, # in um
                np.take(Vacancy_per_ion[i],matching_indexes),
#                np.asarray(No_inc[i][matching_indexes]),
                np.asarray(Single_vacancies[i]),
                np.asarray(Divacancies[i]),
                np.asarray(Trivacancies[i]),
                np.asarray(Tetra_vacancies[i]),
                np.asarray(Penta_vacancies[i]),
                np.asarray(vacancies_6_10[i]),
                np.asarray(vacancies_11_50[i]),
                np.asarray(vacancies_51_100[i]),
                np.asarray(vacancies_101_200[i]),
                np.asarray(vacancies_201_500[i]),
                np.asarray(vacancies_above_500[i]),
                np.asarray(Total_vacancies[i])
            )
        ),
        delimiter="\t",
    )
    # np.savetxt(
    #     "NIEL_OPTICS_" + str(Elements[i-1]) + ".txt",
    #     np.column_stack(
    #         (
    #             np.asarray(Energy[i]) * 1000, # in keV
    #             np.asarray(NIEL_energy[i]) * 1000, # in keV
    #             np.asarray(NIEL_energy_vac[i]) * 1000, # in keV
    #             np.take(Range[i])/10000, # in um
    #             np.take(Vacancy_per_ion[i]),
    #             np.asarray(Single_vacancies[i]),
    #             np.asarray(Divacancies[i]),
    #             np.asarray(Trivacancies[i]),
    #             np.asarray(Tetra_vacancies[i]),
    #             np.asarray(Penta_vacancies[i]),
    #             np.asarray(vacancies_6_10[i]),
    #             np.asarray(vacancies_11_50[i]),
    #             np.asarray(vacancies_51_100[i]),
    #             np.asarray(vacancies_101_200[i]),
    #             np.asarray(vacancies_201_500[i]),
    #             np.asarray(vacancies_above_500[i]),
    #             np.asarray(Total_vacancies[i])
    #         )
    #     ),
    #     delimiter="\t",
    # )
#Energy_Silicon_linhard = np.asarray(Energy[14])
#NIEL_linhard_Si_classic = lndh.niel_si_lindhard_classic(Energy_Silicon_linhard)
# NIEL_linhard[14]_classic = lndh.niel[14]_lindhard_updated(Energy[14]licon_linhard)
# plot_niel_trim.plot_niel_trim(
#     Energy_Silicon_linhard,
#     NIEL_linhard_Si_classic,
#     NIEL_energy,
#     Energy,
#     Total_vacancies,
#     Ion_ionization,
#     Ion_vacancy,
#     Ion_phonon,
#     Recoil_ionization,
#     Recoil_vacancy,
#     Recoil_phonon,
#     Elements,
#     Single_vacancies,
#     Divacancies,
#     Trivacancies,
#     Tetra_vacancies,
#     Penta_vacancies,
#     vacancies_6_10,
#     vacancies_11_50,
#     vacancies_51_100,
#     vacancies_101_200,
#     vacancies_201_500,
#     vacancies_above_500,
#     Energy_vacancies
# )
