import os
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
import lindhardt as lnd  
plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"xtick.labelsize": 15})
plt.rcParams.update({"ytick.labelsize": 15})
plt.rcParams.update({"axes.labelsize": 15})
def log(x, a, b):
    return a + b * np.log(x)
Na    = 6.022e23
mSi   = 28.086
rhoSi = 2.33
Ntarg = (Na*rhoSi)/mSi  # target density (number per cm3)
rhoSi = 2.33
Nthick = 0.01
def Srim_output_ana(item_path):
    vacancy_path = os.path.join(item_path, "VACANCY.txt")
    ioniz_path = os.path.join(item_path, "IONIZ.txt")
    phonon_path = os.path.join(item_path, "PHONON.txt")
    e2rec_path = os.path.join(item_path, "E2RECOIL.txt")
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
    vac_recoil_tot-=vac_ion_tot
    I1 = ioniz_ion_tot * 100 / tot_loss
    I2 = ioniz_recoil_tot * 100 / tot_loss
    V1 = vac_ion_tot * 100 / tot_loss
    V2 = vac_recoil_tot * 100 / tot_loss
    P1 = phonon_ion_tot * 100 / tot_loss
    P2 = phonon_recoil_tot * 100 / tot_loss
    return (I1, I2, V1, V2, P1, P2)
def load_data(root_folder):
    Energy = []
    No_inc = []
    Total_vacancies = []
    Single_vacancies = []
    Divacancies = []
    Trivacancies = []
    Tetra_vacancies = []
    Penta_vacancies = []
    vacancies_6_10 =[]
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
    norm = False
    ###### example of what all is in the dict
    #dict ={"Cluster_populations_counter":Cluster_populations_counter,"Defects_isolated":Defects_isolated, 
    #"Defects_clustered": Defets_clustered, "No_clusters":No_clusters, "Defects_total":Defects_total, 
    #"Cluster_populations":Cluster_populations, "Cluster_volumes_rectangle":Cluster_volumes_rectangle,
    #"Cluster_lengths_rectangle":Cluster_lengths_rectangle, "Cluster_radius":Cluster_radius,
    #"Energy_isolated":Energy_isolated, "Energy_clustered":Energy_clustered, "no_inc":no_inc}
    # Iterate over all the items in the root folder
    for item in sorted(os.listdir(root_folder), key=lambda x: (float(re.search(r'\d+\.\d+|\d+', x).group()) if re.search(r'\d+\.\d+|\d+', x) else 0, x)):
        item_path = os.path.join(root_folder, item)
        if item == ".vscode":
            continue   
        if os.path.isdir(item_path):
            # Iterate over all the files in the subfolder
            for file_name in os.listdir(item_path):
                if file_name.endswith('.json'):
                    energy = float(item.split("_")[0])
                    file_path = os.path.join(item_path, file_name)
                    Energy.append(energy)
                    #i1, i2, v1, v2, p1, p2 = Srim_output_ana(item_path)
                    i1, i2, v1, v2, p1, p2 =0,0,0,0,0,0
                    Ion_ionization.append(i1)
                    Ion_vacancy.append(v1)
                    Ion_phonon.append(p1)
                    Recoil_ionization.append(i2)
                    Recoil_vacancy.append(v2)
                    Recoil_phonon.append(p2)
                    file_path = os.path.join(item_path, file_name)
                    print(file_path)
                    # Open the JSON file and load the data
                    with open(file_path, "r") as json_file:
                        data = json.load(json_file)
                        Cluster_population_counter = data["Cluster_populations_counter"]
                        Total_vacancies.append(data["Defects_total"])
                        Single_vacancies.append(data["Defects_isolated"])
                        # total_vacncies_check = sum(np.asarray(data['Defects_clustered'])+np.asarray(data['Defects_isolated']))
                        Divacancies.append(np.nan_to_num(Cluster_population_counter.get("2")))
                        Trivacancies.append(np.nan_to_num(Cluster_population_counter.get("3")))
                        Tetra_vacancies.append(np.nan_to_num(Cluster_population_counter.get("4")))
                        Penta_vacancies.append(np.nan_to_num(Cluster_population_counter.get("5")))
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
                        #Energy_primary_recoil = np.concatenate(
                        #    data["Energy_primary_recoil"]
                        #).ravel()
                        Energy_TRIM.append(
                            (
                                np.sum(data["Energy_clustered"])
                                + np.sum(data["Energy_isolated"])
                            )
                            / 10 ** 6
                        )  # in MeV
                        #Energy_primary_recoils.append(
                        #    np.sum(Energy_primary_recoil) / 10 ** 6
                        #)
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
        np.asarray(Total_vacancies)/np.asarray(No_inc),
        np.asarray(Single_vacancies)/np.asarray(No_inc),
        np.asarray(Divacancies)/np.asarray(No_inc),
        np.asarray(Trivacancies)/np.asarray(No_inc),
        np.asarray(Tetra_vacancies)/np.asarray(No_inc),
        np.asarray(Penta_vacancies)/np.asarray(No_inc),
        np.asarray(vacancies_6_10)/np.asarray(No_inc),
        np.asarray(vacancies_11_50)/np.asarray(No_inc),
        np.asarray(vacancies_51_100)/np.asarray(No_inc),
        np.asarray(vacancies_101_200)/np.asarray(No_inc),
        np.asarray(vacancies_201_500)/np.asarray(No_inc),
        np.asarray(vacancies_above_500)/np.asarray(No_inc),
        norm
    )

#Ion_ionization_electrons, Ion_vacancy_electrons, Ion_phonon_electrons, Recoil_ionization_electrons, Recoil_vacancy_electrons, Recoil_phonon_electrons, Energy_electrons, Total_vacancies_electrons, Single_vacancies_electrons, Divacancies_electrons, Trivacancies_electrons, Tetra_vacancies_electrons, Penta_vacancies_electrons, vacancies_6_10_electrons, vacancies_11_50_electrons, vacancies_51_100_electrons, vacancies_101_200_electrons, vacancies_201_500_electrons, vacancies_above_500_electrons, norm = load_data('./e-s_QGSP_BIC_HP')
#Ion_ionization_protons, Ion_vacancy_protons, Ion_phonon_protons, Recoil_ionization_protons, Recoil_vacancy_protons, Recoil_phonon_protons, Energy_protons, Total_vacancies_protons, Single_vacancies_protons, Divacancies_protons, Trivacancies_protons, Tetra_vacancies_protons, Penta_vacancies_protons, vacancies_6_10_protons, vacancies_11_50_protons, vacancies_51_100_protons, vacancies_101_200_protons, vacancies_201_500_protons, vacancies_above_500_protons, norm = load_data('./protons')
#Ion_ionization_neutrons, Ion_vacancy_neutrons, Ion_phonon_neutrons, Recoil_ionization_neutrons, Recoil_vacancy_neutrons, Recoil_phonon_neutrons, Energy_neutrons, Total_vacancies_neutrons, Single_vacancies_neutrons, Divacancies_neutrons, Trivacancies_neutrons, Tetra_vacancies_neutrons, Penta_vacancies_neutrons, vacancies_6_10_neutrons, vacancies_11_50_neutrons, vacancies_51_100_neutrons, vacancies_101_200_neutrons, vacancies_201_500_neutrons, vacancies_above_500_neutrons, norm = load_data('./n_QGSP_BIC_HP')
plt.figure()
#norm = 1
#plt.plot(Energy_electrons, np.asarray(Total_vacancies_electrons)*norm, "--x", color = "violet", label="Total vacancies electrons")
#plt.plot(Energy_electrons, np.asarray(Single_vacancies_electrons)*norm, "--^", color = "violet", label="Single vacancies electrons")
#plt.plot(Energy_electrons, np.asarray(Divacancies_electrons)*norm, "--o", color = "violet", label="Divacancies electrons")
#plt.plot(Energy_electrons, np.asarray(Trivacancies_electrons)*norm, "--s", color = "violet", label="Trivacancies electrons")
#plt.plot(Energy_neutrons, np.asarray(Total_vacancies_neutrons)*norm, "--x", color = "cyan", label="Total vacancies neutrons")
#plt.plot(Energy_neutrons, np.asarray(Single_vacancies_neutrons)*norm, "--^", color = "cyan", label="Single vacancies neutrons")
#plt.plot(Energy_neutrons, np.asarray(Divacancies_neutrons)*norm, "--o", color = "cyan", label="Divacancies neutrons")
#plt.plot(Energy_neutrons, np.asarray(Trivacancies_neutrons)*norm, "--s", color = "cyan", label="Trivacancies neutrons")
#print(Energy_neutrons)
#print(np.asarray(Total_vacancies_neutrons))
x_griffin, y_griffin = np.loadtxt('NIELdata/NIEL_Neutrons_Griffin.dat',
                            skiprows = 1, unpack=True )
x_huhtinen, y_huhtinen = np.loadtxt('NIELdata/NIEL_Protons_Huhtinen.dat',
                            skiprows = 1, unpack=True )
x_summers, y_summers = np.loadtxt('NIELdata/NIEL_Protons_Summers.dat',
                            skiprows = 1, unpack=True )
x_konobeyev, y_konobeyev = np.loadtxt('NIELdata/NIEL_Neutrons_Konobeyev.dat',
                            skiprows = 1, unpack=True )
x_griffin_orig, y_griffin_orig = np.loadtxt('NIELdata/Si-28_Griffin_origin_endf.txt',
                            skiprows = 1, unpack=True, delimiter=",")
x_endf_all, y_endf_all = np.loadtxt('NIELdata/endf_neutrons_all.txt',
                            skiprows = 1, unpack=True, delimiter=",")
x_summers_e, y_summers_e = np.loadtxt('NIELdata/NIEL_Electrons_Summers.dat',
                            skiprows = 1, unpack=True )
x_sr_niel_e, y_sr_niel_e, _, _ = np.loadtxt('NIELdata/NIEL_Electrons_sr-niel.dat',
                            skiprows = 1, unpack=True )
x_sr_niel_p, y_sr_niel_p, _, _ = np.loadtxt('NIELdata/NIEL_Protons_sr-niel.dat',
                            skiprows = 1, unpack=True )
x_sr_niel_p_hadron, y_sr_niel_p_hadron, _, _ = np.loadtxt('NIELdata/NIEL_Protons_sr-niel_hadron.dat',
                            skiprows = 1, unpack=True )
norm = 1
plt.plot(x_griffin, y_griffin*norm, color = "blue", label="RD48 collected: priv. comm. Griffin")
#plt.plot(x_griffin_orig/10**6, y_griffin_orig*norm, color = "cyan", label="neutrons Griffin orig")
y_lind = lnd.niel_si_lindhard_classic(x_griffin_orig/10**6)*y_griffin_orig
y_lind_all = lnd.niel_si_lindhard_classic(x_endf_all/10**6)*y_endf_all
#plt.plot(x_griffin_orig/10**6, y_lind*norm, color = "green", label="Lindhard applied on ENDF VII. elastic cross section")
#plt.plot(x_endf_all/10**6, y_lind_all*norm, color = "red", label="Lindhard applied on ENDF VII. total cross section")
#plt.title("Lower energy range neutrons")
# plt.plot(x_summers, y_summers*norm, color = "red", label = "protons Summers")
# plt.plot(x_huhtinen, y_huhtinen*norm, color = "black", label = "protons Huhtinen")
# plt.plot(x_konobeyev, y_konobeyev*norm, color = "green", label = "neutrons konobeyev")
# plt.plot(x_summers_e, y_summers_e*norm, color = "pink", label = "electrons Summers")
#plt.plot(x_sr_niel_e, y_sr_niel_e/(0.00002144*95), color = "orange", label = "electrons sr-niel") # original value in MeV*cm^2/g, converted to keV*cm^2/g, converted to MeV/mb
#plt.plot(x_sr_niel_p, y_sr_niel_p/(0.00002144*95), color = "orange", label = "protons sr-niel")
#plt.plot(x_sr_niel_p_hadron, y_sr_niel_p_hadron/(0.00002144*95), color = "orange", label = "protons hadron sr-niel")
plt.xscale('log')
plt.yscale('log')
plt.ylabel("D(E) /95 MeV mb")
#plt.title("NIEL compared to reference values.")
plt.xlabel("Energy of the incident particle [MeV]")
plt.tight_layout()

# legend_elements_1 =[
#     Line2D([0], [0], marker='x', color='black',markerfacecolor='black',label='Total vacancies', markersize=10),
#     Line2D([0], [0], marker='^', color='black', lw=1, markerfacecolor='black',label='Single vacancies', markersize=10),
#     Line2D([0], [0], marker='o', color='black', lw=1, markerfacecolor='black',label='Divacancies', markersize=10),
#     Line2D([0], [0], marker='s', color='black', lw=1, markerfacecolor='black',label='Trivacancies', markersize=10),
#     Line2D([0], [0], marker='d', color='black', lw=1, markerfacecolor='black',label='Above tri vacancies', markersize=10)]
# legend_elements_2 =[
#     Line2D([0], [0], color='violet',label='electrons'),
#     Line2D([0], [0], color='cyan', label='neutrons')]

# legend_1 = plt.legend(handles=legend_elements_2, loc='center right')
# legend_2 = plt.legend(handles=legend_elements_1, loc='upper right')
# plt.gca().add_artist(legend_1)
plt.legend()
#plt.xlim(0.001, 20)
plt.xlim(0.001,10e5)
plt.ylim(0.0005,230)
plt.show()