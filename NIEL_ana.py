import os
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D 
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
    Hexa_vacancies = []
    Hepta_vacancies = []
    Octa_vacancies = []
    Nona_vacancies = []
    Deca_vacancies = []
    Above_tri_vacancies = []
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
                    file_path = os.path.join(item_path, file_name)
                    print(file_path)
                    energy =float(item.split('_')[0])
                    Energy.append(energy)
                    #i1, i2, v1, v2, p1, p2 = Srim_output_ana(item_path)
                    # Open the JSON file and load the data
                    with open(file_path, 'r') as json_file:
                        data = json.load(json_file)
                        Cluster_population_counter = data['Cluster_populations_counter']
                        Total_vacancies.append(data['Defects_total'])
                        Single_vacancies.append(data['Defects_isolated'])
                        Divacancies.append(Cluster_population_counter.get("2"))
                        Trivacancies.append(Cluster_population_counter.get("3"))
                        Tetra_vacancies.append(Cluster_population_counter.get("4"))
                        Penta_vacancies.append(Cluster_population_counter.get("5"))
                        Hexa_vacancies.append(Cluster_population_counter.get("6"))
                        Hepta_vacancies.append(Cluster_population_counter.get("7"))
                        Octa_vacancies.append(Cluster_population_counter.get("8"))
                        Nona_vacancies.append(Cluster_population_counter.get("9"))
                        Deca_vacancies.append(Cluster_population_counter.get("10"))
                        Above_tri_vacancies.append(sum(value for key, value in Cluster_population_counter.items() if int(key) > 3))
                        No_inc.append(data['no_inc'])                        
                        if energy == 1:
                            norm = 1/(data['Defects_total']/data['no_inc'])
    return Energy, Single_vacancies, Divacancies, Trivacancies, Above_tri_vacancies, Total_vacancies, No_inc, norm
Energy_electrons, Single_vacancies_electrons, Divacancies_electrons, Trivacancies_electrons, Above_tri_vacancies_electrons, Total_vacancies_electrons, no_inc_electrons, norm = load_data('./e-s_QGSP_BIC_HP')
#Energy_protons, Single_vacancies_protons, Divacancies_protons, Trivacancies_protons, Above_tri_vacancies_protons, Total_vacancies_protons, no_inc_protons, norm = load_data('./p-s_QGSP_BIC_HP')
Energy_neutrons, Single_vacancies_neutrons, Divacancies_neutrons, Trivacancies_neutrons, Above_tri_vacancies_neutrons, Total_vacancies_neutrons, no_inc_neutrons, norm = load_data('./n_QGSP_BIC_HP')
print(norm)
print(Total_vacancies_electrons)
plt.figure()
plt.plot(Energy_electrons, np.asarray(Total_vacancies_electrons)/(no_inc_electrons)*norm, "--x", color = "violet", label="Total vacancies electrons")
plt.plot(Energy_electrons, np.asarray(Single_vacancies_electrons)/(no_inc_electrons)*norm, "--^", color = "violet", label="Single vacancies electrons")
plt.plot(Energy_electrons, np.asarray(Divacancies_electrons)/(no_inc_electrons)*norm, "--o", color = "violet", label="Divacancies electrons")
plt.plot(Energy_electrons, np.asarray(Trivacancies_electrons)/(no_inc_electrons)*norm, "--s", color = "violet", label="Trivacancies electrons")
plt.plot(Energy_electrons, np.asarray(Above_tri_vacancies_electrons)/(no_inc_electrons)*norm, "--d", color = "violet", label="Above tri vacancies electrons")
plt.plot(Energy_neutrons, np.asarray(Total_vacancies_neutrons)/(no_inc_neutrons)*norm, "--x", color = "cyan", label="Total vacancies neutrons")
plt.plot(Energy_neutrons, np.asarray(Total_vacancies_neutrons)/(no_inc_neutrons)*norm-np.asarray(Single_vacancies_neutrons)/(no_inc_neutrons)*norm, "--x", color = "cyan", label="Total vacancies neutrons")
plt.plot(Energy_neutrons, np.asarray(Single_vacancies_neutrons)/(no_inc_neutrons)*norm, "--^", color = "cyan", label="Single vacancies neutrons")
plt.plot(Energy_neutrons, np.asarray(Divacancies_neutrons)/(no_inc_neutrons)*norm, "--o", color = "cyan", label="Divacancies neutrons")
plt.plot(Energy_neutrons, np.asarray(Trivacancies_neutrons)/(no_inc_neutrons)*norm, "--s", color = "cyan", label="Trivacancies neutrons")
plt.plot(Energy_neutrons, np.asarray(Above_tri_vacancies_neutrons)/(no_inc_neutrons)*norm, "--d", color = "cyan", label="Above tri vacancies neutrons")
print(Energy_neutrons)
print(np.asarray(Total_vacancies_neutrons)/(no_inc_neutrons))
x_griffin, y_griffin = np.loadtxt('NIELdata/NIEL_Neutrons_Griffin.dat',
                            skiprows = 1, unpack=True )
x_huhtinen, y_huhtinen = np.loadtxt('NIELdata/NIEL_Protons_Huhtinen.dat',
                            skiprows = 1, unpack=True )
x_summers, y_summers = np.loadtxt('NIELdata/NIEL_Protons_Summers.dat',
                            skiprows = 1, unpack=True )
x_konobeyev, y_konobeyev = np.loadtxt('NIELdata/NIEL_Neutrons_Konobeyev.dat',
                            skiprows = 1, unpack=True )

x_summers_e, y_summers_e = np.loadtxt('NIELdata/NIEL_Electrons_Summers.dat',
                            skiprows = 1, unpack=True )
x_sr_niel_e, y_sr_niel_e, _, _ = np.loadtxt('NIELdata/NIEL_Electrons_sr-niel.dat',
                            skiprows = 1, unpack=True )
x_sr_niel_p, y_sr_niel_p, _, _ = np.loadtxt('NIELdata/NIEL_Protons_sr-niel.dat',
                            skiprows = 1, unpack=True )
x_sr_niel_p_hadron, y_sr_niel_p_hadron, _, _ = np.loadtxt('NIELdata/NIEL_Protons_sr-niel_hadron.dat',
                            skiprows = 1, unpack=True )
norm = 1
plt.plot(x_griffin, y_griffin*norm, color = "blue", label="neutrons Griffin")
plt.plot(x_summers, y_summers*norm, color = "red", label = "protons Summers")
plt.plot(x_huhtinen, y_huhtinen*norm, color = "black", label = "protons Huhtinen")
plt.plot(x_konobeyev, y_konobeyev*norm, color = "green", label = "neutrons konobeyev")
plt.plot(x_summers_e, y_summers_e*norm, color = "pink", label = "electrons Summers")
plt.plot(x_sr_niel_e, y_sr_niel_e/(0.00002144*95), color = "orange", label = "electrons sr-niel") # original value in MeV*cm^2/g, converted to keV*cm^2/g, converted to MeV/mb
plt.plot(x_sr_niel_p, y_sr_niel_p/(0.00002144*95), color = "orange", label = "protons sr-niel")
plt.plot(x_sr_niel_p_hadron, y_sr_niel_p_hadron/(0.00002144*95), color = "orange", label = "protons hadron sr-niel")
plt.xscale('log')
plt.yscale('log')
plt.ylabel("D(E) /95 MeV mb")
plt.title("NIEL compared to reference values.")
plt.xlabel("Energy of the incident particle [MeV]")

legend_elements_1 =[
    Line2D([0], [0], marker='x', color='black',markerfacecolor='black',label='Total vacancies', markersize=10),
    Line2D([0], [0], marker='^', color='black', lw=1, markerfacecolor='black',label='Single vacancies', markersize=10),
    Line2D([0], [0], marker='o', color='black', lw=1, markerfacecolor='black',label='Divacancies', markersize=10),
    Line2D([0], [0], marker='s', color='black', lw=1, markerfacecolor='black',label='Trivacancies', markersize=10),
    Line2D([0], [0], marker='d', color='black', lw=1, markerfacecolor='black',label='Above tri vacancies', markersize=10)]
legend_elements_2 =[
    Line2D([0], [0], color='violet',label='electrons'),
    Line2D([0], [0], color='cyan', label='neutrons')]

legend_1 = plt.legend(handles=legend_elements_2, loc='center right')
legend_2 = plt.legend(handles=legend_elements_1, loc='upper right')
plt.gca().add_artist(legend_1)
#plt.legend()
plt.xlim(0.001,10e5)
plt.ylim(0.0005,230)
plt.show()