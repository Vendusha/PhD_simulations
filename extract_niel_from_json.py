import os
import json
import numpy as np
import re
import lindhardt as lndh
import plot_niel_trim
import matplotlib.pyplot as plt
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Load all json files in particular folder structure.")
    parser.add_argument("xi", type=str, help="xi value of the dbscan parameter")
    return parser.parse_args()

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

def sum_vacancies(counter, min_key, max_key):
    return np.nan_to_num(
        sum(value for key, value in counter.items() if min_key < int(key) <= max_key)
    )
def load_json_file(file_path, data_TRIM):
    # Open the JSON file and load the data
    with open(file_path, "r") as json_file:
        data = json.load(json_file)
        Cluster_population_counter = data["Cluster_populations_counter"]
        data_TRIM["Total_vacancies"].append(data["Defects_total"])
        data_TRIM["Single_vacancies"].append(data["Defects_isolated"])
        # data_TRIM["Cluster_lengths_cuboid"].append(np.mean(data["Cluster_lengths_cuboid"]))
        # data_TRIM["Cluster_volumes_cuboid"].append(np.mean(data["Cluster_volumes_cuboid"]))
        # data_TRIM["Cluster_radius"].append(np.mean(data["Cluster_radius"]))
        Cluster_population = data["Cluster_populations"]
        Cluster_diameter = np.array([value * 2000 for value in data["Cluster_radius"]])  # in nm
        Cluster_length_cuboid = np.array([value * 1000 for value in data["Cluster_lengths_cuboid"]])  # in nm
        Cluster_volume_cuboid = np.array(data["Cluster_volumes_cuboid"])
        for i in range(len(Cluster_diameter)):
           if Cluster_diameter[i] < 0.005:
                Cluster_diameter[i] = 0.47
                Cluster_volume_cuboid[i]=0.47**3
                Cluster_length_cuboid[i] = 0.47
        Cluster_desnity_sphere = (data["Cluster_populations"]) / ((4/3) * np.pi * (Cluster_diameter*1/2) ** 3)
        Cluster_density_cuboid = (data["Cluster_populations"]) / (Cluster_volume_cuboid)
        Cluster_population_per_nm=((data["Cluster_populations"]) / Cluster_diameter)
        data_TRIM["Cluster_lengths_cuboid"].append(np.mean(Cluster_length_cuboid))
        data_TRIM["Cluster_volumes_cuboid"].append(np.mean(Cluster_volume_cuboid))
        data_TRIM["Cluster_radius"].append(np.mean(Cluster_diameter)*1/2)
        data_TRIM["Cluster_densities_sphere"].append(np.mean(Cluster_desnity_sphere))
        data_TRIM["Cluster_densities_cuboid"].append(np.mean(Cluster_density_cuboid))
        data_TRIM["Cluster_population_per_nm"].append(np.mean(Cluster_population_per_nm))
        data_TRIM["Divacancies"].append(np.nan_to_num(Cluster_population_counter.get("2")))
        data_TRIM["Trivacancies"].append(np.nan_to_num(Cluster_population_counter.get("3")))
        data_TRIM["Tetra_vacancies"].append(np.nan_to_num(Cluster_population_counter.get("4")))
        data_TRIM["Penta_vacancies"].append(np.nan_to_num(Cluster_population_counter.get("5")))
        data_TRIM["vacancies_6_10"].append(sum_vacancies(Cluster_population_counter, 6, 10))
        data_TRIM["vacancies_11_50"].append(sum_vacancies(Cluster_population_counter, 11, 50))
        data_TRIM["vacancies_51_100"].append(sum_vacancies(Cluster_population_counter, 51, 100))
        data_TRIM["vacancies_101_200"].append(sum_vacancies(Cluster_population_counter, 101, 200))
        data_TRIM["vacancies_201_500"].append(sum_vacancies(Cluster_population_counter, 201, 500))
        data_TRIM["vacancies_above_500"].append(sum_vacancies(Cluster_population_counter, 501, float('inf')))
        data_TRIM["Energy_TRIM"].append((np.sum(data["Energy_clustered"]) + np.sum(data["Energy_isolated"]))/ 10 ** 6)  # in MeV
        Energy_primary_recoil = np.concatenate(data["Energy_primary_recoil"]).ravel()
        data_TRIM["Energy_primary_recoils"].append(np.sum(Energy_primary_recoil) / 10 ** 6)
        data_TRIM["No_inc"].append(data["no_inc"])
    return data_TRIM

def load_data(root_folder, data_keys):
    data_TRIM = {key: [] for key in data_keys}
    vacancy_keys = ["Total_vacancies", "Single_vacancies", "Divacancies", "Trivacancies", "Tetra_vacancies", "Penta_vacancies", 
                    "vacancies_6_10", "vacancies_11_50", "vacancies_51_100", "vacancies_101_200", "vacancies_201_500", 
                    "vacancies_above_500"]
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
        data_TRIM["Energy_TRIM_sim"].append(energy)
        if os.path.isdir(item_path):
            # Iterate over all the files in the subfolder
            for file_name in os.listdir(item_path):
                if file_name =="COLLISON.txt":
                    i1, i2, v1, v2, p1, p2, range, vacancy_per_ion = Srim_output_ana(item_path)
                    data_TRIM["Ion_ionization"].append(i1)
                    data_TRIM["Ion_vacancy"].append(v1)
                    data_TRIM["Ion_phonon"].append(p1)
                    data_TRIM["Recoil_ionization"].append(i2)
                    data_TRIM["Recoil_vacancy"].append(v2)
                    data_TRIM["Recoil_phonon"].append(p2)
                    data_TRIM["Range"].append(range)
                    data_TRIM["Vacancy_per_ion"].append(vacancy_per_ion)
                if file_name == "total_count_ion_included"+str(xi)+".json":
                    data_TRIM["Energy"].append(energy)
                    file_path = os.path.join(item_path, file_name)
                    print(file_path)
                    data_TRIM = load_json_file(file_path, data_TRIM)

    for key in vacancy_keys:
        data_TRIM[key] = [0 if x is None else x for x in data_TRIM[key]]/np.asarray(data_TRIM["No_inc"])
       
    for key in data_keys:
            data_TRIM[key] = np.asarray(data_TRIM[key])
    return data_TRIM

xi = parse_arguments().xi
xi_string = str(xi).replace('.', '_')
output_dir = "SRIM_to_NIEL_21_eV_inf_" + xi_string + "_A"
# List of elements
Elements = [
    "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon",
    "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium",
    "Aluminium", "Silicon", "Phosphorus"
]
#Elements =["Carbon"]
# Dictionaries to hold the data for all elements
data_keys = [
    "Ion_ionization", "Ion_vacancy", "Ion_phonon", "Recoil_ionization", 
    "Recoil_vacancy", "Recoil_phonon", "Energy", "Total_vacancies", 
    "Single_vacancies", "Divacancies", "Trivacancies", "Tetra_vacancies", "Penta_vacancies", 
    "vacancies_6_10", "vacancies_11_50", "vacancies_51_100", "vacancies_101_200", 
    "vacancies_201_500", "vacancies_above_500", "Range", "Vacancy_per_ion", 
    "Energy_TRIM_sim", "No_inc", "Cluster_volumes_cuboid", "Cluster_lengths_cuboid", 
    "Cluster_radius", "Cluster_densities_sphere", "Cluster_densities_cuboid","Cluster_population_per_nm","Energy_TRIM","Energy_primary_recoils","NIEL",
    "NIEL_vac","NIEL_energy","NIEL_energy_vac"
]
data = {key: {} for key in data_keys}
# Load data for all elements and store in the nested dictionaries
for element in Elements:
    data_TRIM = load_data(element, data_keys[:37])
    for key in data_TRIM:
        if key not in ["NIEL","NIEL_vac","NIEL_energy", "NIEL_energy_vac"]:
            data[key][element] = data_TRIM[key]       
    data["NIEL"][element]=np.asarray(data["Ion_vacancy"][element])+np.asarray(data["Recoil_vacancy"][element])+np.asarray(data["Ion_phonon"][element])+np.asarray(data["Recoil_phonon"][element])
    data["NIEL_vac"][element]=np.asarray(data["Ion_vacancy"][element])+np.asarray(data["Recoil_vacancy"][element])
    data["NIEL_energy"][element]=np.multiply(np.asarray(data["NIEL"][element][0:]) / 100, np.asarray(data["Energy_TRIM_sim"][element]))
    data["NIEL_energy_vac"][element]=np.multiply(np.asarray(data["NIEL_vac"][element][0:]) / 100, np.asarray(data["Energy_TRIM_sim"][element]))
    np.savetxt(
        output_dir + "/NIEL_" + element + ".txt",
        np.column_stack(
            (
                np.asarray(data["Energy_TRIM_sim"][element]) * 1000, # in keV
                np.asarray(data["NIEL_energy"][element]) * 1000, # in keV
                np.asarray(data["NIEL_energy_vac"][element]) * 1000, # in keV
                np.asarray(data["Range"][element])/10000, # in um
                np.asarray(data["Vacancy_per_ion"][element])
            )
        ),
        delimiter="\t",
    )
    matching_indexes = np.where(np.isin(data["Energy_TRIM_sim"][element], data["Energy"][element]))[0]
    #principal debugging
    np.savetxt(
        output_dir + "/NIEL_OPTICS_" + element + ".txt",
        np.column_stack(
            (
                np.asarray(data["Energy"][element]) * 1000, # in keV
                np.asarray(data["NIEL_energy"][element][matching_indexes]) * 1000, # in keV
                np.asarray(data["NIEL_energy_vac"][element][matching_indexes]) * 1000, # in keV
                np.take(data["Range"][element],matching_indexes)/10000, # in um
                np.take(data["Vacancy_per_ion"][element],matching_indexes),               
                np.asarray(data["Single_vacancies"][element]),
                np.asarray(data["Divacancies"][element]),
                np.asarray(data["Trivacancies"][element]),
                np.asarray(data["Tetra_vacancies"][element]),
                np.asarray(data["Penta_vacancies"][element]),
                np.asarray(data["vacancies_6_10"][element]),
                np.asarray(data["vacancies_11_50"][element]),
                np.asarray(data["vacancies_51_100"][element]),
                np.asarray(data["vacancies_101_200"][element]),
                np.asarray(data["vacancies_201_500"][element]),
                np.asarray(data["vacancies_above_500"][element]),
                np.asarray(data["Total_vacancies"][element]),
                np.asarray(data["Cluster_volumes_cuboid"][element]),
                np.asarray(data["Cluster_lengths_cuboid"][element]),
                np.asarray(data["Cluster_radius"][element]),
                np.asarray(data["Cluster_densities_sphere"][element]),
                np.asarray(data["Cluster_densities_cuboid"][element]),
                np.asarray(data["Cluster_population_per_nm"][element]),
                np.asarray(data["No_inc"][element]), 
            )
        ),
        delimiter="\t",
    ) 