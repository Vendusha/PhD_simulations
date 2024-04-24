import numpy as np
import os
import sklearn.cluster
import copy
import miniball
OPTICS = sklearn.cluster.OPTICS
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"xtick.labelsize": 15})
plt.rcParams.update({"ytick.labelsize": 15})
plt.rcParams.update({"axes.labelsize": 15})
class RecoilsFrom1PKA:
    def __init__(self):
        self.recoil_proton_numbers = []
        self.recoil_positions = []
        self.recoil_energies = []
        self.recoil_positions_end = []
        self.ion_origin = 0
        self.ion_energy = 0
        self.ion_position = []
        self.last_ion_energy = 0
        self.last_ion_position = []
        self.vacancies_total = 0
        self.vacancies = []
        self.slices_2d = [[] for i in range(100)]
        self.no_clusters = "None"
        self.cluster_populations = []
        self.cluster_energies = []
        self.cluster_volumes_cuboid = []
        self.cluster_radius = []
        self.cluster_lengths_cuboid = []
        self.cluster_densities_sphere = []
        self.cluster_densities_cuboid = []
        self.defects_isolated = "None"
        self.defects_clustered = "None"
        self.energy_isolated = 0
        self.energy_clustered = 0
        self.energy_primary_recoils = []
        self.cluster_centers = []
        self.isolated_defect_distance = []
        self.isolated_defect_distance_x = []
        self.clustered_defect_distance = []
        self.clustered_defect_distance_x = []
        self.defects_total = "None"
        self.track_2Dsection = []
        self.track_3D = []
        self.cluster_colors = []
        self.cluster_colors_2D = []
        self.defects_total_2D = "None"
        self.max_x = "None"
        self.min_x = "None"
        self.max_y = "None"
        self.min_y = "None"
        self.max_z = "None"
        self.min_z = "None"

    def add_event(self, proton_numbers, positions, energies, energy_primary_recoil, positions_end):
        self.recoil_proton_numbers.append(proton_numbers)
        self.recoil_positions.append(positions)
        self.recoil_energies.append(energies)
        self.recoil_positions_end.append(positions_end)
        self.energy_primary_recoils.append(energy_primary_recoil)

    def add_origins(self, primary_event, ion_origin, ion_energy, ion_position):
        self.primary_event = primary_event
        self.ion_origin = ion_origin
        self.ion_energy = ion_energy
        self.ion_position = ion_position


def load_data(folder, particle):
    summary_files = [f for f in os.listdir(str(folder)) if f == "TDATA.txt"]
    summary_file = str(folder) + "/" + summary_files[0]
    txt_files = [f for f in os.listdir(str(folder)) if f == "COLLISON.txt"]
    if len(txt_files) != 1:
        raise ValueError("COLLISON file not found.")
    filename_txt = str(folder) + "/" + txt_files[0]

    # read the .dat file, skip the first 10 lines and extract primary_event, ion_origin, ion_energy, ion_position
    PKA = []
    if particle == "Si":
        with open(summary_file, "r", encoding="latin") as f:
            lines = f.readlines()
            for line in lines:
                if "Total Ions calculated" in line:
                    no_incident = int(line.split()[4])
                    print(no_incident)
        for i in range(0, no_incident):
            recoils = RecoilsFrom1PKA()
            #split the name of the folder by / and take the last element
            energy = float(folder.split("/")[-1])
            primary_event = i
            ion_origin = 14
            ion_energy = energy
            ion_position = [0, 0, 0]
            recoils.add_origins(primary_event, ion_origin, ion_energy, ion_position)
            PKA.append(recoils)
    else:
        print("This ana is not suitable for neutrons and so on. Just ion cascades.")

    skip = False
    with open(filename_txt, "r", encoding="latin") as f:
        # lines = f.readlines()
        recoil_index = 0
        # lines_iter = iter(lines)
        for line in f:
            if skip:
                skip = False
                continue
            if line[2:5] == "Ion":
                recoil_index += 1
            elif line[0]== "³" and ("Start of New Cascade" in line):
                last_ion_energy= line.split("³")[2]
                x_last = line.split("³")[3]
                y_last = line.split("³")[4]
                z_last = line.split("³")[5]
                position_last = [float(x_last) / 10000, float(y_last) / 10000, float(z_last) / 10000]
                PKA[recoil_index-1].last_ion_energy = last_ion_energy
                PKA[recoil_index-1].last_ion_position = position_last
            if line[0] == "Û":
                try:
                    _, _, proton_number, energy, x, y, z, vac, repl, _ = line.split()
                    energy_primary_recoil = 0.0
                except ValueError:
                    _, _, proton_number, energy, x, y, z, vac, repl, _, _ = line.split()
                    energy_primary_recoil = float(energy)
                    
                if repl == "01":
                    skip = True
                    continue
                else:
                    position = [
                        float(x) / 10000,
                        float(y) / 10000,
                        float(z) / 10000,
                    ]  # in um
                    try:
                        PKA[recoil_index - 1].add_event(
                            int(proton_number), position, float(energy), energy_primary_recoil, [0, 0, 0]
                        )
                    except IndexError:
                        pass
                    # try:
                    #     PKA[recoil_index-1].slices_2d[math.floor(position[0])].append(position)
                    # except IndexError:
                    #     pass
                    # PKA[recoil_index-1].vacancies_total += 1
    print("Loaded .txt file")
    return PKA, no_incident
