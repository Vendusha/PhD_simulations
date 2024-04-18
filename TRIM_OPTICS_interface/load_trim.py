import numpy as np
import os
import sklearn.cluster
import copy
#import SimpleHists as sh
import miniball
import trim_plotting_helper as tph
from collections import Counter
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

    def add_max_min(self):
        self.max_x = np.max(self.recoil_positions, axis=0)[0]
        self.min_x = np.min(self.recoil_positions, axis=0)[0]
        self.max_y = np.max(self.recoil_positions, axis=0)[1]
        self.min_y = np.min(self.recoil_positions, axis=0)[1]
        self.max_z = np.max(self.recoil_positions, axis=0)[2]
        self.min_z = np.min(self.recoil_positions, axis=0)[2]

    def flatten(self):
        for i in range(len(self.recoil_positions)):
            self.recoil_positions[i][1] -= self.ion_position[1]
            self.recoil_positions[i][2] -= self.ion_position[2]
        self.ion_position[1] = 0
        self.ion_position[2] = 0
        self.recoil_positions.append(self.ion_position)
        self.max_x = np.max(self.recoil_positions, axis=0)[0]
        self.min_x = np.min(self.recoil_positions, axis=0)[0]
        self.max_y = np.max(self.recoil_positions, axis=0)[1]
        self.min_y = np.min(self.recoil_positions, axis=0)[1]
        self.max_z = np.max(self.recoil_positions, axis=0)[2]
        self.min_z = np.min(self.recoil_positions, axis=0)[2]

    def cluster_analysis(
        self, track, min_s, par, method, showPlot=False, cutoff_1=0, cutoff_2=1
    ):     
        Isolated_defects_distance = []
        Isolated_defects_distance_x = []
        Clustered_defects_distance = []
        Clustered_defects_distance_x = []
        self.defects_total = np.size(track) / 3
        self.divacancy = []
        self.trivacancy = []
        self.tetravacancy = []
        self.pentavacancy = []
        if min_s > np.size(track) / 3:
            self.defects_isolated = np.size(track) / 3
            for index in range(0,(int(np.size(track)/3))):
                Isolated_defects_distance.append(np.sqrt((track[index][0])**2+track[index][1]**2)+track[index][2]**2)
                Isolated_defects_distance_x.append(np.sqrt(track[index][0])**2)
            self.isolated_defect_distance=np.ravel(Isolated_defects_distance)
            self.isolated_defect_distance_x=np.ravel(Isolated_defects_distance_x)
            self.defects_clustered = 0
            self.energy_isolated = np.sum(self.recoil_energies)
            self.energy_clustered = 0
            self.no_clusters = 0
            self.cluster_colors = np.array(["black"] * int(self.defects_total))
            self.cluster_markers = np.array(["."] * int(self.defects_total))
        else:
            Energy = self.recoil_energies
            # Energy.append(self.ion_energy)
            Energy = np.array(Energy)
            if showPlot:
                cluster_centers = []
                cluster_center_annotations = []
            if method == "xi":
                model = OPTICS(
                    min_samples=min_s, cluster_method="xi", metric="minkowski", xi=par
                )
            if method == "dbscan":
                par = par / 10000
                model = OPTICS(
                    min_samples=min_s,
                    cluster_method="dbscan",
                    eps=par,
                    metric="minkowski",
                )
            y_pred = model.fit_predict(track)
            indices_isolated = np.where(y_pred == -1)
            indices_clustered = np.where(y_pred != -1)
            # plot reachability distance            
            self.cluster_colors = np.array(
                ["black"] * int(len(y_pred)), dtype=object
            )
            self.cluster_colors[y_pred == -1] = "red"
            #print(self.cluster_colors)
            self.cluster_markers = np.array(
                ["."] * int(self.defects_total), dtype=object
            )
            self.energy_isolated = np.sum(Energy[indices_isolated])
            self.energy_clustered = np.sum(Energy[indices_clustered])
            self.defects_isolated = np.size(indices_isolated)
            for index in (indices_isolated):
                Isolated_defects_distance.append(np.sqrt((track[index,0])**2+track[index,1]**2)+track[index,2]**2)
                Isolated_defects_distance_x.append(np.sqrt(track[index,0])**2)
            self.isolated_defect_distance=np.ravel(Isolated_defects_distance)
            self.isolated_defect_distance_x=np.ravel(Isolated_defects_distance_x)
            self.defects_clustered = np.size(indices_clustered)
            
            #print(Energy)
            #print(self.recoil_energies)
            clusters = np.unique(y_pred)
            clusters = clusters[clusters != -1]
            self.no_clusters = np.size(clusters)
            for n in clusters:
                #print(clusters)
                current_cluster_idx = np.where(y_pred == n)[0]
                #print(y_pred)
                #print(current_cluster_idx)
                unique_array = copy.deepcopy(track[current_cluster_idx])
                if np.size(current_cluster_idx) == 2:
                    self.cluster_colors[np.where(y_pred == n)] = ["blue"]
                    #print("gotcha")
                    self.cluster_markers[np.where(y_pred == n)] = ["o"]
                    self.divacancy.append((unique_array[0].tolist()))
                elif np.size(current_cluster_idx) == 3:
                    self.cluster_colors[np.where(y_pred == n)] = ["darkgreen"]
                    self.trivacancy.append((unique_array[0].tolist()))
                elif np.size(current_cluster_idx) == 4:
                    self.cluster_colors[np.where(y_pred == n)] = ["black"]
                    self.tetravacancy.append((unique_array[0].tolist()))
                elif np.size(current_cluster_idx) == 5:
                    self.cluster_colors[np.where(y_pred == n)] = ["magenta"]
                    self.pentavacancy.append((unique_array[0].tolist()))
                elif np.size(current_cluster_idx) > 5 and np.size(current_cluster_idx) < 10:
                    self.cluster_colors[np.where(y_pred == n)] = ["cyan"]
                elif np.size(current_cluster_idx) >= 10 and np.size(current_cluster_idx) < 50:
                    self.cluster_colors[np.where(y_pred == n)] = ["orange"]
                elif np.size(current_cluster_idx) >= 50 and np.size(current_cluster_idx) < 100:
                    self.cluster_colors[np.where(y_pred == n)] = ["gray"]
                elif np.size(current_cluster_idx) >= 100:
                    self.cluster_colors[np.where(y_pred == n)] = ["purple"]
                x_max = np.max(track[current_cluster_idx][:, 0])
                x_min = np.min(track[current_cluster_idx][:, 0])
                y_max = np.max(track[current_cluster_idx][:, 1])
                y_min = np.min(track[current_cluster_idx][:, 1])
                z_max = np.max(track[current_cluster_idx][:, 2])
                z_min = np.min(track[current_cluster_idx][:, 2])
                diff_x = np.abs(x_max - x_min)
                diff_y = np.abs(y_max - y_min)
                diff_z = np.abs(z_max - z_min)
                if diff_x == 0:
                    diff_x = 0.0001
                    unique_array[0, 0] = unique_array[0, 0] + 0.005
                if diff_y == 0:
                    diff_y = 0.0001
                    unique_array[0, 1] = unique_array[0, 1] + 0.005
                if diff_z == 0:
                    diff_z = 0.0001
                    unique_array[0, 2] = unique_array[0, 2] + 0.005
                if len(unique_array) == 2:
                    cluster_center = np.mean(unique_array, axis=0)
                    radius = (
                        np.sqrt((diff_x) ** 2 + (diff_y) ** 2 + (diff_z) ** 2)
                    ) / 2
                else:
                    try:
                        cluster_center, radius_2 = miniball.get_bounding_ball(
                            (np.unique(unique_array, axis=0))
                        )
                        radius = np.sqrt(radius_2)
                    except:
                        cluster_center = np.mean(unique_array, axis=0)
                        radius = (
                            np.sqrt((diff_x) ** 2 + (diff_y) ** 2 + (diff_z) ** 2)
                        ) / 2

                volume = np.abs(diff_x) * np.abs(diff_y) * np.abs(diff_z)
                cluster_length = np.sqrt((diff_x) ** 2 + (diff_y) ** 2 + (diff_z) ** 2)
                self.cluster_densities_sphere.append(
                    np.size(current_cluster_idx) / (4*np.pi*radius**3)
                    )
                self.cluster_densities_cuboid.append(
                    np.size(current_cluster_idx) / volume
                )
                self.cluster_populations.append(np.size(current_cluster_idx))
                self.cluster_energies.append(np.sum(Energy[current_cluster_idx]))
                self.cluster_volumes_cuboid.append(volume)
                self.cluster_centers.append(cluster_center)
                Clustered_defects_distance.append(np.sqrt((cluster_center[0])**2+cluster_center[1]**2)+cluster_center[2]**2)
                Clustered_defects_distance_x.append(cluster_center[0])
                self.cluster_lengths_cuboid.append(cluster_length)
                self.cluster_radius.append(radius)
                if showPlot:
                    current_cluster_indeces = np.where(y_pred == n)
                    cluster_center_annotation = int(np.median(current_cluster_indeces))
                    cluster_center_annotations.append(cluster_center_annotation)
            self.clustered_defect_distance= np.ravel(Clustered_defects_distance)
            self.clustered_defect_distance_x= np.ravel(Clustered_defects_distance_x)
            if showPlot:
                fig_3D = plt.figure()
                ax_3D = fig_3D.add_subplot(111, projection="3d")
                #ax_3D.set_title("3d projection plot")
                ax_3D.set_xlabel("X [nm]", labelpad=10)
                ax_3D.set_ylabel("Y [nm]", labelpad=10)
                ax_3D.set_zlabel("Z [nm]", labelpad=10)
                ax_3D.scatter(
                track[:, 0]*1000, track[:,1]*1000, track[:,2]*1000, c=self.cluster_colors, s=20
                ) 
                plt.show()
               
                # print(self.divacancy)
                # print(self.trivacancy)
                # a=track[:,0]
                # arg_2Dsection=np.ravel(np.argwhere((a>cutoff_1)&(a<cutoff_2)))
                # self.track_2Dsection = track[arg_2Dsection]
                # self.defects_total_2D=np.size(track[arg_2Dsection])/3
                # self.track_3D = track
                # norm = mpl.colors.Normalize()
                # cmap = cm.viridis
                # m = cm.ScalarMappable(norm=norm, cmap=cmap)
                # self.cluster_colors = m.to_rgba(y_pred)
                # indeces_isolated = np.where(y_pred == -1)[0]
                # write red color in rgba format
                # red = [1.0, 0.0, 0.0, 1]
        #tph.plot_reachability(np.arange(len(track)), model, self.cluster_colors)

def load_data(folder, particle):
    if particle != "Si" and particle != "Ian":
        dat_files = [f for f in os.listdir(str(folder)) if f.endswith(".dat")]
        if len(dat_files) != 1:
            raise ValueError("should be only one .dat file in the current directory")
        filename_dat = str(folder) + "/" + dat_files[0]
        shist_files = [f for f in os.listdir(str(folder)) if f == "PKA_lin.shist"]
        if len(shist_files) != 1:
            raise ValueError("should be only one .shist file in the current directory")
        filename_shist = str(folder) + "/" + shist_files[0]
    if particle !="Ian":
        summary_files = [f for f in os.listdir(str(folder)) if f == "TDATA.txt"]
        summary_file = str(folder) + "/" + summary_files[0]
    txt_files = [f for f in os.listdir(str(folder)) if f == "COLLISON.txt"]
    if len(txt_files) != 1:
        raise ValueError("COLLISON file not found.")
    filename_txt = str(folder) + "/" + txt_files[0]

    # read the .dat file, skip the first 10 lines and extract primary_event, ion_origin, ion_energy, ion_position
    PKA = []
    if particle != "Si" and particle != "Ian":
        with open(filename_dat, "r") as f:
            lines = f.readlines()
            ion_index = 0
            for line in lines[10:]:
                recoils = RecoilsFrom1PKA()
                primary_event = int(line.split()[0])
                ion_origin = int(line.split()[1])
                ion_energy = float(line.split()[2])
                ion_position = [
                    float(line.split()[3]) / 10000,
                    float(line.split()[4]) / 10000,
                    float(line.split()[5]) / 10000,
                ]
                recoils.add_origins(primary_event, ion_origin, ion_energy, ion_position)
                # try:
                #     PKA[ion_index].slices_2d[math.floor(position[0])].append(ion_position)
                # except IndexError:
                #     pass
                PKA.append(recoils)
                ion_index += 1
        print("Loaded .dat file")
        hist = sh.HistCollection(filename_shist)
        try:
            no_incident = (hist.hist("Primary_neutron")).getIntegral()
        except RuntimeError:
            try:
                no_incident = (hist.hist("Primary_proton")).getIntegral()
            except RuntimeError:
                try:
                    no_incident = (hist.hist("Primary_e-")).getIntegral()
                except RuntimeError:
                    print(
                        "Shist file does not contain any of the following: Primary_neutron, Primary_proton, Primary_e-"
                    )
    elif particle == "Si":
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
    elif particle == "Ian":
        no_incident = 10000
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
