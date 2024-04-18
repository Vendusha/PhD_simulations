import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv
import argparse
import _pickle as cPickle
import optics
import os
plt.rcParams.update({"font.size": 22})
plt.rcParams.update({"xtick.labelsize": 22})
plt.rcParams.update({"ytick.labelsize": 22})
showPlot=False
def analyze_event(filename,e, Track_ion,neighbours=9, r=4,par=0.05, method="dbscan"):
    x_ion, y_ion, z_ion, energy_ion = Track_ion
    X = []
    Y = []
    Z = []
    X_end = []
    Y_end = []
    Z_end = []
    Energy = []
    track = []
    track_vacancies = []
    Energy.append(energy_ion)
    track.append([x_ion, y_ion, z_ion])
    with open(filename) as csvfile:
        for row in csv.reader(csvfile, delimiter=" "):
            Z_number, x, y, z, energy, x_end, y_end, z_end = (
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
                float(row[4]),
                float(row[5]),
                float(row[6]),
                float(row[7]),
            )
            X.append(x)
            Y.append(y)
            Z.append(z)
            X_end.append(x_end)
            Y_end.append(y_end)
            Z_end.append(z_end)
            Energy.append(energy)
            track.append([x, y, z])
            track_vacancies.append([x_end, y_end, z_end])
    track = np.asarray(track)
    track_vacancies = np.asarray(track_vacancies)
    Energy = np.asarray(Energy)
    event = optics.Event(filename, e)
    event.clustering(track, track_vacancies, Energy, neighbours, r, par=par,method=method, showPlot=showPlot)
    return(event)
n = 20
palette = sns.color_palette(None, n)

parser = argparse.ArgumentParser(description="Visualize events.")
parser.add_argument(
    "folder", metavar="str", type=str, nargs="+", help="folder to load."
)
filen = parser.parse_args()
folder = filen.folder[0]
neighbours = int(filen.folder[1])
par = float(filen.folder[2])
method = str(filen.folder[3])
events = 0
e = float(folder)
# neighbours = 9
r = 30
Cluster_populations = np.array([])
Cluster_sizes = np.array([])
Defects_isolated = np.array([])
Defects_clustered = np.array([])
Defects_total = np.array([])
Defects_total_2D = np.array([])
Energy_isolated = np.array([])
Energy_clustered = np.array([])
No_clusters = np.array([])
Energy_init = np.array([])
Parent_event = np.array([])
Ion_origin = np.array([])
# filename_dat="proton_"+str(e)+"_MeV.dat"

txt_files = [f for f in os.listdir(str(folder)) if f.endswith('.dat')]
if len(txt_files) != 1:
    raise ValueError('should be only one txt file in the current directory')
filename_dat=str(folder)+"/"+txt_files[0]
index = 0
track_2Dsection = []
track_3D = []
cluster_colors = []
cluster_colors_2D = []
Ion_no = []
Track_ion = []
Ion_energy = []
with open(filename_dat) as file:
    rows = (line.split() for line in file)
    for row in rows:
        index+=1
        if index>10:
            parent_event, ion_origin, ion_energy, ionx, iony, ionz, _ioncosx, _ioncosy, _ioncosz = (
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
                float(row[4]),
                float(row[5]),
                float(row[6]),
                float(row[7]),
                float(row[8]),
            )
            Ion_origin=np.append(Ion_origin,ion_origin)
            Ion_energy=np.append(Ion_energy,ion_energy)
            Parent_event=np.append(Parent_event,parent_event)
            Track_ion.append([ionx, iony, ionz, ion_energy])
    Track_ion = np.asarray(Track_ion)


for file in os.listdir(folder):
    if file.endswith(".txt") and file[0]=="S":
        txt_split=(file.split("_")[7])
        ion_no=int(txt_split.split(".")[0])-2
        print(file)
        events += 1
        event = analyze_event((os.path.join(folder, file)), e, Track_ion[ion_no], neighbours=neighbours, r=r, par=par, method=method)
        Cluster_populations = np.append(Cluster_populations, event.cluster_populations)
        Cluster_sizes = np.append(Cluster_sizes, event.cluster_sizes)
        No_clusters = np.append(No_clusters, event.no_clusters)
        Energy_isolated = np.append(Energy_isolated, event.energy_isolated)
        Energy_clustered = np.append(Energy_clustered, event.energy_clustered)
        Defects_isolated = np.append(Defects_isolated, event.defects_isolated)
        Defects_clustered = np.append(Defects_clustered, event.defects_clustered)
        Defects_total = np.append(Defects_total, event.defects_total)
        Defects_total_2D = np.append(Defects_total_2D, event.defects_total_2D)
        Energy_init = np.append(Energy_init, event.energy)
        track_2Dsection.append(event.track_2Dsection)
        track_3D.append(event.track_3D)
        cluster_colors.append(event.cluster_colors)
        cluster_colors_2D.append(event.cluster_colors_2D)
        Ion_no.append(ion_no)
        # track_2Dsection=np.append(track_2Dsection, event.track_2Dsection)
        # cluster_colors=np.append(track_2Dsection, event.track_2Dsection)
Events = [events]
#Complete 3D plot
fig_3D = plt.figure()
ax_3D=fig_3D.add_subplot(111, projection="3d")
ax_3D.set_title("3d projection plot")
ax_3D.set_xlabel("X [um]", labelpad=15)
ax_3D.set_ylabel("Y [um]", labelpad=15)
ax_3D.set_zlabel("Z [um]", labelpad=15)
### huhtinen-like 2D plots, origins 2D plots, z is the incident beam
fig_yz_huhtinen = plt.figure()
ax_yz_huhtinen=fig_yz_huhtinen.add_subplot(111)
ax_yz_huhtinen.set_title("Huhtinen projection")
ax_yz_huhtinen.set_xlabel("Y [um]")
ax_yz_huhtinen.set_ylabel("Z [um]")
fig_xy_origins = plt.figure()
ax_xy_origins=fig_xy_origins.add_subplot(111)
ax_xy_origins.set_title("Alpha/Si-recoils origin")
ax_xy_origins.set_xlabel("X [um]")
ax_xy_origins.set_ylabel("Y [um]")
fig_xy_origins_cut = plt.figure()
# ax_xy_origins_cut=fig_xy_origins.add_subplot(111)
# ax_xy_origins_cut.set_title("Alpha/Si-recoils origin")
# ax_xy_origins_cut.set_title("Alpha/Si-recoils origin")
# ax_xy_origins_cut.set_xlabel("X [A]")
# ax_xy_origins_cut.set_ylabel("Y [A]")
fig_yz_origins = plt.figure()
ax_yz_origins=fig_yz_origins.add_subplot(111)
ax_yz_origins.set_title("Alpha/Si-recoils origin")
ax_yz_origins.set_xlabel("Y [um]")
ax_yz_origins.set_ylabel("Z [um]")

# fig1 = plt.figure()
# fig1_a = plt.figure()
#huhtinen projection
# fig3 = plt.figure()
# ax1=fig1.add_subplot(111)
# ax1_a=fig1_a.add_subplot(111)
# ax1_a.set_title("2dxz projection plot")
# ax1.set_title("2dxy projection plot")
# ax3=fig3.add_subplot(111)

for idx, track_i in enumerate((track_3D)):
    # print(track_2Dsection[idx])
    Ion=(Track_ion[Ion_no[idx]])
    # ax_3D.scatter(Ion[0], Ion[1], Ion[2], c="m", s=10)
    # print(str(Ion_no[idx])+"  event is ion number "+str(Ion_origin[Ion_no[idx]]))
    ax_yz_huhtinen.scatter(track_2Dsection[idx][:, 1]/10000, track_2Dsection[idx][:, 2]/10000, c="black", s=0.2)
    # ax_yz_huhtinen.scatter(track_2Dsection[idx][:, 1]/10000, track_2Dsection[idx][:, 2]/10000, c=cluster_colors_2D[idx], s=0.05)
    if Ion_origin[Ion_no[idx]]==2:
        ax_3D.scatter(track_3D[idx][:, 0]/10000, track_3D[idx][:, 1]/10000, track_3D[idx][:, 2]/10000, c="blue", s=0.2)
        ax_xy_origins.scatter(track_3D[idx][:, 0]/10000, track_3D[idx][:, 1]/10000, c="blue", s=0.05)
        ax_xy_origins.scatter(track_2Dsection[idx][:, 0]/10000, track_2Dsection[idx][:, 1]/10000, c="purple", s=0.2)
        ax_yz_origins.scatter(track_2Dsection[idx][:, 1]/10000, track_2Dsection[idx][:, 2]/10000, c="blue", s=0.2)
    else:
        ax_3D.scatter(track_3D[idx][:, 0]/10000, track_3D[idx][:, 1]/10000, track_3D[idx][:, 2]/10000, c="red", s=0.2)
        ax_xy_origins.scatter(track_3D[idx][:, 0]/10000, track_3D[idx][:, 1]/10000, c="red", s=0.05)
        ax_xy_origins.scatter(track_2Dsection[idx][:, 0]/10000, track_2Dsection[idx][:, 1]/10000, c="orange", s=0.2)
        ax_yz_origins.scatter(track_2Dsection[idx][:, 1]/10000, track_2Dsection[idx][:, 2]/10000, c="red", s=0.2)
            # ax_xy_huhtinen.scatter(track_i[:, 0], track_i[:, 1], c="blue", s=0.1)
        # color=cluster_colors[idx]
        # # ax3.scatter(track_i[:, 2], track_i[:, 1], c=color, s=5)
        # if Ion_origin[Ion_no[idx]]==2:
        #     ax1.scatter(track_i[:, 0], track_i[:, 1], c="blue", s=1)
        #     ax1_a.scatter(track_i[:, 0], track_i[:, 2], c="blue", s=1)
        #     ax1_b.scatter(track_i[:, 1], track_i[:, 2], c="blue", s=1)
        #     ax.scatter(track_3D[idx][:, 0], track_3D[idx][:, 1], track_3D[idx][:, 2], c="blue", s=1)
        # else:
        #     ax1.scatter(track_i[:, 0], track_i[:, 1], c="red", s=1)
        #     ax1_a.scatter(track_i[:, 0], track_i[:, 2], c="red", s=1)
        #     ax1_b.scatter(track_i[:, 1], track_i[:, 2], c="red", s=1)
        #     ax.scatter(track_3D[idx][:, 0], track_3D[idx][:, 1], track_3D[idx][:, 2], c="red", s=1)
# for idx, track_i in enumerate((track_2Dsection)):
#         color=cluster_colors[idx]ax_yz_origins.set_xlim(0,10000)
total_defects=np.sum(Defects_total)
total_defects_2D=np.sum(Defects_total_2D)
ax_yz_origins.set_ylim(0,1)
ax_yz_origins.set_xlim(0,1)
ax_yz_origins.set_title(str(total_defects_2D)+" vacancies")
# ax_xy_origins.set_ylim(-5,5)
ax_yz_huhtinen.set_ylim(0,1)
ax_yz_huhtinen.set_xlim(0,1)
ax_yz_huhtinen.set_title(str(total_defects_2D)+" vacancies")
# ax_3D.set_ylim(-2.5,2.5)
# ax_3D.set_zlim(-2.5,2.5)
ax_yz_origins.figure.savefig("Yz_origin.png",dpi=600)
# ax_xy_origins_cut.figure.savefig("Xy_origin_cut.png",dpi=600)
ax_xy_origins.figure.savefig("Xy_origin.png",dpi=600)
ax_yz_huhtinen.figure.savefig("Yz_huhtinen.png",dpi=600)
ax_3D.figure.savefig("3D.png",dpi=600)
# ax_xy_origins.figure.savefig("Xy_origin_zoom_4_5.png",dpi=600)
# ax_xy_origins.set_xlim(3000,4000)
# ax_xy_origins.set_ylim(3000,4000)
# ax_xy_origins.figure.savefig("Xy_origin_zoom_3_4.png",dpi=600)
# ax_xy_origins.set_xlim(2000,3000)
# ax_xy_origins.set_ylim(3000,4000)
# ax_xy_origins.figure.savefig("Xy_origin_zoom_3_4.png",dpi=600)

# ax_xy_origins.set_xlim(0,1000)
# ax_xy_origins.set_ylim(0,1000)
# ax_xy_origins.figure.savefig("Xy_origin_zoom_0_1.png",dpi=600)
# plt.show()
all = [Cluster_populations, Cluster_sizes,
      No_clusters, Energy_isolated, Energy_clustered,Defects_clustered, Defects_total, Defects_isolated, Events]
label = str(folder) + "_" + str(neighbours) + "_neighbours_" + str(r) + "_eps"
cPickle.dump(all, open(label, "wb"))
