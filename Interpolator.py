import matplotlib.pyplot as plt
import matplotlib.cm as cm
import miniball
import matplotlib as mpl
import math
import numpy as np
import sklearn.cluster
from labellines import labelLine, labelLines
import copy
import sklearn.cluster as cluster
OPTICS = sklearn.cluster.OPTICS
mpl.rc("font", size=16)


class Event:
    def __init__(self, name, energy):
        self.name = name
        self.energy = energy
        self.no_clusters = 0
        self.cluster_populations = []
        self.cluster_sizes = []
        self.cluster_sizes_ball = []
        self.defects_isolated = 0
        self.defects_clustered = 0
        self.energy_isolated = 0
        self.energy_clustered = 0
        self.defects_total = 0
        self.track_2Dsection = []
        self.track_3D = []
        self.cluster_colors = []
        self.cluster_colors_2D = []
        self.defects_total_2D=0

    def clustering(
            self, track, track_vacancies, Energy, min_s, threshold, par, method, showPlot=False
    ):
        # track=track_vacancies
        cutoff_1=0
        cutoff_2=50000
        # cutoff_1=400000
        # cutoff_2=500000
        self.defects_total = np.size(track) / 3
        if min_s >= np.size(track) / 3:
            self.defects_isolated = np.size(track) / 3
            self.energy_isolated = np.sum(Energy)
            a=track[:,0]
            arg_2Dsection=np.ravel(np.argwhere((a>cutoff_1)&(a<cutoff_2)))
            self.track_2Dsection = track[arg_2Dsection]
            self.defects_total=np.size(track[arg_2Dsection])/3
            self.track_3D = track
        else:
            # if showPlot:
            #     # cluster_centers = []
            #     cluster_center_annotations = []
            # if method =="xi":
            #     model = OPTICS(min_samples=min_s, cluster_method="xi", metric="minkowski", xi=par)
            # if method=="dbscan":
            #     model = OPTICS(
            #         min_samples=min_s, cluster_method="dbscan", eps=par, metric="minkowski"
            #     )
            # y_pred = model.fit_predict(track)
            # ######everything taken straight from optics ##################
            # indices_isolated = np.where(y_pred == -1)
            # indices_clustered = np.where(y_pred != -1)
            # self.energy_isolated = np.sum(Energy[indices_isolated])
            # self.energy_clustered = np.sum(Energy[indices_clustered])
            # self.defects_isolated = np.size(indices_isolated)
            # self.defects_clustered = np.size(indices_clustered)
            # clusters = np.unique(y_pred)
            # clusters = clusters[clusters != -1]
            # self.no_clusters = np.size(clusters)
            # for n in clusters:
            #     current_cluster_idx = np.where(y_pred == n)[0]
            #     # print(track[current_cluster_idx])
            #     unique_array = copy.deepcopy(track[current_cluster_idx])
            #     try:
            #         cluster_center, r_2 = miniball.get_bounding_ball(
            #         ((np.unique(unique_array, axis=0))))
            #     except:
            #         r_2=0
            #     # r_2 = 1
            #     x_max = (np.max(track[current_cluster_idx][:,0]))
            #     x_min = (np.min(track[current_cluster_idx][:,0]))
            #     y_max = (np.max(track[current_cluster_idx][:,1]))
            #     y_min = (np.min(track[current_cluster_idx][:,1]))
            #     z_max = (np.max(track[current_cluster_idx][:,2]))
            #     z_min = (np.min(track[current_cluster_idx][:,2]))
            #     volume = (x_max-x_min)*(y_max-y_min)*(z_max-z_min)
            #     r = volume
            #     # print(volume)
            #     r_ball = np.sqrt(r_2)
            #     self.cluster_populations.append(np.size(current_cluster_idx))
            #     self.cluster_sizes.append(r)
            #     self.cluster_sizes_ball.append(r_ball)
            #     if showPlot:
            #         current_cluster_indeces = np.where(y_pred == n)
            #         cluster_center_annotation = int(np.median(current_cluster_indeces))
            #         cluster_center_annotations.append(cluster_center_annotation)
            a=track[:,0]
            # arg_2Dsection=np.ravel(np.argwhere(a<cutoff_2))
            arg_2Dsection=np.ravel(np.argwhere((a>cutoff_1)&(a<cutoff_2)))
            self.track_2Dsection = track[arg_2Dsection]
            self.defects_total_2D=np.size(track[arg_2Dsection])/3
            self.track_3D = track
            # norm = mpl.colors.Normalize()
            # cmap = cm.viridis
            # m = cm.ScalarMappable(norm=norm, cmap=cmap)
            # self.cluster_colors = m.to_rgba(y_pred)
            # indeces_isolated = np.where(y_pred == -1)[0]
            # self.cluster_colors[:]=[0.0, 0.494, 1.0, 0.0]
            # self.cluster_colors[indeces_isolated] = [0.501, 0.501, 0.501, 1]
            # self.cluster_colors_2D=self.cluster_colors[arg_2Dsection]

            if showPlot:
                print(self.energy)
                norm = mpl.colors.Normalize()
                cmap = cm.viridis
                m = cm.ScalarMappable(norm=norm, cmap=cmap)
                # colors = m.to_rgba(y_pred_merged)
                # indeces_isolated = np.where(y_pred_merged==-1)[0]
                colors = m.to_rgba(y_pred)
                indeces_isolated = np.where(y_pred == -1)[0]
                colors[indeces_isolated] = [0.501, 0.501, 0.501, 1]
                fig = plt.figure(figsize=(20, 12))
                ax_reachability = fig.add_subplot(221)
                ax_event = fig.add_subplot(222, projection="3d")
                ax_event_xz = fig.add_subplot(223)
                ax_event_yz = fig.add_subplot(224)
                # ax_event = fig.add_subplot(122)
                space = np.arange(len(track))

                for no_cluster in range(0, self.no_clusters):
                    x_pos = cluster_center_annotations[no_cluster] / int(
                        len(model.reachability_)
                    )
                    if np.mod(no_cluster, 2) == 0:
                        y_pos = 0.4
                    else:
                        y_pos = 0.3

                    # ax_reachability.annotate(
                    #     "C "
                    #     + str(no_cluster + 1)
                    #     + "\n"
                    #     + "n="
                    #     + str(self.cluster_populations[no_cluster])
                    #     + "\n"
                    #     + str(np.round(self.cluster_sizes[no_cluster], 2))
                    #     + " nm",
                    #     xy=(x_pos, y_pos),
                    #     xycoords="axes fraction",
                    #     ha="center",
                    #     va="bottom",
                    #     fontsize=16,
                    # )

                # ax_event.scatter(track[:, 0], track[:, 2], c=colors, s=12)
                ######Figure 1
                idx_low=0
                idx_high=1
                ax_event.scatter(track[:, 0], track[:, 1], track[:, 2], c=colors, s=12)
                ax_event_xz.scatter(track[:, 0], track[:, 2], c=colors, s=12)
                ax_event_yz.scatter(track[:, 0], track[:, 1], c=colors, s=12)

                ax_event.scatter(track[:, 0][idx_low:idx_high], track[:, 1][idx_low:idx_high], track[:, 2][idx_low:idx_high], c=colors[idx_low:idx_high], s=12)
                ###to comment out####
                # ax_event.xlim(-70,0)
                # ax_event.ylim(0,200)
                # ax_event.zlim(0,50)
                
                # ax_event_xz.scatter(track[:, 0], track[:, 2], c=colors, s=12)
                # ax_event_yz.scatter(track[:, 0], track[:, 1], c=colors, s=12)

                ax_reachability.bar(
                    space, model.reachability_, align="center", color=colors
                )
                ax_reachability.set_ylabel("Distance [A]")
                ax_reachability.set_xlabel("Ordered silicon recoil [-]")
                # ax_reachability.axhline(threshold, label=str(threshold) + " nm")
                ax_reachability.set_title(
                    str(self.energy*1000)
                    + " keV incident recoil, n (neighbours)= "
                    + str(min_s)
                )
                labelLines(ax_reachability.get_lines(), ha="left")
                ax_reachability.annotate(
                    str(self.defects_isolated)
                    + " single defects\n"
                    + str(self.defects_clustered)
                    + " defects in clusters\n"
                    + str(self.no_clusters)
                    + " clusters formed",
                    xy=(0.4, 1),
                    xytext=(12, -12),
                    va="top",
                    xycoords="axes fraction",
                    textcoords="offset points",
                    fontsize=16,
                )


                # ax_event.tick_params(apars="both", which="both", labelsize=16)
                ax_event.set_xlabel("X [A]", labelpad=15)
                ax_event.set_ylabel("Y [A]", labelpad=15)
                ax_event.set_zlabel("Z [A]", labelpad=15)
                ax_event_xz.set_xlabel("X [A]", labelpad=15)
                ax_event_xz.set_ylabel("Z [A]", labelpad=15)
                ax_event_yz.set_xlabel("Y [A]", labelpad=15)
                ax_event_yz.set_ylabel("Z [A]", labelpad=15)
            plt.show()


def round_up_to_even(f):
    return math.ceil(f / 2.0) * 2
