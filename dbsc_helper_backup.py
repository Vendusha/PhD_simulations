import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import numpy as np
import sklearn.cluster
from matplotlib import gridspec

DBSCAN = sklearn.cluster.DBSCAN
AGGLO = sklearn.cluster.AgglomerativeClustering
OPTICS = sklearn.cluster.OPTICS
import matplotlib.cm as cm
from labellines import labelLine, labelLines
import miniball


def round_up_to_even(f):
    return math.ceil(f / 2.0) * 2


def plt_sphere(ax, list_center, list_radius):
    for c, r in zip(list_center, list_radius):
        u, v = np.mgrid[0 : 2 * np.pi : 50j, 0 : np.pi : 50j]
        x = r * np.cos(u) * np.sin(v)
        y = r * np.sin(u) * np.sin(v)
        z = r * np.cos(v)
        ax.plot_surface(
            x - c[0],
            y - c[1],
            z - c[2],
            color="orange",
            alpha=0.5 * np.random.random() + 0.5,
        )
        return ax

def dbsc_vendula(no_subfig, track, eps, min_s, cols):
    """Density
    based scan, explained here: https://towardsdatascience.com/dbscan-clustering-explained-97556a2ad556"""
    fig1 = plt.figure(figsize=(12, 12))
    rows = math.ceil(no_subfig / cols)
    grid = gridspec.GridSpec(rows, cols)
    idx_subfig = 0
    for epsilon in eps:
        for min_samples in min_s:
            model = DBSCAN(eps=epsilon, min_samples=min_samples)
            model.fit(track)
            y_pred = model.fit_predict(track)
            labels = model.labels_
            ax4 = fig1.add_subplot(grid[idx_subfig], projection="3d")
            colormap = plt.cm.hot
            normalize = plt.Normalize(vmin=min(labels[labels != 0]), vmax=max(labels))
            print(labels)
            # ax4.scatter(track[:, 0], track[:, 1],
            # track[:, 2],c=labels, s=50, alpha=1.0)
            # for idx, label in enumerate(labels):
            #     if label ==-1:
            #         ax4.scatter(track[idx, 0], track[idx, 1],
            #                     track[idx, 2],color="gray", s=50)
            # else:
            #     ax4.scatter(track[idx, 0], track[idx, 1],
            #                 track[idx, 2],c=label/max(labels), cmap=colormap, s=50, alpha=1.0)

            #         labels[idx] = 1
            # ax4.scatter(track[:, 0], track[:, 1],
            # track[:, 2], c=model.labels_, s=4)
            # ax4.scatter(track[:, 0], track[:, 1],
            # track[:, 2], c=y_pred, s=4)
            # idx_l=60
            # idx_h=62
            # ax4.scatter(track[:idx_l, 0], track[:idx_l, 1],
            # track[:idx_l, 2], c="gray", s=50)
            # ax4.scatter(track[idx_h+1:, 0], track[idx_h+1:, 1],
            # track[idx_h+1:, 2], c="gray", s=50)
            ax4.scatter(track[:, 0], track[:, 1], track[:, 2], c="gray", s=50)
            # ax4.scatter(track[idx_h, 0], track[idx_h, 1],
            # track[idx_h, 2], c="pink", s=2500, alpha=0.7)
            # ax4.scatter(track[idx_h, 0], track[idx_h, 1],
            # track[idx_h, 2],color=(1,0,0,1), s=100, marker="x")
            # for idxn in range(idx_h-idx_l+1):
            # idx=idx_l+idxn
            idx = model.core_sample_indices_
            ix = 3
            # idx=np.where(labels!=-1)
            # idx=np.delete(idx[0],ix,0)
            idx = np.delete(idx, ix, 0)
            print(idx)
            ax4.scatter(
                track[ix + 2, 0],
                track[ix + 2, 1],
                track[ix + 2, 2],
                c=labels[ix + 2],
                s=50,
                alpha=1.0,
            )
            ax4.scatter(
                track[idx, 0], track[idx, 1], track[idx, 2], c="red", s=50, alpha=1.0
            )

            # ax4.scatter(track[idx_l, 0], track[idx_l, 1],
            #                     track[idx_l, 2],color=(1,0,0,1), s=50)

            idx_subfig += 1
            ax4.tick_params(axis="both", which="both", labelsize=16)
            ax4.set_xlabel("X [nm]", fontsize=16, labelpad=10)
            ax4.set_ylabel("Y [nm]", fontsize=16, labelpad=10)
            ax4.set_zlabel("Z [nm]", fontsize=16, labelpad=10)
            # id_track=20
            # plt_sphere(ax4, [(track[id_track,0],track[id_track,1],track[id_track,2])], [5])
            # ax.auto_scale_xyz()
            ax4.annotate(
                "min neighbours: "
                + str(min_samples)
                + "\n"
                + "r="
                + str(epsilon)
                + " nm",
                xy=(0.4, 1),
                xytext=(12, -12),
                va="top",
                xycoords="axes fraction",
                textcoords="offset points",
                fontsize=16,
            )
            # print(db.labels_[db.labels_ == -1].size+" point defects.")
    fig1.tight_layout()
    plt.savefig("Tutorial_9")
    

def optics_vendula(
    no_subfig,
    track,
    min_s,
    cols,
    threshold_array,
    track_vacancies,
    automatic_color=True,
    showPlot=False,
):
    """Density
    based scan"""
    # cols=round_up_to_even(no_subfig/5*3)
    # for threshold in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    energy_figure = "50keV"
    if automatic_color == True:
        threshold_array = [4]
    for threshold in threshold_array:
        cols = 2
        if showPlot:
            fig1 = plt.figure(figsize=(20, 12))
        rows = math.ceil(no_subfig / cols)
        grid = gridspec.GridSpec(rows, cols)
        idx_subfig = 0
        space = np.arange(len(track))
        cluster = np.empty(len(space))
        # to do: plot also vacancies maybe?
        for min_samples in min_s:
          if min_samples >= np.size(track)/3:
              no_clusters = 0
              no_clusters_optics = 0
              cluster_centers = []
              cluster_centers_optics = []
              cluster_center_annotations = []
              cluster_populations = []
              cluster_populations_optics = []
              cluster_sizes = []
              cluster_sizes_optics = []
              single_defects = np.size(track)
              defects_in_clusters = 0
          else:

            if showPlot:
                ax3 = fig1.add_subplot(grid[idx_subfig])
            single_defects = 0
            defects_in_clusters = 0
            model = OPTICS(
                min_samples=min_samples, cluster_method="xi", metric="minkowski"
            )
            # model2 = OPTICS(
                # min_samples=min_samples, cluster_method="xi", metric="minkowski"
            # )


            y_pred = model.fit_predict(track)
            y_cluster = y_pred
            y_merged_cluster = y_pred
            norm = mpl.colors.Normalize()
            cmap = cm.viridis
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            colors_y_pred = m.to_rgba(y_pred)
            new_cluster = False
            new_cluster_optics = False
            no_clusters = 0
            no_clusters_merged = 0
            no_clusters_optics = 0
            cluster_centers = []
            cluster_centers_optics = []
            cluster_center_annotations = []
            cluster_populations = []
            cluster_populations_optics = []
            cluster_populations_merged = []
            cluster_sizes = []
            cluster_sizes_optics = []
            cluster_sizes_merged = []
            # cluster_merging = np.array([])
            # indeces = np.where(y_cluster==-1)
            # colors_y_pred[indeces] = [0.501, 0.501, 0.501, 1]
            for idx, colors in enumerate(model.reachability_):
                if model.reachability_[idx] > threshold:
                    if new_cluster == True:
                        new_cluster = False
                        cluster_center_annotation = (
                            cluster_start_idx + (idx - cluster_start_idx) / 2
                        )
                        cluster_center, r_2 = miniball.get_bounding_ball(
                            track[cluster_start_idx:idx]
                        )
                        cluster_size = np.sqrt(r_2)
                        cluster_population = idx - cluster_start_idx
                        cluster_centers.append(cluster_center)
                        cluster_center_annotations.append(cluster_center_annotation)
                        cluster_populations.append(cluster_population)
                        cluster_sizes.append(np.round(cluster_size, 1))
                    colors_y_pred[idx] = [0.501, 0.501, 0.501, 1]
                    cluster[idx] = 0
                    single_defects += 1
                    y_cluster[idx] = -1
                else:
                    if new_cluster == False:
                        new_cluster = True
                        no_clusters += 1
                        cluster_start_idx = idx
                    if automatic_color == False:
                        colors_y_pred[idx] = [0.545, 0, 0, 1]
                    defects_in_clusters += 1
                    cluster[idx] = 1
            if new_cluster == True:
                new_cluster = False
                cluster_center_annotation = (
                    cluster_start_idx + (idx - cluster_start_idx) / 2
                )
                cluster_center, r_2 = miniball.get_bounding_ball(
                    track[cluster_start_idx:idx]
                )
                cluster_size = np.sqrt(r_2)
                cluster_population = idx - cluster_start_idx
                cluster_centers.append(cluster_center)
                cluster_center_annotations.append(cluster_center_annotation)
                cluster_populations.append(cluster_population)
                cluster_sizes.append(np.round(cluster_size, 2))

            clusters_optics = np.unique(y_cluster)
            clusters_optics = clusters_optics[clusters_optics != -1]
            no_clusters_optics = np.size(clusters_optics)
            for n in clusters_optics:
                current_cluster_idx = np.where(y_cluster == n)[0]
                cluster_center, r_2 = miniball.get_bounding_ball(
                    track[current_cluster_idx]
                )
                r = np.sqrt(r_2)
                cluster_populations_optics.append(np.size(current_cluster_idx))
                cluster_sizes_optics.append(r)
            # print(y_cluster)
            y_merged_cluster = y_cluster
            for index,element in enumerate(y_cluster):
                if index > 0:
                    if y_merged_cluster[index-1]!=-1 and y_merged_cluster[index]!=-1 and y_merged_cluster[index]!=y_merged_cluster[index-1]:
                        # print("I got there, my friend")
                        first_cluster = y_merged_cluster[index-1]
                        cluster_to_override = y_merged_cluster[index]
                        y_merged_cluster[y_merged_cluster == cluster_to_override] = first_cluster
            clusters_merged = np.unique(y_merged_cluster)
            clusters_merged = clusters_merged[clusters_merged != -1]
            no_clusters_merged = np.size(clusters_merged)
            for n in clusters_merged:
                current_cluster_idx = np.where(y_cluster == n)[0]
                cluster_center, r_2 = miniball.get_bounding_ball(
                    track[current_cluster_idx]
                )
                r = np.sqrt(r_2)
                cluster_populations_merged.append(np.size(current_cluster_idx))
                cluster_sizes_merged.append(r)
       # print(y_merged_cluster)


        # part for plotting
        if showPlot:
            # print(y_cluster)
            for no_cluster in range(0, no_clusters):
                x_pos = cluster_center_annotations[no_cluster]
                y_pos = threshold
                # y_pos = model.reachability_[int(np.rint(x_pos))]
                length = len(model.reachability_)
                x_pos = int(x_pos) / int(length)
                _, y_max = ax3.get_ylim()
                # y_pos=int(y_pos)/y_max
                y_pos = 0.4
                # print(y_pos)
                if np.mod(no_cluster, 2) == 0:
                    y_pos = 0.4
                else:
                    y_pos = 0.3
                # y_pos=2/no_clusters*0.25*no_cluster
                ax3.annotate(
                    "C "
                    + str(no_cluster + 1)
                    + "\n"
                    + "n="
                    + str(cluster_populations[no_cluster])
                    + "\n"
                    + str(cluster_sizes[no_cluster])
                    + " nm",
                    xy=(x_pos, y_pos),
                    xycoords="axes fraction",
                    ha="center",
                    va="bottom",
                    fontsize=13,
                )
            colors = colors_y_pred
            colors_merged =  m.to_rgba(y_merged_cluster)
            indeces = np.where(y_merged_cluster==-1)
            colors_merged[indeces] = [0.501, 0.501, 0.501, 1]
            # colors=colors_merged
            ax3.bar(space, model.reachability_, align="center", color=colors)
            ax3.tick_params(axis="both", which="both", labelsize=16)
            ax3.annotate(
                str(single_defects)
                + " single defects\n"
                + str(defects_in_clusters)
                + " defects in clusters\n"
                + str(no_clusters)
                + " clusters formed",
                xy=(0.4, 1),
                xytext=(12, -12),
                va="top",
                xycoords="axes fraction",
                textcoords="offset points",
                fontsize=16,
            )

            # ax3.hlines(threshold, 0, len(y_pred), label=str(threshold))
            ax3.axhline(threshold, label=str(threshold))
            print(ax3.hlines)
            labelLines(ax3.get_lines(), ha="left")
            ax3.set_ylabel("Distance [nm]", fontsize=16)
            ax3.set_xlabel("# of silicon recoil [-]", fontsize=16)
            ax4 = fig1.add_subplot(grid[idx_subfig + 1], projection="3d")
            ax4.scatter(track[:, 0], track[:, 1], track[:, 2], c=colors, s=12)
            # ax4.scatter(track_no_isolated[:, 0], track_no_isolated[:, 1], track_no_isolated[:, 2], c=colors_y_pred_no_isolated, s=12)
            # ax4.scatter(track_vacancies[:, 0], track_vacancies[:, 1],
            # track_vacancies[:, 2], c="white", s=12)
            ax4.tick_params(axis="both", which="both", labelsize=16)
            ax4.set_xlabel("X [nm]")
            ax4.set_ylabel("Y [nm]")
            ax4.set_zlabel("Z [nm]")
            idx_subfig += 2
            ax3.set_title(
                "min neighbours:  " + str(min_samples),
                fontdict={"fontsize": 13, "fontweight": "medium"},
            )
            fig1.tight_layout()
            # if automatic_color==True:
            #     plt.savefig(energy_figure+"_Si_colors_"+str(threshold))
            # else:
            #     plt.savefig(energy_figure+"_Si"+str(threshold))
            plt.show()
    return (
        cluster_populations,
        cluster_sizes,
        single_defects,
        defects_in_clusters,
        cluster_populations_optics,
        cluster_sizes_optics,
        no_clusters,
        no_clusters_optics,
        no_clusters_merged,
        cluster_populations_merged,
        cluster_sizes_merged
    )
