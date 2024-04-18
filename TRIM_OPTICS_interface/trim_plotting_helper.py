import matplotlib.pyplot as plt
import numpy as np
import vispy.scene
from vispy.scene import visuals, SceneCanvas
plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"xtick.labelsize": 15})
plt.rcParams.update({"ytick.labelsize": 15})
plt.rcParams.update({"axes.labelsize": 15})
def plot_reachability(space, model, colors):
    fig = plt.figure(figsize=(20, 12))
    ax_reachability = fig.add_subplot()
    ax_reachability.bar(
        space, model.reachability_, align="center", color=colors
        )
                
    ax_reachability.set_ylabel("Distance [A]")
    ax_reachability.set_xlabel("Ordered silicon recoil [-]")
def plot_single_event(x,y,z, colors, total_count):
    color_dict = {'red': (1, 0, 0, 1),
              'green': (0, 1, 0, 1),
              'blue': (0, 0, 1, 1),
              'yellow': (1, 1, 0, 1),
              'black': (0, 0, 0, 1),
              'magenta': (1, 0, 1, 1),
              'cyan': (0, 1, 1, 1),
              'purple': (0.5, 0, 0.5, 1),
              'orange': (1, 0.5, 0, 1),
              'gray': (0.5, 0.5, 0.5, 1),
              'pink': (1, 0.75, 0.8, 1)}
    colors = np.array([color_dict[name] for name in colors])
    # Compute pairwise distances between all points
    dists = np.sqrt((x[:, np.newaxis] - x[np.newaxis, :])**2 +
                    (y[:, np.newaxis] - y[np.newaxis, :])**2 +
                    (z[:, np.newaxis] - z[np.newaxis, :])**2)  
    # Find the indices of the closest neighbors
    #ind = np.argmin(dists, axis=1)
    print("The total distance is:"+str(np.max(dists)/10000)) # in um
    print(np.sqrt(((x[-1]-x[0])**2+(y[-1]-y[0])**2+(z[-1]-z[0])**2))/10000)
    ind = np.argmin(np.where(dists==0, np.inf, dists), axis=1)
    #ind = np.argmin(dists, axis=1)
    canvas = SceneCanvas(keys='interactive', bgcolor='white', size=(1024, 768))
    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    #view.camera = 'arcball'
    #view.camera = "fly"
    #view.camera = 'panzoom'
    scatter = visuals.Markers()
    scatter.set_data(np.column_stack((x, y, z)), edge_color=None, face_color=colors, size=20)
    view.add(scatter)
    lines = visuals.Line(np.array([]), color='black')
    line_data = []
    for i, j in enumerate(ind):
        if i != j:
            line_data.append([x[i], y[i], z[i], x[j], y[j], z[j]])
    lines.set_data(np.array(line_data).reshape(-1, 2, 3))
    view.add(lines)
    
    for i, j in enumerate(ind):
        if i != j:
            dist = round(dists[i, j] * 10000, 2)
            label_pos = np.mean(np.array([[x[i], y[i], z[i]], [x[j], y[j], z[j]]]), axis=0)
            label = visuals.Text(str(dist), pos=label_pos, color='black', font_size=2000)
            view.add(label)
    xmin, ymin, zmin = np.min(np.column_stack((x, y, z)), axis=0)
    xmax, ymax, zmax = np.max(np.column_stack((x, y, z)), axis=0)
    view.camera.set_range(x=(xmin, xmax), y=(ymin, ymax), z=(zmin, zmax))
    view.camera.interactive = True
    axis = visuals.XYZAxis(width=2000, parent=view.scene)
    #grid_lines = visuals.GridLines(parent=view.scene, color="black")
    grid_lines = visuals.GridLines(color='gray', parent=view.scene, scale=(1, 1, 0))
    view.add(grid_lines)
    canvas.show()
    plt.show()
    vispy.app.run()
def plot(
    depth_samples, depth_samples_len, x_section, y_section, z_section, 
    color_origins_section, counts_section, isolated_count_section, section_count, 
    x_3D, y_3D, z_3D, color_array, total_count, count_inside_beam, count_outside_beam, 
    count_outside_beam_upper, count_outside_beam_lower, divacancy_section,trivacancy_section, folder, x_centers_section,
    y_centers_section, z_centers_section, population_section, energy_section, density_ball_section, density_cuboid_section,string_to_append):
    if string_to_append=="_before_correction":
        fig_3D = plt.figure()
        ax_3D = fig_3D.add_subplot(111, projection="3d")
        ax_3D.set_title("3d projection plot")
        ax_3D.set_xlabel("X [um]")
        ax_3D.set_ylabel("Y [um]")
        ax_3D.set_zlabel("Z [um]")
        fig_xy_origins = plt.figure()
        ax_xy_origins = fig_xy_origins.add_subplot(111)
        ax_xy_origins.set_title("Recoils Ion origin")
        ax_xy_origins.set_xlabel("X [um]")
        ax_xy_origins.set_ylabel("Y [um]")
        ax_3D.set_ylim(-2.5, 2.5)
        ax_3D.set_zlim(-2.5, 2.5)
        ax_3D.scatter(
        x_3D, y_3D, z_3D, c=color_array, s=0.2
        )  # time consuming, only do when really interested
        ax_3D.set_title(str(total_count) + " vacancies")
        ax_xy_origins.scatter(
            x_3D, y_3D, c=color_array, s=0.2
        )  # time consuming, only do when really interested
        ax_xy_origins.set_title(
            str(total_count)
            + " vacancies,\n"
            + str(count_inside_beam)
            + " inside beam,"
            + str(count_outside_beam)
            + " outside beam\n"
            + str(count_outside_beam_upper)
            + " in upper pixel, "
            + str(count_outside_beam_lower)
            + " in lower pixel"
        )
        ax_3D.figure.savefig(str(folder) + "/3D"+string_to_append+".png", dpi=800)
        ax_xy_origins.figure.savefig(str(folder) + "/xy_origins"+string_to_append+".png", dpi=800)
        plt.close(fig_3D)
        plt.close(fig_xy_origins)
    fig_yz_energy = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_energy = [
        fig_yz_energy[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    fig_yz_population = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_population = [
        fig_yz_population[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    fig_yz_density_ball = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_density_ball = [
        fig_yz_density_ball[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    fig_yz_density_cuboid = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_density_cuboid = [
        fig_yz_density_cuboid[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    fig_yz_origins_section = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_origins_section = [
        fig_yz_origins_section[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    fig_yz_xvacancy_section = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_xvacancy_section = [
        fig_yz_xvacancy_section[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    fig_yz_huhtinen_section = [plt.figure() for i in range(depth_samples_len)]
    ax_yz_huhtinen_section = [
        fig_yz_huhtinen_section[i].add_subplot(111) for i in range(depth_samples_len)
    ]
    for i in range(0, depth_samples_len):
        if len(x_section[i]) == 0:
            continue
        ax_yz_origins_section[i].set_xlabel("Y [um]")
        ax_yz_origins_section[i].set_ylabel("Z [um]")
        #ax_yz_origins_section[i].set_xlim(0, 1)
        #ax_yz_origins_section[i].set_ylim(0, 1)
        ax_yz_origins_section[i].scatter(
            y_section[i], z_section[i], c=color_origins_section[i], s=0.2
        )
        total_vacancies = 0
        for value, count in counts_section[i].items():
            total_vacancies += value * count
        total_vacancies += isolated_count_section[i]
        counts_bigger_than_3 = sum(filter(lambda x: x > 3, counts_section[i].values()))
        ax_yz_origins_section[i].set_title(
            str(isolated_count_section[i])
            + " isolated ,"
            + str(counts_section[i].get(2))
            + " divacancies ,\n"
            + str(counts_section[i].get(3))
            + " trivacancies ,"
            + str(counts_bigger_than_3)
            + " >3  "
            + str(total_vacancies)
        )

        ax_yz_origins_section[i].figure.savefig(
            str(folder)
            + "/2D_origins_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_origins_section[i].figure)

        ax_yz_huhtinen_section[i].set_xlabel("Y [um]")
        ax_yz_huhtinen_section[i].set_ylabel("Z [um]")
        #ax_yz_huhtinen_section[i].set_xlim(0, 1)
        #ax_yz_huhtinen_section[i].set_ylim(0, 1)
        ax_yz_huhtinen_section[i].scatter(y_section[i], z_section[i], c="black", s=0.2)
        ax_yz_huhtinen_section[i].set_title(str(section_count[i]) + " vacancies")
        ax_yz_huhtinen_section[i].figure.savefig(
            str(folder)
            + "/2D_huhtinen_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_huhtinen_section[i].figure)

        ax_yz_xvacancy_section[i].set_xlabel("Y [um]")
        ax_yz_xvacancy_section[i].set_ylabel("Z [um]")
        #ax_yz_xvacancy_section[i].set_xlim(0, 1)
        #ax_yz_xvacancy_section[i].set_ylim(0, 1)
        divacancy_section[i] = np.asarray(divacancy_section[i])
        trivacancy_section[i] = np.asarray(trivacancy_section[i])
        # print(divacancy_section[i])
        data = [
            (divacancy_section[i], "V2", ".", {"facecolors": "black", "edgecolors": "black"}),
            (trivacancy_section[i], "V3", "o", {"facecolors": "none", "edgecolors": "black"}),
        ]
        title_string = ""
        for section, label, marker, kwargs in data:
            if len(section) > 0:
                x, y, z = zip(*section)
                ax_yz_xvacancy_section[i].scatter(y, z, marker=marker, **kwargs)
                title_string+=(f"{len(x)} {label}  ")
        ax_yz_xvacancy_section[i].set_title(title_string)
        ax_yz_xvacancy_section[i].figure.savefig(
            str(folder)
            + "/2D_xvacancy_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_xvacancy_section[i].figure)

        ax_yz_energy[i].set_xlabel("Y [um]")
        ax_yz_energy[i].set_ylabel("Z [um]")
        #ax_yz_energy[i].set_xlim(0, 1)
        #ax_yz_energy[i].set_ylim(0, 1)
        sc = ax_yz_energy[i].scatter(
            y_centers_section[i], z_centers_section[i], c=energy_section[i], cmap='viridis', s=2
        )
        ax_yz_energy[i].set_title("Cluster energy")
        cbar = ax_yz_energy[i].figure.colorbar(sc, ax=ax_yz_energy[i])
        cbar.ax.set_ylabel("Energy [keV]", rotation=-90, va="bottom")
        ax_yz_energy[i].figure.savefig(
            str(folder)
            + "/2D_energy_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_energy[i].figure)

        ax_yz_population[i].set_xlabel("Y [um]")
        ax_yz_population[i].set_ylabel("Z [um]")
        #ax_yz_population[i].set_xlim(0, 1)
        #ax_yz_population[i].set_ylim(0, 1)
        sc = ax_yz_population[i].scatter(
            y_centers_section[i], z_centers_section[i], c=population_section[i], cmap='viridis', s=2
        )
        ax_yz_population[i].set_title("Cluster population")
        cbar = ax_yz_population[i].figure.colorbar(sc, ax=ax_yz_population[i])
        cbar.ax.set_ylabel("Population", rotation=-90, va="bottom")
        ax_yz_population[i].figure.savefig(
            str(folder)
            + "/2D_population_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_population[i].figure)

        ax_yz_density_ball[i].set_xlabel("Y [um]")
        ax_yz_density_ball[i].set_ylabel("Z [um]")
        #ax_yz_density_ball[i].set_xlim(0, 1)
        #ax_yz_density_ball[i].set_ylim(0, 1)
        sc = ax_yz_density_ball[i].scatter(
            y_centers_section[i], z_centers_section[i], c=density_ball_section[i], cmap='viridis', s=2
        )
        ax_yz_density_ball[i].set_title("Cluster density (sphere)")
        cbar = ax_yz_density_ball[i].figure.colorbar(sc, ax=ax_yz_density_ball[i])  
        cbar.ax.set_ylabel("Density [1/um$^3$]", rotation=-90, va="bottom")
        ax_yz_density_ball[i].figure.savefig(
            str(folder)
            + "/2D_density_ball_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_density_ball[i].figure)
        
        ax_yz_density_cuboid[i].set_xlabel("Y [um]")
        ax_yz_density_cuboid[i].set_ylabel("Z [um]")
        #ax_yz_density_cuboid[i].set_xlim(0, 1)
        #ax_yz_density_cuboid[i].set_ylim(0, 1)
        sc = ax_yz_density_cuboid[i].scatter(
            y_centers_section[i], z_centers_section[i], c=density_cuboid_section[i], cmap='viridis', s=2
        )
        ax_yz_density_cuboid[i].set_title("Cluster density (cuboid)")
        cbar = ax_yz_density_cuboid[i].figure.colorbar(sc, ax=ax_yz_density_cuboid[i])
        cbar.ax.set_ylabel("Density [1/um$^3$]", rotation=-90, va="bottom")
        ax_yz_density_cuboid[i].figure.savefig(
            str(folder)
            + "/2D_density_cuboid_section_"
            + str(depth_samples[i])
            + "-"
            + str(depth_samples[i] + 1)
            +string_to_append
            + ".png",
            dpi=600,
        )
        plt.close(ax_yz_density_cuboid[i].figure)

    return 