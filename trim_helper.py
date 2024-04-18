import numpy as np
import copy
import load_trim
from collections import Counter
def generate_random_tracks(
    n_tracks,
    PKA,
    depth_samples,
    y_min_window,
    y_max_window,
    z_min_window,
    z_max_window,
    x_section,
    y_section,
    z_section,
    color_origins_section,
    divacancy_section,
    trivacancy_section,
    tetravacancy_section,
    pentavacancy_section,
    isolated_count_section,
    color_cluster_section,
    section_count,
    counts_section,
    colors,
    n,
    xi,
    method,
    x_centers_section,
    y_centers_section,
    z_centers_section,
    population_section,
    energy_section,
    density_ball_section,
    density_cuboid_section,
    no_inc
    ):
    n_tracks = int((len(PKA)/no_inc)*n_tracks)
    ion_z = np.random.uniform(-5, 5, size=n_tracks+1)
    ion_y = np.random.uniform(-5, 5, size=n_tracks+1)
    # ion_z = np.concatenate(
    #     [
    #         np.random.uniform(-5, 0, size=int(n_tracks/2)),
    #         np.random.uniform(1, 6, size=int(n_tracks/2)),
    #     ]
    # )
    # ion_y = np.concatenate(
    #     [
    #         np.random.uniform(-5, 0, size=int(n_tracks/2)),
    #         np.random.uniform(1, 6, size=int(n_tracks/2)),
    #     ]
    # )
    ion_x = np.random.uniform(-0.5, 0.5, size=n_tracks + 1)
    for i in range(n_tracks):
        index = i % len(PKA)
        pka = PKA[index]
        if i % 100000 == 0:
            print(str(i) + "/" + str(n_tracks) + " additional events generated.")
        for j in range(0, len(depth_samples)):
            x_min_window = depth_samples[j]
            x_max_window = depth_samples[j] + 1
            x_min = max(pka.min_x+ion_x[i], x_min_window)
            x_max = min(pka.max_x+ion_x[i], x_max_window)
            y_min = max(pka.min_y - pka.ion_position[1] + ion_y[i], y_min_window)
            y_max = min(pka.max_y - pka.ion_position[1] + ion_y[i], y_max_window)
            z_min = max(pka.min_z - pka.ion_position[2] + ion_z[i], z_min_window)
            z_max = min(pka.max_z - pka.ion_position[2] + ion_z[i], z_max_window)

            if x_min > x_max or y_min > y_max or z_min > z_max:
                continue
            else:
                track = pka.recoil_positions
                track_section = []
                x = np.array([row[0] for row in track])+ion_x[i]
                y = np.array([row[1] for row in track]) - pka.ion_position[1] + ion_y[i]
                z = np.array([row[2] for row in track]) - pka.ion_position[2] + ion_z[i]
                x_indices = np.where((x >= x_min_window) & (x <= x_max_window))[0]
                y_indices = np.where((y >= y_min_window) & (y <= y_max_window))[0]
                z_indices = np.where((z >= z_min_window) & (z <= z_max_window))[0]
                indices = np.intersect1d(
                    np.intersect1d(x_indices, y_indices), z_indices
                )
                if indices.size == 0:
                    continue
                x_add, y_add, z_add = zip(*track)
                x_section[j].extend(np.asarray(x_add))
                y_section[j].extend(np.asarray(y_add))
                z_section[j].extend(np.asarray(z_add))
                track_section = np.asarray([x[indices], y[indices], z[indices]]).T
                pka_section = load_trim.RecoilsFrom1PKA()
                pka_section.add_origins(
                    pka.primary_event, pka.ion_origin, pka.ion_energy, pka.ion_position
                )
                indices = np.asarray(indices).ravel()
                color_origins_section[j].extend([colors[pka.ion_origin - 2]] * len(x_add))
                pka_section.recoil_energies = np.take(pka.recoil_energies, indices)
                pka_section.recoil_positions = track_section
                pka_section.cluster_analysis(
                    track_section[0], n, xi, method, showPlot=False, cutoff_1=0, cutoff_2=1
                )
                if len(pka_section.divacancy) > 0:
                    divacancy_section[j].extend(pka_section.divacancy)
                if len(pka_section.trivacancy) > 0:
                    trivacancy_section[j].extend(pka_section.trivacancy)
                if len(pka_section.tetravacancy) > 0:
                    tetravacancy_section[j].extend(pka_section.tetravacancy)
                if len(pka_section.pentavacancy) > 0:
                    pentavacancy_section[j].extend(pka_section.pentavacancy)
                isolated_count_section[j] += pka_section.defects_isolated
                # color_cluster_section[i].append(pka.cluster_colors[i] for i in indices) good for later when doing the whole analysis
                color_cluster_section[j].extend(pka_section.cluster_colors)
                section_count[j] += indices.size
                counts_section[j] += Counter(pka_section.cluster_populations)
                if len(pka_section.cluster_populations) > 0:
                    x_centers, y_centers, z_centers = zip(*(pka_section.cluster_centers))               
                    x_centers = np.asarray(x_centers)
                    y_centers = np.asarray(y_centers)
                    z_centers = np.asarray(z_centers)
                    x_centers_section[j].extend(x_centers.ravel())
                    y_centers_section[j].extend(y_centers.ravel())
                    z_centers_section[j].extend(z_centers.ravel())
                    population_section[j].extend(pka_section.cluster_populations)
                    energy_section[j].extend(np.asarray(pka_section.cluster_energies)/1000)
                    density_ball_section[j].extend(pka_section.cluster_densities_sphere)
                    density_cuboid_section[j].extend(pka_section.cluster_densities_cuboid)

    return (
        x_section,
        y_section,
        z_section,
        divacancy_section,
        trivacancy_section,
        tetravacancy_section,
        pentavacancy_section,
        color_origins_section,
        color_cluster_section,
        section_count,
        counts_section,
        x_centers_section,
        y_centers_section,
        z_centers_section,
        population_section,
        energy_section,
        density_ball_section,
        density_cuboid_section
    )
