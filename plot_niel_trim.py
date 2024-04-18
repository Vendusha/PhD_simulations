import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
import lindhardt as lind
plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"xtick.labelsize": 15})
plt.rcParams.update({"ytick.labelsize": 15})
plt.rcParams.update({"axes.labelsize": 15})
# def sNRT_dpa(i, td, tmax):
#     mSi   = 28.086
#     ZSi = 14.0
#     amu     = 1.6605E-24
#     mass_si = 28.086*amu
#     e = 1.60217663*10**(-19) 
#     def epsilon(Ep, a_TF, M_1, M_2, Z_1, Z_2, e):
#         return (Ep * a_TF * M_2/(Z_1 * Z_2 * e**2*(M_1 + M_2)))
    
#     def NNRT_dpa(Td, Ed):
#         if Td < Ed:
#             return 0
#         elif Td < 2*Ed/0.8 and Td >= Ed:
#             return 1
#         elif Td > 2*Ed/0.8:
#             return 0.8*Td/(2*Ed)
#         else:
#             print("Error errr errr")

#     def integrand(Td, Ed):
#         return NNRT_dpa(Td, Ed) * np.exp(-Td/td)

#     def Coulomb(t, TF, Tmax):
#         a = 1 # set a to 1 for simplicity
#         t1 = 0.5 * (t ** 0.5)
#         t3 = t1 * (Tmax ** 0.5)
#         F = np.log(1 + 1.1383 * t3) / (1 + 1.5863 * t3 + 0.05223 * (t3 ** 2))
#         result = (a ** 2) * (TF ** 2) * F / (2 * t1 * (Tmax ** 1.5))
#         return result

#     Ed = 1 # set the threshold energy to 1
#     integral, _ = quad(integrand, Ed, tmax, args=(Ed,))
#     result = i * integral * Coulomb(tmax)
#     return result

def plot_niel_trim(
    Energy_Silicon_linhard,
    NIEL_linhard_Si_classic,
    NIEL_energy,
    Energy,
    Total_vacancies,
    Ion_ionization,
    Ion_vacancy,
    Ion_phonon,
    Recoil_ionization,
    Recoil_vacancy,
    Recoil_phonon,
    Elements,
    Single_vacancies,
    Divacancies,
    Trivacancies,
    Tetra_vacancies,
    Penta_vacancies,
    vacancies_6_10,
    vacancies_11_50,
    vacancies_51_100,
    vacancies_101_200,
    vacancies_201_500,
    vacancies_above_500,
    Energy_vacancies
):  
#     fig = [plt.figure() for i in range(15)]
#     ax = [fig[i].add_subplot(111) for i in range(15)]
#     labels = ["Single vacancy", "Divacancy", "Trivacancy", "Tetravacancy", "Pentavacancy", "6-10","11-50", "51-100", "101-200", "201-500", "above 500"]
#     for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
# #        fig[i], ax[i] = plt.subplots()
#         for idx, array in enumerate([Single_vacancies, Divacancies, Trivacancies, Tetra_vacancies, Penta_vacancies, vacancies_6_10, vacancies_11_50, vacancies_51_100, vacancies_101_200, vacancies_201_500, vacancies_above_500]):
#             ax[i-1].plot(
#             Energy_vacancies[i] * 1000,
#             array[i] / Total_vacancies[i] * 1000,
#             "-x",
#             label=labels[idx],
#             )
#         ax[i-1].set_xlabel("Si-PKA Energy [keV]")
#         ax[i-1].set_xscale("log")
#         ax[i-1].set_ylabel("# of clusters with x vacancies per 1000 PKA")
#         ax[i-1].legend(fontsize = 8)
#         ax[i-1].set_title(str(Elements[i-1]))

#     fig_vac = [plt.figure() for i in range(11)]
#     ax_vac = [fig_vac[i].add_subplot(111) for i in range(11)]
#     labels = ["Single vacancy", "Divacancy", "Trivacancy", "Tetravacancy", "Pentavacancy", "6-10","11-50", "51-100", "101-200", "201-500", "above 500"]
#     for i, array in enumerate([Single_vacancies, Divacancies, Trivacancies, Tetra_vacancies, Penta_vacancies, vacancies_6_10, vacancies_11_50, vacancies_51_100, vacancies_101_200, vacancies_201_500, vacancies_above_500]):
#         for idx in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:       
#             ax_vac[i-1].plot(
#             Energy_vacancies[idx] * 1000,
#             array[idx] / Total_vacancies[idx] * 1000,
#             "-x",
#             label=Elements[idx-1],
#             )
#         ax_vac[i-1].set_xlabel("Si-PKA Energy [keV]")
#         ax_vac[i-1].set_xscale("log")
#         ax_vac[i-1].set_ylabel("# of clusters with x vacancies per 1000 PKA")
#         ax_vac[i-1].legend(fontsize = 8)
#         ax_vac[i-1].set_title(labels[i])
    fig1, ax1 = plt.subplots()
    ax1.plot(
    np.asarray(Energy_Silicon_linhard) * 1000,
    NIEL_linhard_Si_classic * 1000,
    "--o",
    color="cyan",
    label="Linhard model",
    )
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
        ax1.plot(
            np.asarray(Energy[i]) * 1000,
            NIEL_energy[i] * 1000,
            "--o",
            label="TRIM"+Elements[i-1],
        )

    ax1.set_ylabel("Energy deposited in NIEL [keV]", color="r")
    ax1.tick_params(axis="y", labelcolor="r")
    ax1.set_xlabel("Energy of incident Si-recoil [keV]")
    ax1.legend(loc="upper left", fontsize=8)
    plt.show()
    #ax1.set_xscale("log")
    #ax1.set_yscale("log")
    fig2, ax2 = plt.subplots()
    #fig2a, ax2a = plt.subplots()
    ##########plot of total vacancies from COLLISON file to total vacancies from GUI
    # for i in [1,2,3,4,5,6,7,8,9,10,11,12,13, 14, 15]:
    #     ax2.plot(
    #         np.asarray(Energy[i]) * 1000,
    #         np.asarray(Total_vacancies[i]),
    #         "-x",
    #         label="TRIM"+Elements[i-1],
    #     )
    NRT_DPA = []
    for energy in Energy[14]:
        NRT_DPA.append(lind.NRT_dpa(energy))
    ax2.plot(
        Energy[14] * 10 ** 3,
        NRT_DPA,
        "-o",
        color="c",
        label="NRT DPA Silicon",
    )
    ax2.set_xlabel("Energy [keV]")
    ax2.set_ylabel("Total vacancies")
    ax2.legend(loc="upper left")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
#     for i in range(len(Energy_Si)):
#         print(Total_vacancies_Si[i])
#     #print("Total vacancies from COLLISON file: ", Total_vacancies_Si)
#     print("Energy from COLLISON file: ", Energy_Si)
#     #################################
#     # ax2a.plot(
#     #     np.asarray(Energy_Si) * 1000,
#     #     (np.asarray(Total_vacancies_Si) / 1000-np.asarray(Total_vacancies_GUI))/(np.asarray(Total_vacancies_GUI))*100,
#     #     "-x",
#     #     color="b",
#     #     markersize=10,
#     #     label="Total vacancies from COLLISON file",
#     # )
#     # ax2a.set_xlabel("Energy [keV]")
#     # ax2a.set_ylabel("$(V_{collisionFile}-V_{GUI})/V_{GUI}$ [%]")
#     # #ax2.legend(loc="upper left")
    fig3, ax3 = plt.subplots()
    colors = ["b", "r", "g", "c","m", "y"]
    markers = ["-x", "-o", "-^", "-s"]
    for idx, i in enumerate([14,2, 12,15]):
        ax3.plot(
            np.asarray(Energy[i]) * 1000,
            np.asarray(Ion_ionization[i]),
            markers[idx],
            color = "b",
            label="TRIM"+Elements[i-1],
        )
        ax3.plot(
            np.asarray(Energy[i]) * 1000,
            np.asarray(Ion_vacancy[i]),
            markers[idx],
            color = "r",
            label="TRIM"+Elements[i-1],
        )
        ax3.plot(
            np.asarray(Energy[i]) * 1000,
            np.asarray(Ion_phonon[i]),
            markers[idx],
            color = "g",
            label="TRIM"+Elements[i-1],
        )
        ax3.plot(
            np.asarray(Energy[i]) * 1000,
            np.asarray(Recoil_ionization[i]),
            markers[idx],
            color = "c",
            label="TRIM"+Elements[i-1],
        )
        ax3.plot(
            np.asarray(Energy[i]) * 1000,
            np.asarray(Recoil_vacancy[i]),
            markers[idx],
            color = "m",
            label="TRIM"+Elements[i-1],
        )
        ax3.plot(
            np.asarray(Energy[i]) * 1000,
            np.asarray(Recoil_phonon[i]),
            markers[idx],
            color = "y",
            label="TRIM"+Elements[i-1],
        )

    markers_legend = ax3.legend(handles=[plt.Line2D([], [], marker='x', markersize=10, linestyle='None', color='black'),
                                        plt.Line2D([], [], marker='o', markersize=10, linestyle='None', color='black'),
                                        plt.Line2D([], [], marker='^', markersize=10, linestyle='None', color='black'),
                                        plt.Line2D([], [], marker='s', markersize=10, linestyle='None', color='black')],  
                                labels=['Si', 'He', 'Mg','P'], loc='upper right')

    colors_legend = ax3.legend(handles=[plt.Line2D([], [], linestyle='-', color='blue'),
                                        plt.Line2D([], [], linestyle='-', color='red'),
                                        plt.Line2D([], [], linestyle='-', color='green'),
                                        plt.Line2D([], [], linestyle='-', color='cyan'),
                                        plt.Line2D([], [], linestyle='-', color='magenta'),
                                        plt.Line2D([], [], linestyle='-', color='yellow')],
                                labels=['Ion-ionization', 'Ion-vacancy', 'Ion-phonon', 'Recoil-ionization',
                                        'Recoil-vacancy', 'Recoil-phonon'], loc='center right')

    # Add both legends to the plot
    ax3.add_artist(markers_legend)
    ax3.add_artist(colors_legend)

    # Set axis labels and scale
    ax3.set_xlabel("Energy [keV]")
    ax3.set_ylabel("Energy loss [%]")
    ax3.set_xscale("log")
    ax3.set_xlabel("Energy [keV]")
    ax3.set_ylabel("Energy loss [%]")
    ax3.set_xscale("log")
#     # ax2 = ax1.twinx()
#     #ax3.legend()
    fig4, ax4 = plt.subplots()
    ax4.set_title("NIEL/vac")
    # for i in [1,2,3,4,5,6,7,8,9,10,11,12,13, 14, 15]:
    #     ax4.plot(
    #         Energy[i] * 1000,
    #         NIEL_energy[i] * 10 ** 6 / Total_vacancies[i],
    #         "-x",
    #         label=Elements[i-1],
    #     )
    ax4.set_ylabel("NIEL/VAC [eV]")
    ax4.set_xlabel("Energy [keV]")
    #ax4.set_xlim(0, 50)
    ax4.legend()
    ax4.set_xscale("log")
    fig5, ax5 = plt.subplots()
    ax5.set_title("Total vacancies as a function of NIEL")
    # for i in [1,2,3,4,5,6,7,8,9,10,11,12,13, 14, 15]:
    #     ax5.plot(
    #         NIEL_energy[i] * 10 ** 6,
    #         Total_vacancies[i],
    #         "-o",
    #         label=Elements[i-1],
    #     )
    ax5.set_ylabel("Total vacancies")
    ax5.set_xlabel("NIEL [keV]")
    plt.legend()

    # ax2.set_ylabel("Total vacancies (TRIM)", color='b')
    # ax2.tick_params(axis='y', labelcolor='b')
    # ax1.legend(loc="upper left")
    # ax2.legend(loc="lower right")
    # fig2, ax2 = plt.subplots()
    # ax2.set_xlabel("Energy [keV]")
    # ax2.set_ylabel("NIEL [%]")
    # ax2.plot(np.asarray(Energy_Silicon_linhard)*1000, NIEL_linhard_Si_classic/(np.asarray(Energy_Silicon_linhard))*100, "--o", color = "red", label="Linhard model")
    # ax2.plot(np.asarray(Energy_Si)*1000,  np.asarray(Total_vacancies_Si)/(1000)*32/1000/1/(np.asarray(Energy_Si)*1000)*100, "-^", color = "b", label="Total vacancies/32 eV")
    # plt.legend()
    plt.show()
