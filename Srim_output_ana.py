import numpy as np
e_binding = 2
e_init = 500
depth, vac_ion, vac_recoil = np.loadtxt(open("VACANCY.txt",'rt').readlines()[:-1], skiprows=34, unpack=True)
depth, ioniz_ion, ioniz_recoil = np.loadtxt("IONIZ.txt", skiprows=25, unpack=True)
depth, phonon_ion, phonon_recoil = np.loadtxt("PHONON.txt",skiprows=23, unpack=True)
depth, e2rec_ion, e2rec_recoil = np.loadtxt("E2RECOIL.txt", skiprows=24, unpack=True)
depth_width = depth[1]-depth[0]
vac_ion_tot = np.sum(vac_ion)*depth_width*2
vac_recoil_tot = np.sum(vac_recoil)*depth_width*2
ioniz_ion_tot = np.sum(ioniz_ion)*depth_width
ioniz_recoil_tot = np.sum(ioniz_recoil)*depth_width
phonon_ion_tot = np.sum(phonon_ion)*depth_width
phonon_recoil_tot = np.sum(phonon_recoil)*depth_width
e2rec_ion_tot = np.sum(e2rec_ion)*depth_width
e2rec_recoil_tot = np.sum(e2rec_recoil)*depth_width
tot_loss = vac_ion_tot + ioniz_ion_tot + phonon_ion_tot +vac_recoil_tot + ioniz_recoil_tot + phonon_recoil_tot
I1 = ioniz_ion_tot * 100 / tot_loss
I2 = ioniz_recoil_tot * 100 / tot_loss
V1 = vac_ion_tot * 100 / tot_loss
V2 = vac_recoil_tot * 100 / tot_loss
P1 = phonon_ion_tot * 100 / tot_loss
P2 = phonon_recoil_tot * 100 / tot_loss
print("╔════════════════════════════╦═════════════════════════╗")
print("║           PERCENTAGE LOSS  ║     ENERGY LOSS (keV)   ║")
print("║           IONS    RECOILS  ║    IONS      RECOILS    ║")
print("║          ───────  ───────  ║  ────────── ──────────  ║")
print("║ Ioniz =   {:.3f}    {:.3f}  ║  {:.3f} {:.3f}  ║".format(I1, I2, ioniz_ion_tot/1000, ioniz_recoil_tot/1000))
print("║ Vac.  =   {:.3f}     {:.3f}  ║  {:.3f} {:.3f}  ║".format(V1, V2, vac_ion_tot/1000, vac_recoil_tot/1000))
print("║ Phonon=   {:.3f}   {:.3f}  ║  {:.3f} {:.3f}  ║".format(P1, P2, phonon_ion_tot/1000, phonon_recoil_tot/1000))
TOTAL1 = (ioniz_ion_tot + vac_ion_tot + phonon_ion_tot) * 100 / tot_loss
TOTAL2 = (ioniz_recoil_tot + vac_recoil_tot + phonon_recoil_tot) * 100 / tot_loss
TOTAL3 = (ioniz_ion_tot + vac_ion_tot + phonon_ion_tot) / 1000
TOTAL4 = (ioniz_recoil_tot + vac_recoil_tot + phonon_recoil_tot) / 1000
print("║ Totals=  {:.3f} {:.3f}     ║  {:.3f} {:.3f}  ║".format(TOTAL1, TOTAL2, TOTAL3, TOTAL4))
print("╚════════════════════════════╩═════════════════════════╝")
print("        ╔═══════════════════════════════════════╗")
print("        ║   TOTAL ENERGY LOSS = {:.3f} keV  ║".format(tot_loss/1000))
print("        ╚═══════════════════════════════════════╝")
VAC1 = int(V1/e_binding)
VAC2 = int(V2/e_binding)
print("╔═════════════════════════════════════════════════════╗")
print("║  Energy Transferred to Recoils = {:.3f} keV      ║".format(e2rec_ion_tot/1000))
print("║  Energy Received by Recoils    = {:.3f} keV      ║".format(e2rec_recoil_tot/1000))
print("║  Ion Vacancies = {:f} , Recoil Vacancies = {:f} ║".format(VAC1, VAC2))
print("║  Total Vacancies/Ion = {:f}                     ║".format(VAC1 + VAC2))
print("╚═════════════════════════════════════════════════════╝")

print("NOTE : Totals may be in error if final ion is not yet stopped,")
print("       or if Ions or Recoils leave the Plotting Window.")
print("NOTE : Summation currently valid only for single element targets.")