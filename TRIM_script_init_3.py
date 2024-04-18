import os
import shutil
import subprocess
import glob
import re
import pyautogui
import time
import sys 

def change_energy(file_path, energy, z):
    mass = [1.008, 4.003, 7.016, 9.012, 11.009,12, 14.003, 15.995, 18.998, 19.992, 22.99, 23.985, 26.982, 27.977, 30.974]
    line_rewriting = False
    line_rewriting_energy = False
    lattice_rewriting = False
    target_elements_rewriting = False
    ion_rewriting = False
    surface_rewriting = False
    with open(file_path, "r") as f:
        lines = f.readlines()    
    for i, line in enumerate(lines):
        parts = line.split()
        if ion_rewriting:
           parts[0] = f"{int(z)}"
           parts[1] = f"{float(mass[z-1])}"
           ion_rewriting = False
        if target_elements_rewriting:
            parts[0] = "Atom 1 = Si =       14  28.086"
            parts = [parts[0]]
            target_elements_rewriting = False
        if lattice_rewriting:
            parts[0] = f"{int(2)}"
            lattice_rewriting = False
        if surface_rewriting:
            parts[0] = f"{int(21)}"
            surface_rewriting = False
        if line_rewriting_energy:
                parts[2] = f"{energy:>11.3f}"
                if energy <5:
                    parts[4] = f"{float(10000):>11.3f}"
                    parts[6] = f"{float(10000):>11.3f}"
                elif energy <1000:
                    parts[4] = f"{float(1000):>11.3f}"
                    parts[6] = f"{float(1000):>11.3f}"
                else:
                    parts[4] = f"{float(300):>11.3f}"
                    parts[6] = f"{float(300):>11.3f}" 
                line_rewriting_energy = False
        if line_rewriting:
            parts[0] = f"{int(3)}"
            line_rewriting = False
        if "Ion: Z1" in line:
            ion_rewriting = True
        if "Energy (keV)" in line:
            line_rewriting_energy = True
        if "Cascades" in line:
            line_rewriting = True
        if "lattice binding energies" in line:
            lattice_rewriting = True
        if "surface binding energies" in line:
            surface_rewriting = True
        if "Target Elements:" in line:
            target_elements_rewriting = True
        lines[i] = " ".join(parts) + "\r\n"
    
    with open(file_path, "w", encoding="cp850") as f:
        f.writelines(lines)

#path = '/media/dulinka/SKA_clustered_is/G4_results/RD50_TRIM_folder/Helium_finish'
#trim_dir = '/home/dulinka/SRIM_TRIM/'
path = sys.argv[1]
z_number = int(sys.argv[2])
trim_dir = "/home/dulinka/SRIM_TRIM_2/"
os.chdir("/home/dulinka/SRIM_TRIM_2/")
for folder in sorted(os.listdir(path), key=lambda x: (float(re.search(r'\d+\.\d+|\d+', x).group()) if re.search(r'\d+\.\d+|\d+', x) else 0, x)):
    item_path = os.path.join(path, folder)
    energy =float(folder.split('_')[0])*1000
    # if energy <25:
    #     continue
    for file in os.listdir(item_path):
        if file.endswith(".dat"):
            dat_file = os.path.join(item_path, file)
            break
    #shutil.copy2(dat_file, os.path.join(trim_dir,"TRIM.DAT"))
    change_energy(os.path.join(trim_dir, 'TRIM.IN'), energy, z_number)
    result = subprocess.Popen(["wine", "TRIM.exe"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    title = "End of TRIM.DAT calculation"
    print(energy)
    while True:
        time.sleep(4)
        screenshot = pyautogui.screenshot()
        screenshot.save(os.path.join(item_path, str(energy/1000)+'.png'))  
        time.sleep(4)
        if not(result.poll() is None): #if the process ended
            break   
    for file in ['E2RECOIL.txt', 'IONIZ.txt', 'LATERAL.txt', 'NOVAC.txt', 'PHONON.txt', 'RANGE.txt', 'TDATA.txt', 'VACANCY.txt']:
        shutil.copy2(os.path.join(trim_dir, file), item_path)
    for file in ['COLLISON.txt', 'SPUTTER.txt', 'TRANSMIT.txt', 'TRANSREC.txt', 'TRIMOUT.txt']:
        shutil.copy2(os.path.join(trim_dir, 'SRIM Outputs', file), item_path)