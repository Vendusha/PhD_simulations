import os
import pandas as pd
import shutil
# before running the script conda activate merging_env
# Function to process and merge two files
def process_and_merge_files(file_50um, file_inf):
    # Specify 'None' as a NaN value and replace all NaN values with 0 upon reading the file
    energies_50um = pd.read_csv(file_50um, delimiter='\t', header=None, usecols=[0], na_values='None').fillna(0).values.flatten()
    energies_inf = pd.read_csv(file_inf, delimiter='\t', header=None, usecols=[0], na_values='None').fillna(0).values.flatten()
    
    # Find the first energy value in '50um' that is also in 'inf'
    matching_energy = next((e for e in energies_50um if e in energies_inf), None)

    # Get the index of the matching energy in 'inf'
    if matching_energy is not None:
        cutoff_index = energies_inf.tolist().index(matching_energy)
    else:
        cutoff_index = -1

    # Read full content of the 'inf' file up to (but not including) the matching energy, handle 'None' as NaN
    if cutoff_index != -1:
        df_inf = pd.read_csv(file_inf, delimiter='\t', header=None, nrows=cutoff_index, na_values='None').fillna(0)
    else:
        df_inf = pd.read_csv(file_inf, delimiter='\t', header=None, na_values='None').fillna(0)

    # Read full content of the '50um' file starting from the matching energy, handle 'None' as NaN
    df_50um = pd.read_csv(file_50um, delimiter='\t', header=None, na_values='None').fillna(0)
    df_50um_from_matching = df_50um[df_50um[0] >= matching_energy] if matching_energy is not None else df_50um

    # Merge the two dataframes
    merged_df = pd.concat([df_inf, df_50um_from_matching])

    return merged_df

# Define parameters for folder names
energy_value = "21_eV"
um50_folder_names = ["50_um","100_um","300_um"]
Part_values = ["4_7_A","7_05_A","9_4_A"]
for um50_folder_name in um50_folder_names:
    for part_value in Part_values:
        # Directory paths for 'inf' and '50 um' folders based on parameters
        inf_folder_path = f'./Full_cascades_monolayer_corrected/SRIM_to_NIEL_{energy_value}_inf_{part_value}'  # Replace with the actual path to the 'inf' folder
        um50_folder_path = f'./Ions_{um50_folder_name}_monolayer_corrected/SRIM_to_NIEL_{energy_value}_{um50_folder_name}_{part_value}_part'  # Replace with the actual path to the '50 um' folder

        # Directory path for the '50_um_final' folder based on parameters
        final_folder_path = f'./Ions_{um50_folder_name}_monolayer_corrected/SRIM_to_NIEL_{energy_value}_{um50_folder_name}_{part_value}'  # Replace with the desired path for the final merged files
        # if final folder does not exist, create it
        if not os.path.exists(final_folder_path):
            os.makedirs(final_folder_path)
        # Iterate through files in 'inf' folder and merge corresponding '50 um' files
        for inf_file_name in os.listdir(inf_folder_path):
            if inf_file_name.endswith(".txt"):
                inf_file_path = os.path.join(inf_folder_path, inf_file_name)
                #um50_file_name = inf_file_name.replace("_inf", "_50_um")
                um50_file_path = os.path.join(um50_folder_path, inf_file_name)
                print(um50_file_path)
                
                if os.path.exists(um50_file_path):
                    # Process and merge the files
                    merged_df = process_and_merge_files(um50_file_path, inf_file_path)            
                    # Save the merged file in the '50_um_final' folder
                    output_file_path = os.path.join(final_folder_path, inf_file_name.replace("_inf", "_merged"))
                    merged_df.to_csv(output_file_path, sep='\t', index=False, header=False)
                    print(f"Merged and saved: {inf_file_name}")
                else:
                    shutil.copy(inf_file_path, final_folder_path)
                    print(final_folder_path)
                    print(f"Copied: {inf_file_name}")



        print("Merging process completed.")