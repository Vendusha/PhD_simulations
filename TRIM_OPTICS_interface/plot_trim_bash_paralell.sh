#!/bin/bash
source ~/dgcode_projects/bootstrap.sh
# conda activate dgcode && cd ~/dbverify && sb && cd -
# Define the function to process each folder with a given threshold
process_folder() {
    folder=$1
    threshold=$2
    echo "Processing folder: $folder with threshold: $threshold"
    python plot_trim.py "$folder" 2 "$threshold" "dbscan"
}

export -f process_folder

# Define the parent directories and thresholds
base_path=".."  # Adjust this as needed to be absolute if running from different location
parents=("Hydrogen" "Helium" "Lithium" "Beryllium" "Boron" "Carbon" "Nitrogen" "Oxygen" "Fluorine" "Neon" "Sodium" "Magnesium" "Aluminium" "Silicon" "Phosphorus")
thresholds=(2.4 4.7 7.05 9.4)

# Generate and pass each directory to parallel
for parent in "${parents[@]}"; do
    subdirs=($base_path/$parent/*)  # Expands to all subdirectories within each parent
    for subdir in "${subdirs[@]}"; do
        if [ -d "$subdir" ]; then  # Check if it is a directory
            parallel process_folder ::: "$subdir" ::: "${thresholds[@]}"
        fi
    done
done
