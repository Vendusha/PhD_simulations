#!/bin/bash
for folder in Sodium_to_do/*; do
     if [ -d "$folder" ]; then
        echo "$folder"
        python plot_trim.py "$folder" 2 4.7 "dbscan"
    fi
done && for folder in Magnesium/*; do
     if [ -d "$folder" ]; then
        echo "$folder"
        python plot_trim.py "$folder" 2 4.7 "dbscan"
    fi
done && for folder in Aluminium/*; do
     if [ -d "$folder" ]; then
        echo "$folder"
        python plot_trim.py "$folder" 2 4.7 "dbscan"
    fi
done && for folder in Silicon/*; do
     if [ -d "$folder" ]; then
        echo "$folder"
        python plot_trim.py "$folder" 2 4.7 "dbscan"
    fi
done && for folder in Phosphorus/*; do
     if [ -d "$folder" ]; then
        echo "$folder"
        python plot_trim.py "$folder" 2 4.7 "dbscan"
    fi
done 