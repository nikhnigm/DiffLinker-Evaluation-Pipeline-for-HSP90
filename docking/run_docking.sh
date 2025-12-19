#!/bin/bash

# --- CONFIGURATION ---
RECEPTOR="receptor.pdbqt"
CENTER_X=-1.0
CENTER_Y=35.3
CENTER_Z=-0.8

# Search Box Size (20 Angstroms is standard)
SIZE_X=20
SIZE_Y=20
SIZE_Z=20

EXHAUSTIVENESS=4
# ---------------------

mkdir -p docking_results
echo "Filename,Binding_Affinity" > docking_scores.csv

# Verify receptor exists
if [ ! -f "$RECEPTOR" ]; then
    echo "Error: $RECEPTOR not found!"
    exit 1
fi

echo "Starting docking..."

for ligand in pdbqt_files/*.pdbqt; do
    if [ -f "$ligand" ]; then
        base_name=$(basename "$ligand" .pdbqt)
        echo "Docking $base_name..."
        
        vina \
            --receptor $RECEPTOR \
            --ligand "$ligand" \
            --center_x $CENTER_X \
            --center_y $CENTER_Y \
            --center_z $CENTER_Z \
            --size_x $SIZE_X \
            --size_y $SIZE_Y \
            --size_z $SIZE_Z \
            --exhaustiveness $EXHAUSTIVENESS \
            --out "docking_results/${base_name}_out.pdbqt" \
            --cpu 4 > "docking_results/${base_name}.log"
        
        # Get the score
        score=$(grep "   1 " "docking_results/${base_name}.log" | head -n 1 | awk '{print $2}')
        echo "${base_name},${score}" >> docking_scores.csv
    fi
done

echo "Done! Scores saved to docking_scores.csv"
