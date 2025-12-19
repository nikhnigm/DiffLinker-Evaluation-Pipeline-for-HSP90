import glob
import os
from rdkit import Chem
from rdkit.Chem import RDConfig
import sys

# Add the current directory to python path so we can import sascorer
sys.path.append(os.getcwd())

try:
    import sascorer
except ImportError:
    print("Error: Could not import 'sascorer'.")  #Download from https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score/sascorer.py
    sys.exit(1)

# Setup output file
output_csv = "sa_scores.csv"

# Adjusted width to 45 to fit the folder path in the filename column
print(f"{'Filename':<45} | {'SA Score':<8} | {'Status'}")
print("-" * 75) 

with open(output_csv, "w") as out:
    out.write("Filename,SA_Score\n")  
    
    # Define the folder containing the SDF files
    input_folder = "linker_sdf"
    
    # Process all SDF files in the 'linker_sdf' folder
    search_path = os.path.join(input_folder, "*.sdf")
    sdf_files = glob.glob(search_path)
    
    # Check if files were actually found
    if not sdf_files:
        print("Please ensure the folder exists and contains your files.")
    
    for f in sdf_files:
        mol = Chem.MolFromMolFile(f)
        
        if mol:
            try:
                # Calculate SA Score
                # Ranges from 1 (easy) to 10 (hard)
                sa = sascorer.calculateScore(mol)  
                
                # Interpretation
                if sa < 3.5:
                    status = "Easy"
                elif sa < 6.0:
                    status = "Medium"
                else:
                    status = "Hard"
                    
                print(f"{f:<45} | {sa:.2f}     | {status}")
                out.write(f"{f},{sa:.2f}\n")
                
            except Exception as e:
                print(f"{f:<45} | Error    | {e}")
        else:
            print(f"{f:<45} | Invalid  | -")

print("-" * 75)
print(f"Results saved to {output_csv}")
