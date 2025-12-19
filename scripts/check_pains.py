import glob
import os
from rdkit import Chem
from rdkit.Chem import FilterCatalog

# 1. Initialize the PAINS filter
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog.FilterCatalog(params)

# Setup output
input_folder = "linker_sdf"
output_csv = "pains_results.csv"

print(f"{'Filename':<45} | {'PAINS?':<8} | {'Reason'}")
print("-" * 80)

with open(output_csv, "w") as out:
    out.write("Filename,Has_PAINS,Reason\n")
    
    sdf_files = glob.glob(os.path.join(input_folder, "*.sdf"))
    
    if not sdf_files:
        print(f"No SDF files found in {input_folder}")
    
    for f in sdf_files:
        mol = Chem.MolFromMolFile(f)
        
        if mol:
            # Check for matches
            if catalog.HasMatch(mol):
                # Get the specific PAINS family name
                entry = catalog.GetFirstMatch(mol)
                reason = entry.GetDescription()
                status = "YES"
                print(f"{f:<45} | \033[91m{status:<8}\033[0m | {reason}") # Red text for danger
            else:
                status = "No"
                reason = "-"
                print(f"{f:<45} | {status:<8} | {reason}")
                
            out.write(f"{f},{status},{reason}\n")
        else:
            print(f"{f:<45} | Error    | Invalid Molecule")

print("-" * 80)
print(f"Results saved to {output_csv}")
