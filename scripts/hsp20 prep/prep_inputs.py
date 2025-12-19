from biopandas.pdb import PandasPdb
from rdkit import Chem
import sys
import os

# CONFIGURATION
PDB_FILE = "2xab.pdb"
LIGAND_CODE = "VHD"  # Confirmed from your debug steps
CHAIN_ID = "A"       # Isolate only one side of the dimer

def prep_structure():
    print(f"--- Processing {PDB_FILE} ---")
    
    # Check if file exists
    if not os.path.exists(PDB_FILE):
        print(f"Error: {PDB_FILE} not found in current directory.")
        return

    # Load the PDB structure
    try:
        ppdb = PandasPdb().read_pdb(PDB_FILE)
    except Exception as e:
        print(f"Error reading PDB file: {e}")
        return

    # --- 1. SAVE CLEAN PROTEIN (Chain A Only) ---
    print(f"Isolating Protein Chain {CHAIN_ID}...")
    
    # Filter ATOMs for Chain A
    protein = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == CHAIN_ID]
    ppdb.df['ATOM'] = protein
    
    # Clear HETATMs (Ligands/Water) for the protein file
    # We assign a dummy DataFrame to clear it out safely
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] == 'placeholder'] 
    
    ppdb.to_pdb(path="protein_clean.pdb", records=['ATOM'], gz=False)
    print("-> Saved 'protein_clean.pdb' (Receptor)")

    # --- 2. EXTRACT LIGAND (Chain A Only) ---
    print(f"Extracting Ligand {LIGAND_CODE} from Chain {CHAIN_ID}...")
    
    # Reload PDB to get full data back
    ppdb = PandasPdb().read_pdb(PDB_FILE)
    
    # Filter for Ligand Code AND Chain ID
    ligand_df = ppdb.df['HETATM'][
        (ppdb.df['HETATM']['residue_name'] == LIGAND_CODE) & 
        (ppdb.df['HETATM']['chain_id'] == CHAIN_ID)
    ]
    
    # Sanity Check
    if len(ligand_df) == 0:
        print(f"\n[ERROR] Could not find ligand '{LIGAND_CODE}' in Chain '{CHAIN_ID}'.")
        print("Debug Info - Found these residues in HETATM:")
        print(ppdb.df['HETATM']['residue_name'].unique())
        return

    # Save temporarily as PDB so RDKit can read it
    temp_lig_pdb = "ligand_temp.pdb"
    ligand_pdb = PandasPdb()
    ligand_pdb.df['HETATM'] = ligand_df
    # Add dummy ATOM record to satisfy BioPandas structure requirements
    ligand_pdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == 'XXX'] 
    
    ligand_pdb.to_pdb(path=temp_lig_pdb, records=['HETATM'], gz=False)

    # --- 3. CONVERT TO SDF ---
    # removeHs=False is critical to keep the chemistry valid
    mol = Chem.MolFromPDBFile(temp_lig_pdb, removeHs=False)
    
    if mol:
        print(f"-> Ligand loaded successfully: {mol.GetNumAtoms()} atoms")
        writer = Chem.SDWriter("ligand_full.sdf")
        writer.write(mol)
        writer.close()
        print("-> Saved 'ligand_full.sdf' (Ground Truth)")
        
        # Cleanup temp file
        if os.path.exists(temp_lig_pdb):
            os.remove(temp_lig_pdb)
    else:
        print("Error: RDKit could not parse the extracted PDB.")

if __name__ == "__main__":
    prep_structure()
