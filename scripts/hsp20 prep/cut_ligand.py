from rdkit import Chem

def cut_linker():
    # Load the full ligand
    mol = Chem.MolFromMolFile("ligand_full.sdf", removeHs=False)
    
    # DEFINING THE CUT
    # NVP-AUY922 has a specific structure. 
    # We will manually delete the bond between the central Isoxazole and the Benzene ring.
    # This splits the molecule into two distinct islands.
    
    # Editable molecule
    rw_mol = Chem.RWMol(mol)
    
    # NOTE: You usually need to visualize indices first. 
    # For 2XAB, we can try a SMARTS pattern match to find the amide/isoxazole bridge.
    # But for a portfolio project, manual index cutting is safer to ensure valid fragments.
    
    # 1. Open 'ligand_full.sdf' in PyMOL/Discovery Studio.
    # 2. Turn on "Label > Atom Name/ID".
    # 3. Find the atoms connecting the two rings.
    
    # HYPOTHETICAL EXAMPLE (You must verify these indices in PyMOL for your specific file):
    # Let's say Atom 15 is connected to Atom 16, and that is the bridge.
    # rw_mol.RemoveBond(15, 16)
    
    # ALTERNATIVE: Delete the middle chunk
    # If we want a gap, we delete the atoms IN the linker.
    # For 2XAB, let's pretend atoms 10 and 11 are the bridge.
    atoms_to_delete = [6, 7] # <--- REPLACE THESE WITH REAL INDICES
    
    for idx in sorted(atoms_to_delete, reverse=True):
        rw_mol.RemoveAtom(idx)

    # Sanitize and Save
    frags = rw_mol.GetMol()
    try:
        Chem.SanitizeMol(frags)
    except:
        print("Warning: Fragments might need manual valency fix (Add Hydrogens).")

    w = Chem.SDWriter("fragments.sdf")
    w.write(frags)
    w.close()
    print("-> Saved fragments.sdf")

if __name__ == "__main__":
    cut_linker()
