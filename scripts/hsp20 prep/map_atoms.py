from rdkit import Chem
from rdkit.Chem import Draw

def map_atom_indices():
    # Load the full ligand you just made
    mol = Chem.MolFromMolFile("ligand_full.sdf", removeHs=False)
    
    # Label every atom with its Index ID
    for atom in mol.GetAtoms():
        atom.SetProp('atomLabel', str(atom.GetIdx()))

    # Draw and save the map
    Draw.MolToFile(mol, 'ligand_map.png', size=(800, 600))
    print("-> Generated 'ligand_map.png'. Open this image to see atom numbers.")

if __name__ == "__main__":
    map_atom_indices()
