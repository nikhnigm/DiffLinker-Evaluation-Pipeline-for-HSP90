import pandas as pd
import glob
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, RDLogger

# 1. Setup & Environment
# Suppress RDKit warnings to keep the terminal output clean
RDLogger.DisableLog('rdApp.*') 

try:
    # Load docking results to calculate Ligand Efficiency (LE)
    scores = pd.read_csv("docking_scores.csv")
    # Ensure ID matches the SDF filename without extensions or paths
    scores['ID'] = scores['Filename'].apply(lambda x: os.path.basename(x).replace('.pdbqt', '').replace('.log', ''))
except Exception as e:
    print(f"Warning: docking_scores.csv not found or error parsing: {e}. LE will be 0.")
    scores = pd.DataFrame(columns=['ID', 'Binding_Affinity'])

print(f"{'ID':<30} | {'Valid?':<6} | {'Unique?':<7} | {'Lipinski':<8} | {'QED':<6} | {'LE':<5}")
print("-" * 85)

# 2. Processing Pipeline
results = []
seen_inchis = set()
# Search for all SDF files in the current working directory
sdf_files = glob.glob("*.sdf")

for f in sdf_files:
    filename = os.path.basename(f).replace('.sdf', '')
    mol = Chem.MolFromMolFile(f)
    
    # Tier 1: Validity (The Gatekeeper)
    if not mol:
        print(f"{filename:<30} | Fail   | -       | -        | -    | -")
        continue
    
    try:
        Chem.SanitizeMol(mol)
        is_valid = "Yes"
    except:
        print(f"{filename:<30} | No     | -       | -        | -    | -")
        continue

    # Tier 2: Uniqueness (Diversity Check)
    inchi = Chem.MolToInchiKey(mol)
    is_unique = "Yes" if inchi not in seen_inchis else "No"
    seen_inchis.add(inchi)

    # Tier 3: Drug-likeness (The Pill Check)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    qed_val = round(QED.qed(mol), 2)
    
    # Lipinski's Rule of 5: Max 1 violation allowed
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    lip_status = "Pass" if violations <= 1 else f"Fail({violations})"

    # Tier 4: Ligand Efficiency (Potency per Atom)
    match = scores[scores['ID'] == filename]
    le = 0.0
    if not match.empty:
        affinity = float(match.iloc[0]['Binding_Affinity'])
        heavy_atoms = mol.GetNumHeavyAtoms()
        le = round((-1 * affinity) / heavy_atoms, 2) if heavy_atoms > 0 else 0

    # Output to Console
    print(f"{filename:<30} | {is_valid:<6} | {is_unique:<7} | {lip_status:<8} | {qed_val:<6} | {le:<5}")
    
    # Store for CSV Report
    results.append({
        'ID': filename,
        'Valid': is_valid,
        'Unique': is_unique,
        'Lipinski_Violations': violations,
        'QED': qed_val,
        'LE': le,
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'Verdict': 'Keep' if is_unique == 'Yes' and violations <= 1 and is_valid == 'Yes' else 'Drop'
    })

# 3. Final Export
df = pd.DataFrame(results)
df.to_csv("comprehensive_analysis.csv", index=False)
print("-" * 85)
print(f"Done. Processed {len(sdf_files)} molecules. Final report: comprehensive_analysis.csv")