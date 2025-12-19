# DiffLinker Evaluation Pipeline for HSP90 : De Novo Linker Generation and Evaluation
A deep learning-based drug discovery pipeline that utilizes DiffLinker (Diffusion-based conditional generative model) to design molecular linkers. The project features a rigorous "Surgical Cut" benchmark on HSP90 and a multi-tier evaluation framework integrating AutoDock Vina and RDKit.

### Project Overview
This project demonstrates the ability to "recover" a known drug molecule by breaking it into fragments and tasking an AI to generate chemically valid, high-affinity linkers. We specifically target HSP90 (Heat Shock Protein 90) using the crystal structure *PDB: 2XAB*.

---
### Phase 1 : Data Preparation (Surgical Cut)

1. To benchmark the model, we established a "Ground Truth" recovery mission.
2. Target Selection: Isolated Chain A of HSP90 from PDB: 2XAB.
3. Ligand Extraction: Extracted the inhibitor NVP-AUY922 and cleaned the environment (removed water/ non-functional ions).
4. The Cut: Using cut_ligand.py, we identified and deleted the Carbonyl bridge (C=O), splitting the drug into two distinct islands: the Resorcinol ring and the Isoxazole ring.

// image

### Phase 2: AI Generation (DiffLinker)

1. We deployed DiffLinker in a restricted CPU environment (WSL2) to generate 50 potential linker candidates.
2. Infrastructure Highlights:
Environment: WSL2 (Ubuntu) with a Python 3.8/Miniforge3 (Mamba) stack.
3. CPU Optimization: Patched src/egnn.py to force tensor operations onto the CPU, overcoming CUDA-only hardcoding.
4. Dependency Management: Downgraded pip < 24.0 to handle legacy pytorch-lightning metadata requirements.


**The Generation Pipeline:**

1. Input: 3D fragments of the original drug.
2. Model: geom_difflinker.ckpt (Pre-trained on GEOM dataset).
3. Process: Utilized 3D Equivariant Graph Neural Networks (EGNN) to predict atom types and positions in the 3D pocket.
4. Result: 50 unique .xyz candidates generated in ~38 minutes.


### Phase 3: Evaluation & Docking Pipeline

We implemented a "Fail Fast" hierarchical funnel to ensure only the most viable candidates were selected.

**1. Chemical Perception & Validity**
    We converted raw point clouds (.xyz) into chemical graphs (.sdf) using OpenBabel.
    
    Sanitization: 100% validity rate across 50 molecules.
    Drug-Likeness: Average QED score of ~0.70, indicating high drug-likeness.
    
**2. Molecular Docking (AutoDock Vina)**
Candidates were docked back into the HSP90 binding pocket to measure binding affinity ($\Delta G$).

    Receptor Prep: Created a PDBQT file with Gasteiger charges.
    Grid Box: Centered at (-1.0, 35.3, -0.8) with a 20Å dimension.
    Top Performance: The best candidates reached -9.9 kcal/mol, matching the potency of the original inhibitor.

### Evaluation Hierarchy

| Tier | Metric | Importance | Action |
| :--- | :--- | :--- | :--- |
| **1** | **RDKit Sanitization** |  **Critical** | Instant rejection of invalid valency (e.g., pentavalent carbon). |
| **2** | **Synthetic Accessibility (SA)** |  **High** | Filter out "hallucinations" (Target: **SA < 6.0**). |
| **3** | **Binding Affinity** |  **High** | Primary ranking factor based on calculated $\Delta G$ (kcal/mol). |
| **4** | **Ligand Efficiency (LE)** |  **Medium** | Tie-breaker favoring smaller, more potent binders. |



### Project Structure

```
.
├── data/                # Raw PDB and fragment files
├── checkpoints/         # DiffLinker model weights
├── scripts/             # Evaluation and patching scripts
├── docking/             # Vina config and receptor PDBQT
├── generated_molecules/ # Raw AI outputs (.xyz, .sdf)
└── evaluation/          # RDKit analysis and docking logs

```



## Installation and Environment Setup

### Environment Initialization
We recommend using Miniforge3 (Mamba) for high-performance dependency solving in WSL2 (Windows Subsystem for Linux 2).

```
# Create the environment with Python 3.8 (required for DiffLinker compatibility)

mamba create -n linker_design python=3.8 -y
mamba activate linker_design
```

### Core Dependencies (CPU-Only Configuration)
Because this project was optimized for CPU hardware without CUDA drivers, PyTorch must be installed via the cpuonly channel to prevent runtime errors.

```
# Install PyTorch CPU and legacy pip to handle older library metadata

mamba install pytorch cpuonly -c pytorch -y
pip install "pip<24.0" # Required for legacy Pytorch Lightning compatibility
```

### Geometric & Chemical Libraries
Geometric libraries must match the specific PyTorch build. Use the following binaries:
```
# Install Torch Geometric dependencies

pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
pip install torch-geometric pytorch-lightning==1.6.0

# Install Chemistry tools

mamba install -c conda-forge rdkit openbabel vina -y
pip install imageio scikit-learn
```

### Critical CPU Patches
The default DiffLinker scripts are hardcoded for CUDA. If you are running on a CPU, you must apply the following patches:

1. EGNN Tensor Routing: Open src/egnn.py and ensure all tensors are explicitly sent to the CPU.
2. Change: Locate line 462 and add .to('cpu') to edge tensors.
3. Change: Hardcode self.device = torch.device('cpu') in the EGNN class initialization.
4. Model Loading: In generate.py, ensure the checkpoint is loaded using map_location='cpu'.


## Results

### Tier 1: Gatekeepers (PAINS)

All 50 candidates passed the PAINS filter (no PAINS alerts detected) and were chemically valid.

### Tier 2: Reality Check (SA Score $\le 6.0$)

```Result: This was the most significant filter. Out of 50 candidates, only 5 molecules passed.```

The highest-performing molecules by docking score (e.g., output_31 at -9.9 kcal/mol) were discarded here because their SA Scores were $> 7.3$, indicating they are likely "model hallucinations" that are impossible to synthesize.


### Tier 3: Performance & Quality Ranking

The remaining 5 candidates were ranked by Binding Affinity ($\Delta G$).

| Filename | Binding Affinity | SA Score | Ligand Efficiency | Lipinski Violations |
| :--- | :---: | :---: | :---: | :---: |
| output_0_hsp90_fragments_ | -8.7 | 5.62 | 0.378 | 0 |
| output_38_hsp90_fragments_ | -8.7 | 5.62 | 0.378 | 0 |
| output_20_hsp90_fragments_ | -8.5 | 5.62 | 0.370 | 0 |
| output_37_hsp90_fragments_ | -8.0 | 5.62 | 0.348 | 0 |
| output_47_hsp90_fragments_ | -7.4 | 5.83 | 0.308 | 0 |