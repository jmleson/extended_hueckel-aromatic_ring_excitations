# Extended Hückel Model for Excitations in Monocyclic Aromatics

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-blue.svg)](https://creativecommons.org/licenses/by-nc/4.0)

> **Author**: [Judith M. Leson](https://orcid.org/0009-0002-5627-2673)  
> **Repository**: [https://github.com/jmleson/extended_hückel-aromatic_ring_excitations](https://github.com/jmleson/extended_hückel-aromatic_ring_excitations)  
> **Related Work**: Dissertation: *A Quantum-Chemical Analysis of Long-Range Dimer Interactions Arising From Triplet Excited States of Monocyclic Aromatics* (University of Duisburg-Essen, 2026)  
> **Data DOI**: [10.71955/DUEDATA-2026-MR4YY63J](https://doi.org/10.71955/DUEDATA-2026-MR4YY63J)  
> 
> **References**:
> - Hapka et al., *J. Chem. Phys.* **2019**, *150*, 164107. 
> - Jiemchooroj et al., *J. Chem. Phys.* **2005**, *122*, 164104.
---

## 📌 Overview

This repository contains a model based on the **Hückel Theory** to derive analytical expressions for dispersion coefficients in molecules with planar **π-systems**. 
The focus is on doubly excited states in the dimer, allowing for varying substitution patterns.

While open-chain systems are supported, the model was designed for **monocyclic aromatic systems**, such as:
- Cyclopentadiene
- Benzene
- Pyridine
- Pyrazine
- Chlorobenzene
- Dichlorobenzenes
- ...

> ✅ The output is a **traceable, and human-readable** LaTeX (or PDF) file for the given molecule.
---

## ✅ Key Features

- **Symbolic derivation** of dispersion coefficients for dimer interactions
- **Inclusion of different excited states** namely ground state interactions, and quintet states in dimers
- **Step-by-step derivation** of Hückel-based molecular orbitals, Hamiltonian matrices, and secular equations  
- **Support for diverse molecular geometries**, including cyclic and open-chain systems, with or without substituents bearing p-orbitals 
- **Visualization of molecular geometries** in png files
- **Automated LaTeX generation** for model results  
- **Human-readable derivations in compiled PDFs**, showing intermediate steps


---
## 🛠️ Basic Concept 
We consider the p-orbitals that are orthogonal to the π-system of a given molecule.  
Excitations are modeled as electron transitions from a p-orbital to a higher, unoccupied s-orbital.

Hückel theory provides access to the molecular orbitals and their energies.  
We label orbitals derived from p-orbitals as $\psi$ and those derived from s-orbitals as $\eta$.

The code computes the dispersion energy for different combinations of monomer states using the formula (Hapka 2019, Jiemchooroj 2005): 

$$
E_{\text{dispersion}}({0}, {\tilde{0}}) = \frac{1}{R^6} \cdot C_6
$$

$$
= \frac{1}{R^6} \cdot
        {\sum_{i,a}} {\sum_{j,b} }
        \frac{
          \left(\quad
            {⟨\psi_0| \hat{\mu}_z |\eta_i^a}⟩ \quad T_{zz} \quad
            {⟨\tilde{\psi_0}| \hat{\mu}_z |\tilde{\eta}_j^b⟩}
          \quad\right)^2
        }{
          {\left(E_i^a  - E_0 \right)} + {\left( \tilde{E}_j^b - \tilde{E}_0 \right)}
        }
$$

This expression accounts for the dispersion interaction between two monomers ($0$, or $\tilde{0}$) and excited states (transition from $i$ into $a$, or from $j$ into $b$).  
The sum runs over all excited states accessible from the initial states of the monomers.

Our code is able to compute the resulting dispersion coefficients $C_6$ for various substituted aromatic systems.



___

## 🔧 How to Use
### 📂 Structure
- `requirements.txt`: Needed python libraries
- `run.py`: Main script that constructs the model for exemplary molecular systems
- `src/`: Folder for source files (Python)
- `RESULTS`: Folder for resulting LaTeX and png files
  - `RESULTS/run_latex.sh`: Compiles `.tex` files into PDFs (output in `RESULTS/out/`)
  - `RESULTS/out/`: Folder for LaTeX compilation files
- `testing`: Folder for unittests, including manually checked derivation parts 


### 📦 Dependencies

This project requires libraries such as:
- `pypointgroup` 
- `sympy` 
- `fractions`
To install the needed requirements, use:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


### 💻 Running the Code
Exemplary usage is provided in `run.py`. 
For instance, to generate a LaTeX file for benzene:
```python
from src.main.MoleculeState import MoleculeState

m = MoleculeState(n=6, bound_cl_to_c_positions=[],n_instead_of_c=[])
m.latex_datei_erstellen("Benzene")
```
This creates a traceable derivation of dispersion coefficients in `RESULTS/Benzene.tex`, that can be compiled via `./RESULTS/run_latex.sh` which yields the corresponding file `RESULTS/out/Benzene.pdf`.  

> ❗ The code prioritizes traceability over speed. Solving the secular equations may take considerable time or even fail, depending on the system. 





### 📊 Print Results to Console
To print results directly to the console instead of generating a file, use:
```python
from src.main.MoleculeState import MoleculeState

m = MoleculeState(n=6, bound_cl_to_c_positions=[],n_instead_of_c=[])
m.calculate_result_for_all_transitions(print_active=print_active)
m.calculate_result_for_all_transitions_results(print_active=print_active)
m.compare_ground_state_assumption(print_active=print_active)
``` 


### 🧪 Run Unit Tests
To run the unit tests, execute the following command from the project root:
```bash
python3 -m testing.TestMoleculeState
python3 -m testing.TestTransitionIntegral
python3 -m testing.TestLinearButadiene
python3 -m testing.TestCyclobutadiene
python3 -m testing.TestBenzene
python3 -m testing.TestChlorobenzene
python3 -m testing.TestPyrazine
python3 -m testing.TestPyridine
```
This runs all test cases and verifies the correctness of the model.
