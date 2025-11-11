# 🧬 msprime and simuPOP Experiments

Just me experimenting with the [**msprime**](https://tskit.dev/msprime/docs/stable/intro.html) and [**simuPOP**](https://simupop.readthedocs.io/en/latest/) libraries for my senior capstone.  
The goal is to generate and analyze simulated genetic data under different demographic and evolutionary scenarios.

---

## 🧠 Project Overview

This project uses **`msprime`** to simulate realistic genetic variation across a cohort of individuals, then derives downstream phenotypes and covariates for statistical modeling.  
All generated files are stored in the **`datasets/`** directory and can be used directly for regression, classification, or GWAS-style analyses.

---

## msprime

### 📂 Output Files (msprime)

| File                                         | Description |
|----------------------------------------------|--------------|
| **`datasets/msprime_sim_cohort.csv`**        | The main dataset containing one row per simulated individual. Includes demographic covariates (`sex`, `age`, `env_index`), genetic predictors (`polygenic_score`, optional `PC1…PCk`), and response variables (`quant_trait`, `disease_status`, `disease_prob`). Ideal for regression or classification experiments. |
| **`datasets/msprime_effect_sizes.csv`**      | Metadata for all retained variants after MAF filtering. Contains each variant’s genomic position (`position`), minor allele frequency (`maf`), assigned effect size (`beta`), and a binary flag (`is_causal`) indicating whether it contributes to the simulated trait. Enables reproducibility and variant-level analysis. |
| **`datasets/msprime_sim_cohort.trees`** *(optional)* | The full `msprime` tree sequence containing the complete genealogical and mutational history of the simulation. Can be reloaded with `tskit.load()` for downstream population-genetics analyses or to regenerate genotype matrices. |

---

### ⚙️ How `make_msprime_dataset.py` Works

- **Coalescent Simulation:** Uses `msprime` to generate a tree sequence representing ancestral lineages and mutations under a Wright–Fisher coalescent model with configurable population size (`--n`), mutation rate (`--mu`), recombination rate (`--r`), and sequence length (`--seq-len`).  
- **Genotype Matrix:** Converts the tree sequence into a genotype matrix (0/1/2 dosages per locus) using the built-in `tskit` genotype API.  
- **MAF Filtering:** Filters loci to retain variants within a specified minor allele frequency range (`--maf-min`, `--maf-max`).  
- **Polygenic Scores:** Assigns random effect sizes to a proportion of loci (`--prop-causal`) to generate a standardized polygenic score (`polygenic_score`).  
- **Covariates:** Adds simulated demographic and environmental variables (`sex`, `age`, `env_index`).  
- **Traits:** Derives both continuous (`quant_trait`) and binary (`disease_status`) outcomes using linear and logistic models based on the polygenic score and covariates.  
- **PCA (optional):** Performs principal component analysis of the genotype matrix to create population-structure covariates (`--pcs` components).  
- **Outputs:** Saves a per-individual CSV dataset and an effect-size CSV file, plus an optional `.trees` file for full genealogical reconstruction.

---

### ⚙️ Example Command

```bash
python make_msprime_dataset.py \
  --name A \
  --n 5000 \
  --seed 42 \
  --out sim \
  --save-effects effect \
  --save-ts tree
```

This produces:
- `datasets/A_msprime_sim.csv`
- `datasets/A_msprime_effect.csv`
- `datasets/A_msprime_tree.csv`

---

## simuPOP

### 📂 Output Files (simuPOP)

| File                                         | Description |
|----------------------------------------------|--------------|
| **`datasets/simupop_sim_cohort.csv`**        | The main dataset simulated via **simuPOP**. Each row represents a diploid individual with demographic covariates (`sex`, `age`, `env_index`), genetic predictors (`polygenic_score`, optional `PC1…PCk`), and response variables (`quant_trait`, `disease_status`, `disease_prob`). Mirrors the msprime dataset format for consistency across modeling experiments. |
| **`datasets/simupop_effect_sizes.csv`**      | Variant-level metadata produced alongside the cohort. Contains the locus position along a synthetic chromosome (`position`), filtered minor allele frequency (`maf`), simulated effect size (`beta`), and causal flag (`is_causal`). Useful for evaluating how polygenic signal and allele frequencies relate to phenotype generation. |

---

### ⚙️ How `make_simupop_dataset.py` Works

- **Forward Simulation:** Uses `simuPOP` to evolve a structured diploid population with configurable mutation (`--mu`), recombination (`--r`), migration (`--mig`), and subpopulation (`--subpops`) parameters.  
- **MAF Filtering:** Filters loci outside the desired minor-allele-frequency range (`--maf-min`, `--maf-max`).  
- **Polygenic Scores:** Randomly designates a subset of SNPs as causal (`--prop-causal`) and generates a standardized polygenic score (`polygenic_score`).  
- **Covariates:** Adds synthetic environmental and demographic covariates (`sex`, `age`, `env_index`).  
- **Traits:** Produces both continuous (`quant_trait`) and binary (`disease_status`, with `disease_prob`) phenotypes using logistic and linear models incorporating the polygenic score and covariates.  
- **PCA (optional):** Performs lightweight genotype PCA using NumPy SVD (`--pcs` controls component count).  
- **Outputs:** Two CSVs written to `datasets/`: one for the cohort, one for variant metadata.

---

### ⚙️ Example Command

```bash
python make_simupop_dataset.py \
  --name B \
  --n 5000 \
  --seed 42 \
  --out sim \
  --save-effects effects
```

This produces:
- `datasets/B_simupop_sim.csv`
- `datasets/B_simupop_effects.csv`

---

## 📂 Public Datasets

For the data found in the `public_datasets` folder here is how the data was generated:

msprime
```bash
python make_msprime_dataset.py --name public --n 5000 --seed 77 --out sim --save-effects effects --save-ts tree
```

simupop
```bash
python make_simupop_dataset.py --name public --n 5000 --seed 42 --out sim --save-effects effects
```

At the moment, simuPop still generates data with a random seed, I will fix this later. For the run the seed was:
`Random Number Generator is set to mt19937 with random seed 0x54f111254a5d51ee.`

---

## ⚙️ Installation (Conda)

### 1️⃣ Clone the repository
```bash
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>
```

### 2️⃣ Create and activate the environment
```bash
conda env create -f environment.yml
conda activate capstone-msprime
```

### 3️⃣ Run the simulation
```bash
python make_msprime_dataset.py --name public --trees yes --n 5000

python make_simupop_dataset.py --name public --n 5000
```

---

## 🐍 Installation (pip)

### 1️⃣ Clone the repository
```bash
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>
```

### 2️⃣ Create and activate a virtual environment
```bash
python3 -m venv .venv
source .venv/bin/activate  # macOS/Linux
# venv\\Scripts\\activate     # Windows
```

### 3️⃣ Install dependencies
```bash
pip install -r requirements.txt
```

### 4️⃣ Run the simulation
```bash
python make_msprime_dataset.py --name public --n 5000

python make_simupop_dataset.py --name public --n 5000
```

