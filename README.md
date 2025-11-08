# 🧬 msprime and simuPOP Experiments

Just me experimenting with the [**msprime**](https://tskit.dev/msprime/docs/stable/intro.html) and [**simuPOP**](https://simupop.readthedocs.io/en/latest/) libraries for my senior capstone.  
The goal is to generate and analyze simulated genetic data under different demographic and evolutionary scenarios.

---

## 🧠 Project Overview

This project uses **`msprime`** to simulate realistic genetic variation across a cohort of individuals, then derives downstream phenotypes and covariates for statistical modeling.  
All generated files are stored in the **`datasets/`** directory and can be used directly for regression, classification, or GWAS-style analyses.

---

### 📂 Output Files

| File | Description |
|------|--------------|
| **`datasets/sim_cohort.csv`** | The main dataset containing one row per simulated individual. Includes demographic covariates (`sex`, `age`, `env_index`), genetic predictors (`polygenic_score`, optional `PC1…PCk`), and response variables (`quant_trait`, `disease_status`, `disease_prob`). Ideal for regression or classification experiments. |
| **`datasets/effect_sizes.csv`** | Metadata for all retained variants after MAF filtering. Contains each variant’s genomic position (`position`), minor allele frequency (`maf`), assigned effect size (`beta`), and a binary flag (`is_causal`) indicating whether it contributes to the simulated trait. Enables reproducibility and variant-level analysis. |
| **`datasets/sim_cohort.trees`** *(optional)* | The full `msprime` tree sequence containing the complete genealogical and mutational history of the simulation. Can be reloaded with `tskit.load()` for downstream population-genetics analyses or to regenerate genotype matrices. |

---

### ⚙️ Example Command

```bash
python make_msprime_dataset.py \
  --name A \
  --n 5000 \
  --seed 42 \
  --out sim_cohort.csv \
  --save-effects effect_sizes.csv \
  --save-ts sim_cohort.trees
```

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
python make_msprime_dataset.py --name A --n 5000 --seed 42 --out sim_cohort.csv
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
# venv\Scripts\activate     # Windows
```

### 3️⃣ Install dependencies
```bash
pip install -r requirements.txt
```

### 4️⃣ Run the simulation
```bash
python make_msprime_dataset.py --name A --n 5000 --seed 42 --out sim_cohort.csv
```