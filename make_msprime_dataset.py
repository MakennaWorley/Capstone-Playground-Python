import os, re, argparse

import numpy as np
import pandas as pd
import msprime
import tskit
from sklearn.decomposition import PCA

DEFAULT_OUT = "sim_cohort"
DEFAULT_EFF = "effect_sizes"
DEFAULT_TREE = "sim_cohort"

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1", "y"):
        return True
    elif v.lower() in ("no", "false", "f", "0", "n"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def slugify(s: str) -> str:
    # keep letters, digits, dot, dash, underscore; collapse others to "_"
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")
    return s or "run"

def simulate_tree_sequence(n_samples: int,
                           seq_len: float = 1e6,
                           mu: float = 1e-8,
                           r: float = 1e-8,
                           Ne: int = 10_000,
                           seed: int = 42) -> tskit.TreeSequence:
    """
    Simulate a simple neutral diploid cohort using msprime.
    - n_samples: number of diploid individuals
    """
    ts = msprime.sim_ancestry(
        samples=n_samples,
        recombination_rate=r,
        sequence_length=seq_len,
        population_size=Ne,
        ploidy=2,
        model="dtwf",
        random_seed=seed
    )
    ts = msprime.sim_mutations(
        ts,
        rate=mu,
        model=msprime.SLiMMutationModel(type=1),  # neutral, infinite sites-like
        random_seed=seed + 1
    )
    return ts

def ts_to_genotype_matrix(ts: tskit.TreeSequence, n_diploid: int):
    """
    Build (variants x samples) genotype matrix with genotypes in {0,1,2}.
    """
    # tskit.genotype_matrix returns haploid counts; for diploids we reshape
    G_hap = ts.genotype_matrix().T  # shape: (num_samples_hap, num_variants)
    # Combine every two haploid rows into one diploid row by summing
    if G_hap.shape[0] != n_diploid * 2:
        raise ValueError("Haploid count does not match requested diploids.")
    G_dip = G_hap.reshape(n_diploid, 2, G_hap.shape[1]).sum(axis=1)  # (n_diploid, num_variants)
    return G_dip  # each entry is 0,1,2

def maf_filter(G: np.ndarray, maf_min: float = 0.01, maf_max: float = 0.5):
    """
    Filter variants by minor allele frequency.
    """
    n = G.shape[0]
    p = G.sum(axis=0) / (2.0 * n)          # allele frequency
    maf = np.minimum(p, 1 - p)
    keep = (maf >= maf_min) & (maf <= maf_max)
    return G[:, keep], keep, maf[keep]

def build_polygenic_score(G: np.ndarray, prop_causal: float = 0.05, rng=None):
    """
    Randomly choose a subset of variants as causal and draw effect sizes.
    Returns polygenic score vector and effect-size vector aligned to columns of G.
    """
    if rng is None:
        rng = np.random.default_rng()
    m = G.shape[1]
    m_causal = max(1, int(np.round(prop_causal * m)))
    causal_idx = rng.choice(m, size=m_causal, replace=False)
    betas = np.zeros(m)
    # Heuristic: smaller effects for polygenicity
    betas[causal_idx] = rng.normal(loc=0.0, scale=0.05, size=m_causal)
    prs = G @ betas
    # Z-score PRS for stability
    prs = (prs - prs.mean()) / (prs.std(ddof=1) + 1e-8)
    return prs, betas, causal_idx

def make_covariates(n: int, rng=None):
    """
    sex (binary), age (numeric), env_index (numeric). All realistic ranges.
    """
    if rng is None:
        rng = np.random.default_rng()
    sex = rng.integers(0, 2, size=n)                  # 0/1
    age = rng.integers(20, 71, size=n).astype(float)  # 20–70
    env_index = rng.normal(0, 1, size=n)              # standardized environmental exposure
    return sex, age, env_index

def make_quant_trait(prs, sex, age, env, rng=None):
    """
    Build a quantitative trait from PRS + covariates.
    """
    if rng is None:
        rng = np.random.default_rng()
    # Coefficients (tuneable)
    beta0 = 0.0
    b_prs = 1.0
    b_sex = 0.3
    b_age = 0.02
    b_env = 0.5
    eps = rng.normal(0, 1, size=prs.shape[0])

    y = beta0 + b_prs * prs + b_sex * sex + b_age * (age - age.mean())/age.std(ddof=1) + b_env * env + eps
    # Standardize for convenience
    y = (y - y.mean()) / (y.std(ddof=1) + 1e-8)
    return y

def make_binary_outcome(prs, sex, age, env, rng=None):
    """
    Logistic disease risk from PRS + covariates; returns 0/1.
    """
    if rng is None:
        rng = np.random.default_rng()
    # Logistic coefficients (tuneable)
    alpha = -0.2                # intercept => prevalence around ~0.4–0.6 with chosen betas
    b_prs = 0.9
    b_sex = 0.4
    b_age = 0.03
    b_env = 0.6

    age_std = (age - age.mean())/age.std(ddof=1)
    lin = alpha + b_prs*prs + b_sex*sex + b_age*age_std + b_env*env
    p = 1 / (1 + np.exp(-lin))
    y_bin = (rng.uniform(0, 1, size=p.shape[0]) < p).astype(int)
    return y_bin, p

def compute_pcs(G: np.ndarray, n_components=2):
    """
    Quick PCA on standardized genotype matrix (centered per SNP).
    """
    # Center columns
    Gc = G - G.mean(axis=0, keepdims=True)
    # Scale to unit variance per SNP (avoid div by 0)
    sd = Gc.std(axis=0, ddof=1, keepdims=True) + 1e-8
    Gz = Gc / sd
    pca = PCA(n_components=n_components, svd_solver=" randomized ".strip())
    PCs = pca.fit_transform(Gz)
    return PCs, pca.explained_variance_ratio_

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", type=str, default=None,
                    help="Base name for outputs: produces <name>_msprime_sim_cohort.csv, <name>_msprime_effect_sizes.csv, (optional) <name>_msprime_sim_cohort.trees")
    ap.add_argument("--trees", type=str2bool, nargs="?", const=True, default=False,
                    help="If set, also save the tree sequence. (With --name, becomes <name>_msprime_sim_cohort.trees)")

    ap.add_argument("--n", type=int, default=5000, help="Number of diploid individuals")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (if omitted, a random seed will be chosen and printed)")

    # Keep explicit paths for backwards compatibility; they’ll override auto names if provided.
    ap.add_argument("--out", type=str, default=DEFAULT_OUT, help="Name for output simulation data file")
    ap.add_argument("--save-effects", type=str, default=DEFAULT_EFF, help="Name for output variant effects file")
    ap.add_argument("--save-ts", type=str, default=None, help="(Optional) Name for output .trees file")

    ap.add_argument("--pcs", type=int, default=2, help="Number of PCs to include (0 to skip)")
    ap.add_argument("--maf-min", type=float, default=0.01, help="Minor allele frequency lower bound")
    ap.add_argument("--maf-max", type=float, default=0.5, help="Minor allele frequency upper bound")
    ap.add_argument("--prop-causal", type=float, default=0.05, help="Proportion of variants set as causal")

    args = ap.parse_args()

    # seed logic
    if args.seed is None:
        args.seed = int(np.random.SeedSequence().entropy % (2**32))
        print(f"[info] No --seed provided. Using randomly generated seed: {args.seed}")
    else:
        print(f"[info] Using user-specified seed: {args.seed}")
    rng = np.random.default_rng(args.seed)
    np.random.seed(args.seed)

    # naming logic
    base = slugify(args.name) if args.name else None

    # If user didn't override the defaults, auto-rename from --name
    if base:
        prefix = f"{base}_msprime"
    else:
        prefix = "msprime"

    if args.out == DEFAULT_OUT:
        args.out = f"{prefix}_{DEFAULT_OUT}.csv"
    else:
        args.out = f"{prefix}_{args.out}.csv"

    if args.save_effects == DEFAULT_EFF:
        args.save_effects = f"{prefix}_{DEFAULT_EFF}.csv"
    else:
        args.save_effects = f"{prefix}_{args.save_effects}.csv"

    if args.save_ts is not None:
        args.save_ts = f"{prefix}_{args.save_ts}.trees"
    elif args.trees:
        args.save_ts = f"{prefix}_{DEFAULT_TREE}.trees"

    # output dir
    output_dir = os.path.join(os.getcwd(), "datasets")
    os.makedirs(output_dir, exist_ok=True)

    # Redirect output paths
    args.out = os.path.join(output_dir, args.out)
    args.save_effects = os.path.join(output_dir, args.save_effects)
    if args.save_ts:
        args.save_ts = os.path.join(output_dir, args.save_ts)

    # 1) Simulate ancestry & mutations
    ts = simulate_tree_sequence(n_samples=args.n, seed=args.seed)
    if args.save_ts:
        ts.dump(args.save_ts)

    # 2) Genotype matrix and MAF filter
    G = ts_to_genotype_matrix(ts, n_diploid=args.n)             # (N, M_all)
    Gf, keep_mask, kept_maf = maf_filter(G, maf_min=args.maf_min, maf_max=args.maf_max)

    # 3) Polygenic score
    prs, betas, causal_idx = build_polygenic_score(Gf, prop_causal=args.prop_causal, rng=rng)

    # 4) Covariates
    sex, age, env = make_covariates(args.n, rng=rng)

    # 5) Traits
    quant_trait = make_quant_trait(prs, sex, age, env, rng=rng)
    disease_status, disease_prob = make_binary_outcome(prs, sex, age, env, rng=rng)

    # 6) Optional PCs
    df = pd.DataFrame({
        "individual_id": np.arange(args.n, dtype=int),
        "sex": sex.astype(int),
        "age": age.astype(float),
        "env_index": env.astype(float),
        "polygenic_score": prs.astype(float),
        "quant_trait": quant_trait.astype(float),
        "disease_status": disease_status.astype(int),
        "disease_prob": disease_prob.astype(float),
    })

    if args.pcs > 0:
        PCs, var_ratio = compute_pcs(Gf, n_components=args.pcs)
        for k in range(args.pcs):
            df[f"PC{k+1}"] = PCs[:, k].astype(float)

    # 7) Save dataset
    df.to_csv(args.out, index=False)

    # 8) Save effects/variant metadata for reproducibility
    kept_sites = [v.site.position for v in ts.variants() if keep_mask[v.site.id]]
    kept_sites = np.array(kept_sites)
    causal_flags = np.zeros(Gf.shape[1], dtype=int)
    causal_flags[causal_idx] = 1
    eff = pd.DataFrame({
        "variant_index": np.arange(Gf.shape[1], dtype=int),
        "position": kept_sites,
        "maf": kept_maf,
        "beta": betas,
        "is_causal": causal_flags
    })
    eff.to_csv(args.save_effects, index=False)

    # 9) Minimal data dictionary (print to console)
    print(f"\nSaved simulation to: {args.out}")
    print(f"Saved variant effects to: {args.save_effects}")
    if args.save_ts:
        print(f"Saved tree sequence to: {args.save_ts}")

if __name__ == "__main__":
    main()
