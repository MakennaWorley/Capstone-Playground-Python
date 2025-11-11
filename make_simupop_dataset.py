import os, re, argparse, sys, math
import numpy as np
import pandas as pd
import simuPOP as sim

DEFAULT_OUT = "sim_cohort"
DEFAULT_EFF = "effect_sizes"

def slugify(s: str) -> str:
    # keep letters, digits, dot, dash, underscore; collapse others to "_"
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")
    return s or "run"

# ---------- Optional PCA without scikit-learn (numpy SVD) ----------
def pca_np(X: np.ndarray, n_components=2):
    """
    Compute PCA using numpy SVD. X should be standardized (mean=0 per column).
    Returns PCs (n_samples x n_components) and explained variance ratio.
    """
    # Center columns
    Xc = X - X.mean(axis=0, keepdims=True)
    # Scale to unit variance to avoid dominance by high-MAF loci
    sd = Xc.std(axis=0, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    Z = Xc / sd
    # SVD on (n x m) standardized matrix
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    # Principal components are columns of U * S
    PCs = U[:, :n_components] * S[:n_components]
    # Explained variance per component
    # Total variance equals sum of singular values squared divided by (n-1)
    var_total = (S**2).sum() / (Z.shape[0] - 1)
    var_comp = (S[:n_components]**2) / (Z.shape[0] - 1)
    evr = var_comp / var_total if var_total > 0 else np.zeros_like(var_comp)
    return PCs, evr

# ---------- Traits ----------
def build_polygenic_score(G: np.ndarray, prop_causal: float = 0.05, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    m = G.shape[1]
    if m == 0:
        # Return zeros so downstream code still runs; caller should consider increasing loci or relaxing MAF
        return np.zeros(G.shape[0]), np.zeros(0), np.array([], dtype=int)
    m_causal = max(1, int(round(prop_causal * m)))
    m_causal = min(m_causal, m)
    causal_idx = rng.choice(m, size=m_causal, replace=False)
    betas = np.zeros(m, dtype=float)
    betas[causal_idx] = rng.normal(loc=0.0, scale=0.05, size=m_causal)
    prs = G @ betas
    prs = (prs - prs.mean()) / (prs.std(ddof=1) + 1e-8)
    return prs, betas, causal_idx

def make_covariates(n: int, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    sex = rng.integers(0, 2, size=n)                  # 0/1 (binary categorical)
    age = rng.integers(20, 71, size=n).astype(float)  # 20–70
    env_index = rng.normal(0, 1, size=n)              # standardized env exposure
    return sex, age, env_index

def make_quant_trait(prs, sex, age, env, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    beta0 = 0.0
    b_prs = 1.0
    b_sex = 0.3
    b_age = 0.02
    b_env = 0.5
    age_std = (age - age.mean())/(age.std(ddof=1) + 1e-8)
    eps = rng.normal(0, 1, size=prs.shape[0])
    y = beta0 + b_prs*prs + b_sex*sex + b_age*age_std + b_env*env + eps
    y = (y - y.mean())/(y.std(ddof=1) + 1e-8)
    return y

def make_binary_outcome(prs, sex, age, env, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    alpha = -0.2
    b_prs = 0.9
    b_sex = 0.4
    b_age = 0.03
    b_env = 0.6
    age_std = (age - age.mean())/(age.std(ddof=1) + 1e-8)
    lin = alpha + b_prs*prs + b_sex*sex + b_age*age_std + b_env*env
    p = 1.0 / (1.0 + np.exp(-lin))
    y_bin = (rng.uniform(0, 1, size=p.shape[0]) < p).astype(int)
    return y_bin, p

# ---------- simuPOP pipeline ----------
def simulate_with_simupop(n: int,
                          loci: int = 20000,
                          seq_len: float = 1e6,
                          mu: float = 1e-8,
                          r: float = 1e-8,
                          generations: int = 50,
                          subpops: int = 2,
                          mig: float = 0.01,
                          seed: int = 42):
    """
    Simulate a diploid cohort using simuPOP with SNPs and simple structure.
    Returns genotype matrix G (n x loci) in {0,1,2} and locus positions along a
    1D chromosome of length seq_len.
    """
    rng = np.random.default_rng(seed)
    sim.setOptions(seed=seed)

    # Split individuals across subpops
    sizes = [n // subpops] * subpops
    sizes[-1] += n - sum(sizes)

    # Randomly distribute loci along [0, seq_len]
    lociPos = np.sort(rng.uniform(0, seq_len, size=loci)).tolist()

    # Recombination intensity in Morgans ~= r * seq_len per meiosis (approx.)
    # simuPOP's Recombinator(intensity=...) uses Morgans per genome.
    rec_intensity = r * seq_len

    pop = sim.Population(size=sizes, loci=[loci], ploidy=2,
                         lociPos=lociPos,
                         infoFields=['ind_id', 'migrate_to'])

    # Tag with individual IDs
    sim.IdTagger().apply(pop)

    # Init sexes randomly across subpops
    pop.evolve(
        initOps=[
            sim.InitSex(),
            # Seed initial polymorphism so SNPMutator has variation to work with
            sim.InitGenotype(freq=[0.95, 0.05])
        ],
        # Build preOps without None values
        preOps=(
            [sim.SNPMutator(u=mu)] +
            ([sim.Migrator(rate=[[0 if i==j else (mig if subpops==2 else mig/(subpops-1)) for j in range(subpops)] for i in range(subpops)])] if (subpops > 1 and mig > 0) else [])
        ),
        # Random mating with recombination
        matingScheme=sim.RandomMating(ops=[sim.Recombinator(intensity=rec_intensity)]),
        postOps=[],
        gen=generations
    )

    # Collect genotypes as 0/1/2 dosages per locus
    # For each individual, genotype() returns a flat array of length ploidy*loci
    # e.g., [a1_l1, a1_l2, ..., a1_lL, a2_l1, a2_l2, ..., a2_lL]
    G = np.empty((n, loci), dtype=np.int8)
    idx = 0
    for sp in range(pop.numSubPop()):
        for i in range(pop.subPopSize(sp)):
            ind = pop.individual(idx)
            g = np.array(ind.genotype(), dtype=np.int8)
            g = g.reshape((2, loci)).sum(axis=0)  # dosage 0,1,2 per locus
            G[idx, :] = g
            idx += 1

    # Minor allele frequency spectrum
    p = G.sum(axis=0) / (2.0 * n)
    maf = np.minimum(p, 1 - p)
    return G, np.array(lociPos), maf

def maf_filter(G: np.ndarray, maf_min=0.01, maf_max=0.5):
    n = G.shape[0]
    if G.size == 0:
        return G, np.array([], dtype=bool), np.array([])
    p = G.sum(axis=0) / (2.0 * n)
    maf = np.minimum(p, 1 - p)
    keep = (maf >= maf_min) & (maf <= maf_max)
    # Adaptive relaxation if nothing passes
    if keep.sum() == 0:
        # keep any polymorphic site
        poly = (maf > 0)
        if poly.sum() > 0:
            keep = poly
        else:
            # all monomorphic; return as-is and let caller handle
            keep = np.zeros(G.shape[1], dtype=bool)
    return G[:, keep], keep, maf[keep]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", type=str, default=None,
                    help="Base name for outputs: makes <name>_simupop_sim_cohort.csv, <name>_simupop_effect_sizes.csv")

    ap.add_argument("--n", type=int, default=5000, help="Number of diploid individuals")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (if omitted, a random seed will be chosen and printed)")

    ap.add_argument("--pcs", type=int, default=2, help="Number of genotype PCs to include (0 to skip)")
    ap.add_argument("--maf-min", type=float, default=0.01, help="Minor allele frequency lower bound")
    ap.add_argument("--maf-max", type=float, default=0.5, help="Minor allele frequency upper bound")
    ap.add_argument("--prop-causal", type=float, default=0.05, help="Proportion of variants considered causal for PRS")

    # Output paths (similar to msprime script)
    ap.add_argument("--out", type=str, default=DEFAULT_OUT, help="Name for output simulation data file")
    ap.add_argument("--save-effects", type=str, default=DEFAULT_EFF, help="Name for output variant effects file")

    ap.add_argument("--loci", type=int, default=20000, help="Number of SNP loci")
    ap.add_argument("--seq-len", type=float, default=1e6, help="Chromosome length (bp) to place loci for recombination scaling")
    ap.add_argument("--mu", type=float, default=1e-8, help="Per-locus per-generation mutation rate for SNPMutator")
    ap.add_argument("--r", type=float, default=1e-8, help="Per-base recombination rate; intensity approximated as r*seq_len")
    ap.add_argument("--gens", type=int, default=50, help="Forward generations to evolve")
    ap.add_argument("--subpops", type=int, default=2, help="Number of subpopulations to induce structure")
    ap.add_argument("--mig", type=float, default=0.01, help="Symmetric migration rate between subpops")

    args = ap.parse_args()

    # seed logic
    if args.seed is None:
        args.seed = int(np.random.SeedSequence().entropy % (2 ** 32))
        print(f"[info] No --seed provided. Using randomly generated seed: {args.seed}")
    else:
        print(f"[info] Using user-specified seed: {args.seed}")
    rng = np.random.default_rng(args.seed)
    np.random.seed(args.seed)

    # naming logic
    base = slugify(args.name) if args.name else None

    # If user didn't override the defaults, auto-rename from --name
    if base:
        prefix = f"{base}_simupop"
    else:
        prefix = "simupop"

    if args.out == DEFAULT_OUT:
        args.out = f"{prefix}_{DEFAULT_OUT}.csv"
    else:
        args.out = f"{prefix}_{args.out}.csv"

    if args.save_effects == DEFAULT_EFF:
        args.save_effects = f"{prefix}_{DEFAULT_EFF}.csv"
    else:
        args.save_effects = f"{prefix}_{args.save_effects}.csv"

    # output dir
    output_dir = os.path.join(os.getcwd(), "datasets")
    os.makedirs(output_dir, exist_ok=True)

    # Redirect output paths
    args.out = os.path.join(output_dir, args.out)
    args.save_effects = os.path.join(output_dir, args.save_effects)

    # 1) Simulate genotype matrix with simuPOP
    try:
        G_all, loci_pos, maf_all = simulate_with_simupop(
            n=args.n, loci=args.loci, seq_len=args.seq_len, mu=args.mu, r=args.r,
            generations=args.gens, subpops=args.subpops, mig=args.mig, seed=args.seed
        )
    except Exception as e:
        # Bubble up to user with clear message
        print(f"\nFATAL: Could not run simuPOP simulation: {e}", file=sys.stderr)
        sys.exit(1)

    # 2) MAF filter
    Gf, keep_mask, kept_maf = maf_filter(G_all, maf_min=args.maf_min, maf_max=args.maf_max)

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
        PCs, evr = pca_np(Gf, n_components=args.pcs)
        for k in range(args.pcs):
            df[f"PC{k+1}"] = PCs[:, k].astype(float)

    # 7) Save dataset
    df.to_csv(args.out, index=False)

    # 8) Save effects and locus metadata
    kept_positions = np.array(loci_pos)[keep_mask]
    causal_flags = np.zeros(Gf.shape[1], dtype=int)
    causal_flags[causal_idx] = 1
    eff = pd.DataFrame({
        "variant_index": np.arange(Gf.shape[1], dtype=int),
        "position": kept_positions.astype(float),
        "maf": kept_maf.astype(float),
        "beta": betas.astype(float),
        "is_causal": causal_flags.astype(int)
    })
    eff.to_csv(args.save_effects, index=False)

    # 9) Minimal data dictionary (stdout)
    print(f"\nSaved simulation to: {args.out}")
    print(f"Saved variant effects to: {args.save_effects}")

if __name__ == "__main__":
    main()
