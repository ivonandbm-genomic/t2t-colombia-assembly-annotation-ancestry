
____________________________________________________________________________________
#!/usr/bin/env python3

"""
Python script for processing multiple DeepVariant VCF files generated from PacBio HiFi sequencing:
- Filters only PASS biallelic SNPs
- Applies quality filters
- Normalizes variants using a reference genome
- Assigns unique IDs to each variant
- Indexes VCF files
"""

import os
import subprocess
from pathlib import Path
from glob import glob

# ──────────────────────────────
# USER CONFIGURATION
# ──────────────────────────────
VCF_DIR = Path("/mnt/diskrare/arlenb/08/small_variants/hg38")
BAM_DIR = Path("/mnt/diskrare/arlenb/08/aligned_reads/hg38")
OUT_DIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil")
REF = "/home/rare/ivon/data/hg38.fa")  # Path to the reference genome
THREADS = 32


# Create output directory if it does not exist
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Retrieve list of VCF files
vcfs = sorted(VCF_DIR.glob("*.vcf.gz"))
print(f"{len(vcfs)} VCF files detected for processing.")

# File processing loop
for vcf in vcfs:
    basename = vcf.stem.replace("_aligned", "")
    print(f"\nProcessing {basename}...")

    # Step 1: Filter PASS biallelic SNPs only
    step1 = OUT_DIR / f"{basename}.pass.snps.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-f", "PASS", "-v", "snps",
        "-Oz", "-o", str(step1), str(vcf)
    ], check=True)

    # Step 2: Quality filtering (QUAL >= 20, DP >= 12, GQ >= 20, no missing genotypes)
    step2 = OUT_DIR / f"{basename}.filtered.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-i", 'QUAL >= 20 & FORMAT/DP >= 12 & FORMAT/GQ >= 20 & GT!="mis"',
        "-Oz", "-o", str(step2), str(step1)
    ], check=True)

    # Step 3: Normalize variants using the reference genome
    step3 = OUT_DIR / f"{basename}.norm.vcf.gz"
    subprocess.run([
        "bcftools", "norm",
        "-f", REF,
        "-Oz", "-o", str(step3), str(step2)
    ], check=True)

    # Step 4: Assign a unique ID to each SNP
    step4 = OUT_DIR / f"{basename}.norm.id.vcf.gz"
    subprocess.run([
        "bcftools", "annotate",
        "--set-id", "+'%CHROM\_%POS\_%REF\_%ALT'",
        "-Oz", "-o", str(step4), str(step3)
    ], check=True)

    # Step 5: Index the final VCF file
    subprocess.run(["bcftools", "index", str(step4)], check=True)

print("\nProcessing completed successfully for all VCF files.")



______________________________________________________________________________REFERENCE PANEL 


#2_REFERENCE PANEL STEP 1_GENERATION OF panel_5superpopulations.tsv#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sys

BASE = "/mnt/diskrare/ivonb/refamerindios"

PANEL_TSV = os.path.join(BASE, "panel_5superpopulations.tsv")
PANEL_IDS = os.path.join(BASE, "panel_5super.ids")

VCF_BASE = BASE
OUT_BASE = os.path.join(BASE, "vcfs_panel")

BCFTOOLS = "bcftools"          # or absolute path if required
THREADS = 32                   # number of threads


def run_cmd(cmd, **kwargs):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)


def main():
    os.makedirs(OUT_BASE, exist_ok=True)

    if not os.path.isfile(PANEL_TSV):
        print(f"Could not find {PANEL_TSV}", file=sys.stderr)
        sys.exit(1)

    print("Generating panel ID list...")
    ids = []
    with open(PANEL_TSV) as f_in:
        header = next(f_in, None)
        for line in f_in:
            if not line.strip():
                continue
            sample = line.split("\t")[0]
            ids.append(sample)

    with open(PANEL_IDS, "w") as f_out:
        for s in ids:
            f_out.write(s + "\n")

    print(f"{len(ids)} IDs written to {PANEL_IDS}")

    for chr_num in range(1, 23):
        chr_str = str(chr_num)
        vcf_in = os.path.join(
            VCF_BASE,
            f"gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr_str}.vcf.bgz",
        )
        out_vcf = os.path.join(
            OUT_BASE,
            f"panel_5super.chr{chr_str}.vcf.gz",
        )

        print(f"\nChr{chr_str}:")
        print(f"Input : {vcf_in}")
        print(f"Output: {out_vcf}")

        if not os.path.isfile(vcf_in):
            print(f"Could not find {vcf_in}, skipping.")
            continue

        # bcftools view with 32 threads
        cmd_view = [
            BCFTOOLS,
            "view",
            "--threads", str(THREADS),
            "-S", PANEL_IDS,
            "-Oz",
            "-o", out_vcf,
            vcf_in,
        ]
        run_cmd(cmd_view)

        # bcftools index with 32 threads
        cmd_index = [
            BCFTOOLS,
            "index",
            "--threads", str(THREADS),
            "-t",
            out_vcf,
        ]
        run_cmd(cmd_index)

    print(f"\nPanel-filtered VCFs generated in: {OUT_BASE}")


if __name__ == "__main__":
    main()


#REFERENCE PANEL STEP 2_SELECTION OF 300 SAMPLES PER SUPERPOPULATION#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Builds a five-superpopulation reference panel
from gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz

Superpopulations:
AFR, AMR, EUR, EAS, SAS
(AMI is grouped into AMR).

Output:
- panel_5superpopulations.tsv : table with sample, project, population, superpopulation
- panel_AFR.ids, panel_AMR.ids, ... : sample ID lists by superpopulation
"""

import os
import sys
import json
import argparse
import pandas as pd
from pandas import json_normalize

# ─────────────────────────────────────────────
# Default parameters
# ─────────────────────────────────────────────
DEFAULT_META = "/mnt/diskrare/ivonb/refamerindios/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz"
DEFAULT_OUTDIR = "/mnt/diskrare/ivonb/refamerindios"
N_PER_SUPERPOP = 300
SEED = 12345

# Column names
SAMPLE_COL   = "s"
PROJECT_COL  = "project"
POP_COL      = "population"
SUPERPOP_COL = "genetic_region"

SUPERPOPS_5 = ["AFR", "AMR", "EUR", "EAS", "SAS"]


# ─────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_cols(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(
            f"Missing columns in metadata: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def map_to_5_superpops(genetic_region):
    """
    Maps genetic_region into five groups:
    AFR, AMR, EUR, EAS, SAS.

    AMI is grouped into AMR.
    Other groups are excluded.
    """
    if pd.isna(genetic_region):
        return None
    gr = str(genetic_region).upper()
    if gr in SUPERPOPS_5:
        return gr
    if gr == "AMI":
        return "AMR"
    return None


def parse_json_field(x):
    """
    Parses JSON fields from gnomAD metadata safely.
    """
    if pd.isna(x):
        return {}
    s = str(x).strip()
    if s == "" or s.upper() == "NA":
        return {}
    try:
        return json.loads(s)
    except Exception:
        eprint(f"Could not parse JSON: {s[:80]}...")
        return {}


# ─────────────────────────────────────────────
# Main logic
# ─────────────────────────────────────────────
def build_panel(meta_path: str, outdir: str):
    print(f"Reading metadata: {meta_path}")
    df_raw = pd.read_csv(
        meta_path,
        sep="\t",
        compression="gzip",
        low_memory=False
    )

    print("Raw metadata columns:", list(df_raw.columns))

    if "hgdp_tgp_meta" not in df_raw.columns:
        raise SystemExit("Could not find column 'hgdp_tgp_meta'.")

    if SAMPLE_COL not in df_raw.columns:
        raise SystemExit(f"Could not find sample ID column '{SAMPLE_COL}'.")

    # Parse hgdp_tgp_meta
    meta_parsed = df_raw["hgdp_tgp_meta"].apply(parse_json_field)
    meta_flat = json_normalize(meta_parsed)

    print("Subcolumns in hgdp_tgp_meta:", meta_flat.columns.tolist())

    cols_from_raw = [SAMPLE_COL]

    for qc_col in ["gnomad_high_quality", "high_quality"]:
        if qc_col in df_raw.columns:
            cols_from_raw.append(qc_col)

    df = pd.concat([df_raw[cols_from_raw], meta_flat], axis=1)

    check_cols(df, [SAMPLE_COL, PROJECT_COL, POP_COL, SUPERPOP_COL])

    print("Columns selected for panel construction:",
          [SAMPLE_COL, PROJECT_COL, POP_COL, SUPERPOP_COL])

    # Optional project filtering
    accepted_projects = ["1000 Genomes", "HGDP"]
    df = df[df[PROJECT_COL].isin(accepted_projects)]

    print(f"Samples after project filtering {accepted_projects}: {len(df)}")

    # Quality filtering
    if "gnomad_high_quality" in df.columns:
        before = len(df)
        df = df[df["gnomad_high_quality"] == True]
        print(f"After gnomad_high_quality filtering: {before} → {len(df)}")

    if "high_quality" in df.columns:
        before = len(df)
        df = df[df["high_quality"] == True]
        print(f"After high_quality filtering: {before} → {len(df)}")

    # Map to five superpopulations
    df["superpop_5"] = df[SUPERPOP_COL].apply(map_to_5_superpops)

    before = len(df)
    df = df[df["superpop_5"].notna()]

    print(f"After five-superpopulation mapping: {before} → {len(df)}")

    df = df[df["superpop_5"].isin(SUPERPOPS_5)]

    print("Sample counts by superpopulation before sampling:")
    print(df["superpop_5"].value_counts())

    # Balanced sampling
    os.makedirs(outdir, exist_ok=True)

    panel_rows = []

    for sp in SUPERPOPS_5:
        sub = df[df["superpop_5"] == sp]

        n_available = len(sub)

        if n_available == 0:
            eprint(f"No samples available for superpopulation {sp}.")
            continue

        n_take = min(N_PER_SUPERPOP, n_available)

        sub_sel = sub.sample(n_take, random_state=SEED)

        out_ids = os.path.join(outdir, f"panel_{sp}.ids")

        sub_sel[SAMPLE_COL].to_csv(out_ids, index=False, header=False)

        print(f"Saved {n_take} IDs for {sp} in {out_ids}")

        for _, row in sub_sel.iterrows():
            panel_rows.append({
                "sample": row[SAMPLE_COL],
                "project": row[PROJECT_COL],
                "population": row[POP_COL],
                "superpop": row["superpop_5"],
            })

    if not panel_rows:
        raise SystemExit("No samples were selected for the final panel.")

    panel_df = pd.DataFrame(panel_rows)

    out_panel = os.path.join(outdir, "panel_5superpopulations.tsv")

    panel_df.to_csv(out_panel, sep="\t", index=False)

    print(f"Combined panel saved in: {out_panel}")

    print("Final sample counts by superpopulation:")
    print(panel_df["superpop"].value_counts())


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Builds a five-superpopulation reference panel "
                    "from HGDP + 1000 Genomes gnomAD metadata."
    )

    p.add_argument(
        "--meta",
        default=DEFAULT_META,
        help=f"Path to metadata file (.tsv.bgz) [default: {DEFAULT_META}]"
    )

    p.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help=f"Output directory [default: {DEFAULT_OUTDIR}]"
    )

    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_panel(args.meta, args.outdir)


#ALL CHROMOSOMES#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sys

BASE = "/mnt/diskrare/ivonb/refamerindios"

PANEL_TSV = os.path.join(BASE, "panel_5superpopulations.tsv")
PANEL_IDS = os.path.join(BASE, "panel_5super.ids")

VCF_BASE = BASE
OUT_BASE = os.path.join(BASE, "vcfs_panel")

BCFTOOLS = "bcftools"
THREADS = 32


def run_cmd(cmd, **kwargs):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)


def main():
    os.makedirs(OUT_BASE, exist_ok=True)

    if not os.path.isfile(PANEL_TSV):
        print(f"Could not find {PANEL_TSV}", file=sys.stderr)
        sys.exit(1)

    print("Generating panel ID list...")
    ids = []

    with open(PANEL_TSV) as f_in:
        header = next(f_in, None)

        for line in f_in:
            if not line.strip():
                continue

            sample = line.split("\t")[0]
            ids.append(sample)

    with open(PANEL_IDS, "w") as f_out:
        for s in ids:
            f_out.write(s + "\n")

    print(f"{len(ids)} IDs written in {PANEL_IDS}")

    for chr_num in range(1, 23):

        chr_str = str(chr_num)

        vcf_in = os.path.join(
            VCF_BASE,
            f"gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr_str}.vcf.bgz",
        )

        out_vcf = os.path.join(
            OUT_BASE,
            f"panel_5super.chr{chr_str}.vcf.gz",
        )

        print(f"\nChr{chr_str}:")
        print(f"Input : {vcf_in}")
        print(f"Output: {out_vcf}")

        if not os.path.isfile(vcf_in):
            print(f"Could not find {vcf_in}, skipping.")
            continue

        # bcftools view using 32 threads
        cmd_view = [
            BCFTOOLS,
            "view",
            "--threads", str(THREADS),
            "-S", PANEL_IDS,
            "-Oz",
            "-o", out_vcf,
            vcf_in,
        ]

        run_cmd(cmd_view)

        # bcftools index using 32 threads
        cmd_index = [
            BCFTOOLS,
            "index",
            "--threads", str(THREADS),
            "-t",
            out_vcf,
        ]

        run_cmd(cmd_index)

    print(f"\nPanel-filtered VCFs generated in: {OUT_BASE}")


if __name__ == "__main__":
    main()
___________________________________________________________________________________________#PLINK#


#PLINK#

#!/usr/bin/env bash
set -euo pipefail

PLINK="/home/rare/programs/Plink/plink"

WORKDIR="/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
cd "$WORKDIR"

# 1. Common SNP intersection between the reference panel and the 18 genomes
cut -f2 panel_5super_autosomes_fixid.bim | sort -u > ids.panel.txt
cut -f2 18genomes_auto_fixid.bim        | sort -u > ids.18g.txt
comm -12 ids.panel.txt ids.18g.txt > ids.common.txt

# 2. Extract common SNPs
$PLINK \
  --bfile panel_5super_autosomes_fixid \
  --extract ids.common.txt \
  --make-bed \
  --out panel_5super_autosomes_common

$PLINK \
  --bfile 18genomes_auto_fixid \
  --extract ids.common.txt \
  --make-bed \
  --out 18genomes_auto_common

# 3. Merge reference panel and 18 genomes
$PLINK \
  --bfile panel_5super_autosomes_common \
  --bmerge 18genomes_auto_common.bed 18genomes_auto_common.bim 18genomes_auto_common.fam \
  --make-bed \
  --out merged_18_1200_common_raw

# 4. Basic quality control
$PLINK \
  --bfile merged_18_1200_common_raw \
  --geno 0.05 \
  --make-bed \
  --out merged_18_1200_common_qc05

# 5. Minor allele frequency filtering (MAF 0.05)
$PLINK \
  --bfile merged_18_1200_common_qc05 \
  --maf 0.05 \
  --make-bed \
  --out merged_18_1200_common_qc05_maf05_nomind

# 6. Linkage disequilibrium pruning
$PLINK \
  --bfile merged_18_1200_common_qc05_maf05_nomind \
  --indep-pairwise 50 5 0.2 \
  --out merged_18_1200_common_qc05_maf05_nomind_ld

$PLINK \
  --bfile merged_18_1200_common_qc05_maf05_nomind \
  --extract merged_18_1200_common_qc05_maf05_nomind_ld.prune.in \
  --make-bed \
  --out merged_18_1200_common_qc05_maf05_nomind_pruned


___________________________________________________________________________________________________________#PCA#
#PCA#

#!/usr/bin/env bash
set -euo pipefail

PLINK="/home/rare/programs/Plink/plink"
WORKDIR="/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"

cd "$WORKDIR"

echo "Running PCA"

$PLINK \
  --bfile merged_18_1200_common_qc05_maf05_nomind_pruned \
  --pca 20 \
  --out PCA_m05

echo "PCA completed"


___________________________________________________________________________________________________________#PLOT_PCA#
#PLOT_PCA#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# === 0) Paths and Configuration ===
BASE_DIR = "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
OUT_DIR = Path("/home/rare/ivon/figpaper")
OUT_DIR.mkdir(parents=True, exist_ok=True)

os.chdir(BASE_DIR)

# === 1) Load PCA eigenvec file ===
evec = pd.read_csv("PCA_m05.eigenvec", sep=r"\s+", header=None)
evec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, 21)]

# === 2) Load superpopulation map ===
pop = pd.read_csv(
    "/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map",
    sep=r"\s+",
    header=None,
    names=["IID", "POP"]
)

# Merge PCA + population information
df = evec.merge(pop, on="IID", how="left")
df["POP"] = df["POP"].fillna("UNK")

# === 3) Identify samples 017/018 ===
mask_own = df["IID"].astype(str).str.startswith("PB000696_")
df.loc[mask_own, "SHORT"] = (
    df.loc[mask_own, "IID"].str.split("_").str[1].str[:3]
)

mask_label = (mask_own & df["SHORT"].isin(["017", "018"])) | df["IID"].astype(str).isin(["017", "018"])

# === 4) Unified color palette ===
color_map = {
    "AMR": "#1F77B4",  # Blue
    "EUR": "#FF7F0E",  # Orange
    "AFR": "#2CA02C",  # Green
    "EAS": "#D62728",  # Red
    "SAS": "#9467BD",  # Purple
    "UNK": "#D3D3D3",  # Gray
}
df["color"] = df["POP"].map(color_map).fillna("#D3D3D3")

# === 5) Plot ===
# Figure size adjusted for effective 400 dpi resolution
fig, ax = plt.subplots(figsize=(10, 8))

legend_pops = [p for p in ["AFR", "AMR", "EUR", "EAS", "SAS"] if p in set(df["POP"])]

# Plot reference population points
for pop_label in legend_pops:
    sub = df[df["POP"] == pop_label]
    if not sub.empty:
        ax.scatter(
            sub["PC1"], sub["PC2"],
            s=30,  # Slightly larger points for 400 dpi
            c=sub["color"],
            label=pop_label,
            alpha=0.7,
            edgecolors='white',
            linewidths=0.2
        )

# Plot UNK samples with low opacity
sub_unk = df[df["POP"] == "UNK"]
if not sub_unk.empty:
    ax.scatter(sub_unk["PC1"], sub_unk["PC2"], s=15, c="#D3D3D3", alpha=0.3, zorder=1)

# Labels for samples 017 and 018
for _, row in df[mask_label].iterrows():
    label_txt = row["SHORT"] if pd.notna(row.get("SHORT")) else str(row["IID"])
    if label_txt.startswith("PB000696_"):
        label_txt = label_txt.split("_")[1][:3]

    ax.text(
        row["PC1"], row["PC2"], label_txt,
        fontsize=10, fontweight='bold',
        ha="center", va="center",
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1),
        color="black", zorder=10
    )

# Axes and visual formatting
ax.set_xlabel("PC1", fontsize=12, fontweight='bold')
ax.set_ylabel("PC2", fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
ax.legend(
    title="Superpopulation",
    title_fontsize=10,
    fontsize=9,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout()

# === 6) Save at 400 DPI ===
out_name = "PCA_PC1_PC2_400DPI_017_018"

plt.savefig(OUT_DIR / f"{out_name}.png", dpi=400, bbox_inches="tight")
plt.savefig(OUT_DIR / f"{out_name}.pdf", bbox_inches="tight")

plt.close()

print(f"Files generated at 400 DPI in: {OUT_DIR}")
________________________________________________________________________________________________________________________#ADMIXURE#

# ADMIXTURE#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
import sys
import re

# =========================
# CONFIG
# =========================
WORKDIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual")

ADMIXTURE = "/usr/local/bin/admixture"
PREFIX = "merged_18_1200_common_qc05_maf05_nomind_pruned"

K_MIN = 1
K_MAX = 6
CV_FOLDS = 10
THREADS = 32
SEED = 12345

# =========================
# CHECKS
# =========================
def check_inputs():
    if not Path(ADMIXTURE).exists():
        sys.exit(f"ADMIXTURE not found at: {ADMIXTURE}")

    for ext in ("bed", "bim", "fam"):
        f = WORKDIR / f"{PREFIX}.{ext}"
        if not f.exists() or f.stat().st_size == 0:
            sys.exit(f"Missing or empty file: {f}")

def run_one_k(k: int) -> Path:
    log_file = WORKDIR / f"cv_K{k}.log"
    cmd = [
        ADMIXTURE,
        f"--cv={CV_FOLDS}",
        f"-j{THREADS}",
        f"--seed={SEED}",
        f"{PREFIX}.bed",
        str(k),
    ]

    print("\n======================================")
    print(f"Running ADMIXTURE K={k}")
    print("CMD:", " ".join(cmd))
    print("LOG:", log_file.name)
    print("======================================")

    # Execute and simultaneously print/save output
    with open(log_file, "w") as log:
        proc = subprocess.Popen(
            cmd,
            cwd=WORKDIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            log.write(line)

    rc = proc.wait()
    if rc != 0:
        sys.exit(f"ADMIXTURE failed at K={k} (return code={rc}). Check {log_file}")

    return log_file

def extract_cv_error(log_file: Path) -> str | None:
    # Search for lines such as: "CV error (K=4): 0.38712"
    pat = re.compile(r"CV error\s*\(K=\d+\)\s*:\s*([0-9.eE+-]+)")
    with open(log_file) as f:
        for line in f:
            m = pat.search(line)
            if m:
                return m.group(1)
    return None

def main():
    check_inputs()

    print("WORKDIR:", WORKDIR)
    print("PREFIX :", PREFIX)
    print("ADMIXTURE:", ADMIXTURE)
    print(f"Parameters: --cv={CV_FOLDS} -j{THREADS} --seed={SEED}  |  K={K_MIN}..{K_MAX}")

    results = []

    for k in range(K_MIN, K_MAX + 1):
        logf = run_one_k(k)
        cv = extract_cv_error(logf)

        if cv is None:
            print(f"'CV error' was not found in {logf.name}")
        else:
            results.append((k, cv))

    # Save summary table
    out_tsv = WORKDIR / "cv_error.tsv"

    with open(out_tsv, "w") as out:
        out.write("K\tCV_error\n")
        for k, cv in results:
            out.write(f"{k}\t{cv}\n")

    print("\nAnalysis completed.")
    print("Logs: cv_K1.log ... cv_K6.log (same directory)")
    print("Summary:", out_tsv)

    # Print CV errors
    if results:
        print("\nCV errors:")
        for k, cv in results:
            print(f"K={k}\tCV={cv}")

if __name__ == "__main__":
    main()
    
    
    
    ___________________________________________________________________________________________________#PLOT_ADMIXURE#

#PLOT ADMIXTURE#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# CONFIG
# =========================
BASE_DIR = "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
OUT_DIR  = "/home/rare/ivon/figpaper"
PREFIX   = "merged_18_1200_common_qc05_maf05_nomind_pruned"
MAP_FILE = "/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map"

K_MIN = 1
K_MAX = 6
K_MAIN = 4  # Main barplot + pie charts

# Equal-width blocks
BLOCKS_ORDER = ["AFR", "EUR", "EAS", "AMR", "COLM", "COLF"]
BLOCK_WIDTH        = 200
GAP_BETWEEN_BLOCKS = 60

# Project blocks
KEEP_COLM = {"017"}
KEEP_COLF = {"018"}

# If manual component labels are desired (length must equal K_MAIN), uncomment:
# COMP_LABELS_K4 = ["EUR", "AMR", "EAS", "AFR"]
COMP_LABELS_K4 = None  # If None, labels will be inferred from AFR/EUR/EAS/AMR blocks

os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# HELPERS
# =========================
def check_file(path: str):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Missing or empty file: {path}")

def load_fam(path: str) -> pd.DataFrame:
    fam = pd.read_csv(
        path, sep=r"\s+", header=None,
        names=["FID","IID","PID","MID","SEX","PHENO"]
    )
    return fam[["FID","IID"]].copy()

def load_map(path: str) -> pd.DataFrame:
    m = pd.read_csv(path, sep=r"\s+", header=None, names=["IID","SUPERPOP"])
    m["SUPERPOP"] = m["SUPERPOP"].astype(str).str.strip()
    return m

def load_Q(path: str) -> np.ndarray:
    Q = np.loadtxt(path)
    if Q.ndim == 1:
        Q = Q.reshape(-1, 1)
    return Q

def is_project(iid: str) -> bool:
    return str(iid).startswith("PB000696_")

def short_id(iid: str) -> str:
    # PB000696_017C_... -> 017
    try:
        return str(iid).split("_")[1][:3]
    except Exception:
        return str(iid)

def build_df_and_layout(fam: pd.DataFrame, pop: pd.DataFrame):

    """
    Build dataframe with IID, SHORT, GROUP and block layout:
    - x_positions (left border)
    - widths per row (equal-width blocks)
    - block_centers, block_ranges, used_blocks
    """

    df = fam.merge(pop, on="IID", how="left")
    df["SHORT"] = df["IID"].apply(short_id)

    # GROUP: default SUPERPOP; if project sample and SHORT matches -> COLM/COLF
    df["GROUP"] = df["SUPERPOP"]
    df.loc[df["IID"].apply(is_project) & df["SHORT"].isin(KEEP_COLM), "GROUP"] = "COLM"
    df.loc[df["IID"].apply(is_project) & df["SHORT"].isin(KEEP_COLF), "GROUP"] = "COLF"

    # Filtering
    df = df[df["GROUP"].notna()].copy()
    df = df[df["GROUP"] != "UNK"].copy()
    df = df[df["GROUP"].isin(BLOCKS_ORDER)].copy()

    # Sort by blocks; within each block sort by SHORT
    df["BLOCK_ORDER"] = pd.Categorical(df["GROUP"], categories=BLOCKS_ORDER, ordered=True)
    df.sort_values(["BLOCK_ORDER", "SHORT"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Count samples per block
    sample_counts = {b: int((df["GROUP"] == b).sum()) for b in BLOCKS_ORDER}
    sample_counts = {b: n for b, n in sample_counts.items() if n > 0}

    # Width per block
    bar_widths_block = {b: (BLOCK_WIDTH / sample_counts[b]) for b in sample_counts}

    # x positions
    x_positions = np.zeros(len(df), dtype=float)
    block_centers = []
    block_ranges = []
    used_blocks = []

    current_x = 0.0
    idx = 0

    for b in BLOCKS_ORDER:
        n = sample_counts.get(b, 0)

        if n == 0:
            continue

        if used_blocks:
            current_x += GAP_BETWEEN_BLOCKS

        start = current_x
        w = bar_widths_block[b]

        x_positions[idx:idx+n] = start + np.arange(n) * w
        end = start + BLOCK_WIDTH

        block_centers.append(start + BLOCK_WIDTH / 2)
        block_ranges.append((start, end))
        used_blocks.append(b)

        current_x = end
        idx += n

    widths = np.array([bar_widths_block[g] for g in df["GROUP"].tolist()], dtype=float)

    return df, x_positions, widths, block_centers, block_ranges, used_blocks

def align_Q_to_df(Q: np.ndarray, fam: pd.DataFrame, df: pd.DataFrame) -> np.ndarray:
    fam_pos = {iid: i for i, iid in enumerate(fam["IID"].astype(str).tolist())}
    rows = np.array([fam_pos[str(iid)] for iid in df["IID"].astype(str).tolist()], dtype=int)
    Qk = Q[rows, :].astype(float)
    Qk = Qk / np.clip(Qk.sum(axis=1, keepdims=True), 1e-12, None)
    return Qk

def add_block_labels_and_separators(ax, block_centers, block_ranges, used_blocks):

    # Labels
    for c, lab in zip(block_centers, used_blocks):
        ax.text(
            c, 1.02, lab,
            ha="center", va="bottom",
            transform=ax.get_xaxis_transform(),
            fontsize=11, fontweight="bold"
        )

    # Separators
    for j in range(1, len(block_ranges)):
        sep_pos = block_ranges[j-1][1] + (GAP_BETWEEN_BLOCKS / 2)
        ax.axvline(sep_pos, color="gray", linestyle="--", linewidth=1)

def infer_component_labels_from_blocks(df: pd.DataFrame, Qk: np.ndarray, ref_blocks=("AFR","EUR","EAS","AMR")):

    """
    Infer ancestry component labels from reference blocks.
    """

    K = Qk.shape[1]
    block_means = {}

    for b in ref_blocks:
        idx = df.index[df["GROUP"] == b].to_list()

        if len(idx) == 0:
            continue

        block_means[b] = Qk[idx, :].mean(axis=0)

    if len(block_means) < 2:
        return [f"C{j+1}" for j in range(K)]

    labels = [None] * K
    used = set()

    clarity = []

    for j in range(K):
        scores = [(b, block_means[b][j]) for b in block_means]
        scores.sort(key=lambda x: x[1], reverse=True)

        top1 = scores[0][1]
        top2 = scores[1][1] if len(scores) > 1 else 0.0

        clarity.append((j, top1 - top2, scores))

    clarity.sort(key=lambda x: x[1], reverse=True)

    for j, _, scores in clarity:

        chosen = None

        for b, _v in scores:
            if b not in used:
                chosen = b
                break

        if chosen is None:
            chosen = scores[0][0]

        labels[j] = chosen
        used.add(chosen)

    if K > len(ref_blocks):

        counts = {}

        for lab in labels:
            counts[lab] = counts.get(lab, 0) + 1

        seen = {}
        new_labels = []

        for lab in labels:

            if counts.get(lab, 0) > 1:
                seen[lab] = seen.get(lab, 0) + 1
                new_labels.append(f"{lab}{seen[lab]}")
            else:
                new_labels.append(lab)

        labels = new_labels

    return labels

def get_component_colors(K: int):

    """
    Consistent color palette for barplots and pie charts.
    """

    cmap = plt.get_cmap("tab10") if K <= 10 else plt.get_cmap("tab20")
    return [cmap(i) for i in range(K)]

# =========================
# FIGURES
# =========================
def plot_multipanel_Ks(df, x_positions, widths, fam, block_centers, block_ranges, used_blocks, out_png, K_min, K_max):

    nrows = K_max - K_min + 1

    fig, axes = plt.subplots(
        nrows=nrows, ncols=1,
        figsize=(24, 2.4 * nrows),
        sharex=True
    )

    if nrows == 1:
        axes = [axes]

    for i, K in enumerate(range(K_min, K_max + 1)):

        q_path = os.path.join(BASE_DIR, f"{PREFIX}.{K}.Q")
        check_file(q_path)

        Q = load_Q(q_path)
        Qk = align_Q_to_df(Q, fam, df)

        colors = get_component_colors(Qk.shape[1])

        ax = axes[i]

        bottom = np.zeros(Qk.shape[0], dtype=float)

        for kcol in range(Qk.shape[1]):

            ax.bar(
                x_positions, Qk[:, kcol],
                bottom=bottom,
                width=widths,
                align="edge",
                linewidth=0,
                color=colors[kcol]
            )

            bottom += Qk[:, kcol]

        ax.set_ylim(0, 1)
        ax.set_ylabel(f"K={K}", rotation=0, labelpad=30, va="center")
        ax.set_yticks([0, 0.5, 1.0])

        add_block_labels_and_separators(ax, block_centers, block_ranges, used_blocks)

        ax.set_xticks([])

    plt.tight_layout()
    plt.savefig(out_png, dpi=400, bbox_inches="tight")
    plt.close()

    print("Saved:", out_png)

def plot_barplot_oneK_K4(df, x_positions, widths, Qk, K, out_png, out_pdf, block_centers, block_ranges, used_blocks, comp_labels):

    """
    Main K=4 barplot.
    """

    fig, ax = plt.subplots(figsize=(24, 3.2))

    colors = get_component_colors(Qk.shape[1])

    bottom = np.zeros(Qk.shape[0], dtype=float)

    for kcol in range(Qk.shape[1]):

        ax.bar(
            x_positions, Qk[:, kcol],
            bottom=bottom,
            width=widths,
            align="edge",
            linewidth=0,
            color=colors[kcol]
        )

        bottom += Qk[:, kcol]

    ax.set_ylim(0, 1)
    ax.set_ylabel("Ancestry proportion")
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xticks([])

    add_block_labels_and_separators(ax, block_centers, block_ranges, used_blocks)

    if comp_labels is not None and len(comp_labels) == Qk.shape[1]:

        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i]) for i in range(Qk.shape[1])]

        ax.legend(
            handles, comp_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.18),
            ncol=min(6, len(comp_labels)),
            frameon=False
        )

    plt.tight_layout()

    plt.savefig(out_png, dpi=400, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")

    plt.close()

    print("Saved:", out_png)
    print("Saved:", out_pdf)

def plot_pies_K4_017_018(df, Qk, K, out_png, comp_labels):

    """
    Pie charts for samples 017 (COLM) and 018 (COLF).
    """

    df_colm = df[df["GROUP"] == "COLM"]
    df_colf = df[df["GROUP"] == "COLF"]

    if df_colm.shape[0] == 0 or df_colf.shape[0] == 0:
        print("COLM and/or COLF samples were not found.")
        return

    i_colm = int(df_colm.index[0])
    i_colf = int(df_colf.index[0])

    if comp_labels is None or len(comp_labels) != Qk.shape[1]:
        comp_labels = [f"C{j+1}" for j in range(Qk.shape[1])]

    colors = get_component_colors(Qk.shape[1])

    def autopct_fmt(pct):
        return f"{pct:.1f}%" if pct >= 1 else ""

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))

    wedges0, _texts0, _auto0 = axes[0].pie(
        Qk[i_colm, :],
        startangle=90,
        counterclock=False,
        colors=colors,
        autopct=autopct_fmt,
        textprops={"fontsize": 9},
        wedgeprops={"linewidth": 0.6, "edgecolor": "white"}
    )

    axes[0].set_title("017 (COLM)")
    axes[0].set_aspect("equal")

    wedges1, _texts1, _auto1 = axes[1].pie(
        Qk[i_colf, :],
        startangle=90,
        counterclock=False,
        colors=colors,
        autopct=autopct_fmt,
        textprops={"fontsize": 9},
        wedgeprops={"linewidth": 0.6, "edgecolor": "white"}
    )

    axes[1].set_title("018 (COLF)")
    axes[1].set_aspect("equal")

    fig.legend(
        wedges0, comp_labels,
        loc="lower center",
        ncol=min(6, len(comp_labels)),
        frameon=False,
        bbox_to_anchor=(0.5, -0.02)
    )

    plt.tight_layout()
    plt.savefig(out_png, dpi=400, bbox_inches="tight")

    plt.close()

    print("Saved:", out_png)

# =========================
# MAIN
# =========================
def main():

    os.chdir(BASE_DIR)

    fam_path = os.path.join(BASE_DIR, f"{PREFIX}.fam")

    check_file(fam_path)
    check_file(MAP_FILE)

    fam = load_fam(fam_path)
    pop = load_map(MAP_FILE)

    df, x_positions, widths, block_centers, block_ranges, used_blocks = build_df_and_layout(fam, pop)

    # ---------- 1) Multipanel K=1..6 ----------
    out_multi = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MIN}-K{K_MAX}.equalblocks_COLM_COLF.png")

    plot_multipanel_Ks(
        df, x_positions, widths, fam,
        block_centers, block_ranges, used_blocks,
        out_multi, K_MIN, K_MAX
    )

    # ---------- 2) Main K=4 barplot ----------
    q_main = os.path.join(BASE_DIR, f"{PREFIX}.{K_MAIN}.Q")

    check_file(q_main)

    Q = load_Q(q_main)
    Qk_main = align_Q_to_df(Q, fam, df)

    if COMP_LABELS_K4 is not None:

        if len(COMP_LABELS_K4) != Qk_main.shape[1]:
            sys.exit(f"COMP_LABELS_K4 length is {len(COMP_LABELS_K4)} but K={Qk_main.shape[1]}")

        comp_labels = COMP_LABELS_K4

    else:
        comp_labels = infer_component_labels_from_blocks(
            df, Qk_main,
            ref_blocks=("AFR","EUR","EAS","AMR")
        )

    out_k4_png = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MAIN}.barplot_equalblocks_COLM_COLF.png")
    out_k4_pdf = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MAIN}.barplot_equalblocks_COLM_COLF.pdf")

    plot_barplot_oneK_K4(
        df, x_positions, widths, Qk_main, K_MAIN,
        out_k4_png, out_k4_pdf,
        block_centers, block_ranges, used_blocks,
        comp_labels=comp_labels
    )

    # ---------- 3) K=4 pie charts ----------
    out_pies = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MAIN}.pies_017_COLM_018_COLF.png")

    plot_pies_K4_017_018(
        df,
        Qk_main,
        K_MAIN,
        out_pies,
        comp_labels=comp_labels
    )

    print("\nAnalysis completed. Review output files in:", OUT_DIR)
    print("   - K1..K6 multipanel plot")
    print("   - K4 barplot")
    print("   - K4 pie charts")

if __name__ == "__main__":
    main()
_________________________________________________________________________________________________________________________
#RFMIX#

Phased chromosome-specific VCFs
↓
10_rfmix_intersect_harmonize_merge_panel_query.sh
↓
Merged panel1199 + sample
↓
11_prepare_rfmix_inputs_query_ref_labels_maps.py
↓
chrN.query.vcf.gz
chrN.ref.vcf.gz
chrN.superpopulation_labels.txt
chrN.snp_locations
↓
12_run_rfmix_all_samples.sh


#PHASING#

#Phase genomes#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from pathlib import Path

# ───────────────────────────────
# CONFIGURATION FOR THE 18 GENOMES
# ───────────────────────────────
VCF_DIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes")
BAM_DIR = Path("/mnt/diskrare/arlenb/08/aligned_reads/hg38")
OUT_DIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/18genfased")
REF = "/mnt/diskrare/arlenb/reference/hg38.fasta"  # reference genome path
WHATSHAP = "/home/rare/.local/bin/whatshap"

# Create output directory
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Log file for issues
log = open("log_phasing_18.txt", "w")

# Exact list of 18 samples
SAMPLES = [
    "08_1_A01_bc2043_001P",
    "08_1_A01_bc2044_002P",
    "08_1_B01_bc2045_003P",
    "08_1_C01_bc2046_004P",
    "08_1_A01_bc2047_005P",
    "08_1_B01_bc2048_006P",
    "08_1_C01_bc2049_007P",
    "08_1_C01_bc2050_008D",
    "08_1_B01_bc2051_009D",
    "08_1_D01_bc2052_010D",
    "08_1_D01_bc2053_011D",
    "08_1_D01_bc2054_012D",
    "08_1_A01_bc2055_013A",
    "08_1_A01_bc2056_014A",
    "08_1_B01_bc2057_015A",
    "08_1_B01_bc2058_016D",
    "08_1_C01_bc2059_017C",
    "08_1_D01_bc2060_018C"
]

print("Starting phasing of the 18 genomes using WhatsHap...\n")

for sample in SAMPLES:
    vcf_path = VCF_DIR / f"{sample}.autosomes.vcf.gz"
    bam_path = BAM_DIR / f"{sample}.bam"
    bai_path = BAM_DIR / f"{sample}.bam.bai"

    print(f"Processing: {sample}")

    if not vcf_path.exists():
        msg = f"Missing VCF: {vcf_path}"
        print(msg)
        log.write(msg + "\n")
        continue

    if not bam_path.exists():
        msg = f"Missing BAM: {bam_path} (verify whether the .bam.bai index exists)"
        print(msg)
        log.write(msg + "\n")
        continue

    if not bai_path.exists():
        print("Missing .bai index → creating index...")
        subprocess.run(["samtools", "index", str(bam_path)], check=True)

    phased_vcf = OUT_DIR / f"{sample}.autosomes.phased.vcf.gz"

    cmd = [
        WHATSHAP, "phase",
        "--reference", REF,
        "--indels",
        "--ignore-read-groups",
        "-o", str(phased_vcf),
        str(vcf_path),
        str(bam_path)
    ]

    try:
        subprocess.run(cmd, check=True)
        print("   Phasing completed")

        # Index phased VCF
        subprocess.run(["tabix", "-p", "vcf", str(phased_vcf)], check=True)

        # Optional: retain only fully phased variants (both alleles phased with |)
        final_vcf = OUT_DIR / f"{sample}.autosomes.fully_phased.vcf.gz"

        subprocess.run([
            "bcftools", "view",
            "-i", 'GT~"|"',
            "-Oz", "-o", str(final_vcf),
            str(phased_vcf)
        ], check=True)

        subprocess.run(["tabix", "-p", "vcf", str(final_vcf)], check=True)

        print(f"Completed: {sample} → {final_vcf.name}\n")

    except subprocess.CalledProcessError as e:
        msg = f"Error in {sample}: {e}"
        print(msg)
        log.write(msg + "\n")

log.close()

print("All samples processed")
print(f"Phased VCFs saved in: {OUT_DIR}")
print("Check log_phasing_18.txt for potential issues.")

___________________________________________________________________________________________________________________________________#SPLIT#

#SPLIT PHASED 18 GENOMES BY CHROMOSOME#

#!/bin/bash
set -euo pipefail

# Directory containing fully phased VCFs
PHASED_DIR="18genfased"
SPLIT_DIR="18gen_by_chr"

mkdir -p $SPLIT_DIR
cd $PHASED_DIR

echo "Starting chromosome-level split for the 18 genomes..."

# Mapping: long sample name → short code (001 to 018)
declare -A CODES
CODES["08_1_A01_bc2043_001P"]=001
CODES["08_1_A01_bc2044_002P"]=002
CODES["08_1_B01_bc2045_003P"]=003
CODES["08_1_C01_bc2046_004P"]=004
CODES["08_1_A01_bc2047_005P"]=005
CODES["08_1_B01_bc2048_006P"]=006
CODES["08_1_C01_bc2049_007P"]=007
CODES["08_1_C01_bc2050_008D"]=008
CODES["08_1_B01_bc2051_009D"]=009
CODES["08_1_D01_bc2052_010D"]=010
CODES["08_1_D01_bc2053_011D"]=011
CODES["08_1_D01_bc2054_012D"]=012
CODES["08_1_A01_bc2055_013A"]=013
CODES["08_1_A01_bc2056_014A"]=014
CODES["08_1_B01_bc2057_015A"]=015
CODES["08_1_B01_bc2058_016D"]=016
CODES["08_1_C01_bc2059_017C"]=017
CODES["08_1_D01_bc2060_018C"]=018

for full_vcf in *.autosomes.fully_phased.vcf.gz; do
    [[ -f "$full_vcf" ]] || continue
    
    base=${full_vcf%.autosomes.fully_phased.vcf.gz}
    code=${CODES[$base]:-UNKNOWN}
    
    if [[ "$code" == "UNKNOWN" ]]; then
        echo "Warning: Unrecognized sample $base → skipping"
        continue
    fi
    
    sample_dir="../$SPLIT_DIR/$code"
    mkdir -p "$sample_dir"
    
    echo "Splitting genome $code ($base)..."
    
    # Chromosomes 1-22
    for chr in {1..22}; do
        out="${sample_dir}/chr${chr}.vcf.gz"
        bcftools view -Oz -o "$out" "$full_vcf" $chr && tabix -p vcf "$out" &
    done
    
    # Chromosomes X and Y (if present in the VCF; does not fail if absent)
    bcftools view -Oz -o "${sample_dir}/chrX.vcf.gz" "$full_vcf" X 2>/dev/null && tabix -p vcf "${sample_dir}/chrX.vcf.gz" || true &
    bcftools view -Oz -o "${sample_dir}/chrY.vcf.gz" "$full_vcf" Y 2>/dev/null && tabix -p vcf "${sample_dir}/chrY.vcf.gz" || true &
    
    wait  # wait until processing for the current genome is completed
    echo "Genome $code successfully split and indexed"
done

echo "All 18 genomes were successfully split by chromosome into $SPLIT_DIR"
___________________________________________________________________________________________________________________________________________________________#MERGED#
#panel1199 chrN
+
sample chrN
↓
common site intersection
↓
exact CHROM POS REF ALT harmonization
↓
final panel
↓
panel + sample merge
↓
merged VCF with 1200 samples
___________________________________________________________________________________________________________________________________________________________#MERGED 017 (COLM)#


#MERGED 017#

#!/usr/bin/env bash
set -euo pipefail

# =========================
# PATHS FOR SAMPLE 017
# =========================
SAMPLE="017"

BASE_SAMPLE="/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/18genfased/${SAMPLE}"
OUTBASE="${BASE_SAMPLE}/merged"
PANELDIR="/mnt/diskrare/ivonb/refamerindios/panel_hgdp1kg_1200"

mkdir -p "$OUTBASE"

for chr in {1..22}; do
  echo "==================== ${SAMPLE} chr${chr} ===================="

  odir="${OUTBASE}/isec_chr${chr}"
  mkdir -p "$odir"
  cd "$odir"

  PANEL="${PANELDIR}/hgdp1kgp_chr${chr}.SNV_ONLY.BIAL.panel1199.vcf.gz"
  QUERY="${BASE_SAMPLE}/chr${chr}.vcf.gz"

  [[ -s "$PANEL" ]] || { echo "PANEL file not found: $PANEL"; exit 1; }
  [[ -s "$QUERY" ]] || { echo "QUERY file not found: $QUERY"; exit 1; }

  # Remove previous outputs
  rm -f chr${chr}.${SAMPLE}.common.vcf.gz chr${chr}.${SAMPLE}.common.vcf.gz.tbi \
        chr${chr}.panel1199.final.vcf.gz chr${chr}.panel1199.final.vcf.gz.tbi \
        chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz \
        chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz.tbi \
        sites.pos sites${SAMPLE}.4col.sorted panel.header.vcf panel.body.filtered \
        0000.vcf.gz 0000.vcf.gz.tbi 0001.vcf.gz 0001.vcf.gz.tbi README.txt sites.txt

  # 1. Intersection between panel and sample 017
  bcftools isec -n=2 -w1 -O z -p "$odir" "$PANEL" "$QUERY"

  [[ -s sites.txt ]] || { echo "sites.txt was not generated"; exit 1; }
  nsites=$(wc -l < sites.txt)
  echo "sites.txt: $nsites"

  # 2. Filter sample 017 by common positions
  cut -f1,2 sites.txt > sites.pos

  bcftools view -T sites.pos -Oz \
    -o chr${chr}.${SAMPLE}.common.vcf.gz \
    "$QUERY"

  tabix -f -p vcf chr${chr}.${SAMPLE}.common.vcf.gz

  nquery=$(bcftools view -H chr${chr}.${SAMPLE}.common.vcf.gz | wc -l)
  echo "${SAMPLE}.common: $nquery"

  [[ "$nquery" -eq "$nsites" ]] || { echo "${SAMPLE}.common != sites.txt"; exit 1; }

  # 3. Extract exact CHROM POS REF ALT sites from sample 017
  bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' chr${chr}.${SAMPLE}.common.vcf.gz \
    | sort -u > sites${SAMPLE}.4col.sorted

  # 4. Build final panel with exact allele matching
  bcftools view -h "$PANEL" > panel.header.vcf

  bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT[\t%GT]\n' "$PANEL" \
  | awk 'BEGIN{FS=OFS="\t"}
      NR==FNR { key[$1"\t"$2"\t"$3"\t"$4]=1; next }
      {
        k=$1"\t"$2"\t"$3"\t"$4
        if(k in key){
          for(i=5;i<=NF;i++) printf "%s%s", $i, (i==NF?ORS:OFS)
        }
      }' sites${SAMPLE}.4col.sorted - \
  > panel.body.filtered

  cat panel.header.vcf panel.body.filtered | bgzip -c > chr${chr}.panel1199.final.vcf.gz
  tabix -f -p vcf chr${chr}.panel1199.final.vcf.gz

  npanel=$(bcftools view -H chr${chr}.panel1199.final.vcf.gz | wc -l)
  nsamp_panel=$(bcftools query -l chr${chr}.panel1199.final.vcf.gz | wc -l)

  echo "panel.final: $npanel | samples(panel)=$nsamp_panel"

  [[ "$npanel" -eq "$nsites" ]] || { echo "panel.final != sites.txt"; exit 1; }
  [[ "$nsamp_panel" -eq 1199 ]] || { echo "panel.final samples != 1199"; exit 1; }

  # 5. Merge panel + sample 017
  bcftools merge -m none -Oz \
    -o chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz \
    chr${chr}.panel1199.final.vcf.gz \
    chr${chr}.${SAMPLE}.common.vcf.gz

  tabix -f -p vcf chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz

  nmerged=$(bcftools view -H chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz | wc -l)
  nsamp_merged=$(bcftools query -l chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz | wc -l)

  echo "merged: $nmerged | samples(merged)=$nsamp_merged"

  [[ "$nmerged" -eq "$nsites" ]] || { echo "merged != sites.txt"; exit 1; }
  [[ "$nsamp_merged" -eq 1200 ]] || { echo "merged samples != 1200"; exit 1; }

  echo "chr${chr} OK"
done

echo "======== FINAL SUMMARY ${SAMPLE} ========"

for chr in {1..22}; do
  dir="${OUTBASE}/isec_chr${chr}"
  nsites=$(wc -l < "${dir}/sites.txt" 2>/dev/null || echo 0)
  nmerged=$(bcftools view -H "${dir}/chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz" 2>/dev/null | wc -l)
  echo "chr${chr} sites=${nsites} merged=${nmerged}"
done



___________________________________________________________________________________________________________________________________________________________#MERGED 018 (COLF)#
#MERGED GENOME_018#

#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
from pathlib import Path

BASE_018 = Path("/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/18genfased/018")
OUTBASE = BASE_018 / "with_ref_panel1199"
PANELDIR = Path("/mnt/diskrare/ivonb/refamerindios/panel_hgdp1kg_1200")

CHRS = list(range(1, 23))

def run(cmd, cwd=None):
    """Run a shell command and raise an error if it fails."""
    print("  $", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)

def capture(cmd, cwd=None) -> str:
    """Run a shell command and capture stdout."""
    return subprocess.check_output(cmd, cwd=cwd, text=True).strip()

def count_sites_txt(path: Path) -> int:
    return sum(1 for _ in path.open())

def count_vcf_records(vcf_gz: Path) -> int:
    out = capture(["bash", "-lc", f"bcftools view -H {vcf_gz} | wc -l"])
    return int(out)

def count_samples(vcf_gz: Path) -> int:
    out = capture(["bash", "-lc", f"bcftools query -l {vcf_gz} | wc -l"])
    return int(out)

def main():
    OUTBASE.mkdir(parents=True, exist_ok=True)

    for chr_ in CHRS:
        print(f"==================== chr{chr_} ====================")
        odir = OUTBASE / f"isec_chr{chr_}"
        odir.mkdir(parents=True, exist_ok=True)

        panel = PANELDIR / f"hgdp1kgp_chr{chr_}.SNV_ONLY.BIAL.panel1199.vcf.gz"
        query = BASE_018 / f"chr{chr_}.vcf.gz"

        if not panel.exists() or panel.stat().st_size == 0:
            print(f"PANEL file not found: {panel}", file=sys.stderr)
            sys.exit(1)

        if not query.exists() or query.stat().st_size == 0:
            print(f"QUERY file not found: {query}", file=sys.stderr)
            sys.exit(1)

        # Remove previous outputs
        to_remove = [
            odir / f"chr{chr_}.018.common.vcf.gz",
            odir / f"chr{chr_}.018.common.vcf.gz.tbi",
            odir / f"chr{chr_}.panel1199.final.vcf.gz",
            odir / f"chr{chr_}.panel1199.final.vcf.gz.tbi",
            odir / f"chr{chr_}.panel1199_plus_018.common.merged.vcf.gz",
            odir / f"chr{chr_}.panel1199_plus_018.common.merged.vcf.gz.tbi",
            odir / "sites.pos",
            odir / "sites018.4col.sorted",
            odir / "panel.header.vcf",
            odir / "panel.body.filtered",
            odir / "0000.vcf.gz",
            odir / "0000.vcf.gz.tbi",
            odir / "0001.vcf.gz",
            odir / "0001.vcf.gz.tbi",
            odir / "README.txt",
            odir / "sites.txt",
        ]

        for p in to_remove:
            if p.exists():
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()

        # 1) Intersection
        run([
            "bcftools", "isec", "-n=2", "-w1", "-O", "z",
            "-p", str(odir),
            str(panel), str(query)
        ])

        sites_txt = odir / "sites.txt"

        if not sites_txt.exists() or sites_txt.stat().st_size == 0:
            print("sites.txt was not generated", file=sys.stderr)
            sys.exit(1)

        nsites = count_sites_txt(sites_txt)
        print(f"sites.txt: {nsites}")

        # 2) Filter sample 018 by positions
        sites_pos = odir / "sites.pos"
        run(["bash", "-lc", f"cut -f1,2 {sites_txt} > {sites_pos}"])

        vcf_018_common = odir / f"chr{chr_}.018.common.vcf.gz"

        run([
            "bcftools", "view",
            "-T", str(sites_pos),
            "-Oz",
            "-o", str(vcf_018_common),
            str(query)
        ])

        run(["tabix", "-f", "-p", "vcf", str(vcf_018_common)])

        n018 = count_vcf_records(vcf_018_common)
        print(f"018.common: {n018}")

        if n018 != nsites:
            print("018.common != sites.txt", file=sys.stderr)
            sys.exit(1)

        print(f"chr{chr_} OK")

if __name__ == "__main__":
    main()



__________________________________________________________________________________________________________________________________#RFMIX FILES#

#RFMIX FILE PREPARATION#

_________________________________________________________________________________________________________________________________#017 (COLM)#  
#query_ref_locations_017#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from cyvcf2 import VCF

# ================= CONFIG =================
CHRS = range(1, 23)

BASE = Path(
    "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/"
    "18genomes/18genfased/017/merged"
)

PANEL_IDS = Path(
    "/mnt/diskrare/ivonb/refamerindios/panel_hgdp1kg_1200/panel_1199.ids"
)

MAPDIR = Path("/home/rare/ivon/outputs/rfmix/rfmix_map")
# ==========================================


def run(cmd):
    print(">", " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_genetic_map(chr_str):
    gmap = MAPDIR / f"genetic_map_chr{chr_str}_rfmix.txt"

    if not gmap.exists():
        raise FileNotFoundError(gmap)

    df = pd.read_csv(
        gmap,
        sep=r"\s+",
        header=None,
        names=["pos", "cM", "id"]
    ).dropna()

    df = df.sort_values("pos")
    return df["pos"].values, df["cM"].values


def make_snp_locations(vcf_path, chr_str, out_path):
    pos_arr, cm_arr = load_genetic_map(chr_str)
    vcf = VCF(str(vcf_path))

    out = []

    for v in vcf:
        if not v.is_snp or len(v.ALT) != 1:
            continue

        pos = v.POS
        cm = np.interp(pos, pos_arr, cm_arr)
        vid = f"{v.CHROM}:{pos}:{v.REF}:{v.ALT[0]}"
        out.append(f"{vid}\t{pos}\t{cm}\n")

    if not out:
        raise RuntimeError(f"No SNPs found in {vcf_path}")

    out_path.write_text("".join(out))


def main():
    for chrn in CHRS:
        chr_str = str(chrn)
        idir = BASE / f"isec_chr{chr_str}"
        merged = idir / f"chr{chr_str}.panel1199_plus_017.common.merged.vcf.gz"

        if not merged.exists():
            print(f"chr{chr_str}: merged file not found, skipping")
            continue

        print(f"\n=========== chr{chr_str} ===========")

        # Detect query sample
        res = subprocess.check_output(
            ["bcftools", "query", "-l", str(merged)],
            text=True
        ).strip().splitlines()

        query_id = res[-1]
        print("Query ID:", query_id)

        # REF
        ref_vcf = idir / f"chr{chr_str}.ref.vcf.gz"

        run([
            "bcftools", "view",
            "-S", str(PANEL_IDS),
            "-Oz", "-o", str(ref_vcf),
            str(merged)
        ])

        run(["tabix", "-p", "vcf", str(ref_vcf)])

        # QUERY
        query_vcf = idir / f"chr{chr_str}.query.vcf.gz"

        run([
            "bcftools", "view",
            "-s", query_id,
            "-Oz", "-o", str(query_vcf),
            str(merged)
        ])

        run(["tabix", "-p", "vcf", str(query_vcf)])

        # SNP locations
        snp_loc = idir / f"chr{chr_str}.snp_locations"
        make_snp_locations(merged, chr_str, snp_loc)

        print(f"chr{chr_str} completed: ref | query | snp_locations")


if __name__ == "__main__":
    main()
_____________________________________________________________________________________________________________________________________________________#018 (COLF)#
#018 (COLF)#
#.query_ref_location_018#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from cyvcf2 import VCF

# ================= CONFIG =================
CHRS = range(1, 23)

# Base directory for sample 018 (contains isec_chr1..isec_chr22)
BASE = Path(
    "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/"
    "18genomes/18genfased/018/with_ref_panel1199"
)

# Panel IDs list (1199)
PANEL_IDS = Path(
    "/mnt/diskrare/ivonb/refamerindios/panel_hgdp1kg_1200/panel_1199.ids"
)

# Map file (IDs -> superpopulation), 2 columns: ID <tab/space> SUPERPOP
SUPERPOP_MAP = Path("/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map")

# Directory containing genetic_map_chrN_rfmix.txt
MAPDIR = Path("/home/rare/ivon/outputs/rfmix/rfmix_map")
# ==========================================


def run(cmd):
    print(">", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True)


def load_genetic_map(chr_str: str):
    gmap = MAPDIR / f"genetic_map_chr{chr_str}_rfmix.txt"
    if not gmap.exists():
        raise FileNotFoundError(f"Genetic map not found: {gmap}")

    df = pd.read_csv(
        gmap,
        sep=r"\s+",
        header=None,
        names=["pos", "cM", "id"],
        dtype={"pos": float, "cM": float, "id": str},
    ).dropna(subset=["pos", "cM"])

    df = df.sort_values("pos")
    return df["pos"].to_numpy(dtype=float), df["cM"].to_numpy(dtype=float)


def make_snp_locations_3col(vcf_path: Path, chr_str: str, out_path: Path):
    """
    Output in 3-column format:
      chrN <TAB> POS <TAB> cM
    """
    pos_arr, cm_arr = load_genetic_map(chr_str)
    vcf = VCF(str(vcf_path))

    out_lines = []
    total = 0
    kept = 0

    for v in vcf:
        total += 1
        if (not v.is_snp) or (len(v.ALT) != 1):
            continue

        pos = int(v.POS)
        cm = float(np.interp(pos, pos_arr, cm_arr))
        out_lines.append(f"{v.CHROM}\t{pos}\t{cm}\n")
        kept += 1

    if not out_lines:
        raise RuntimeError(f"No biallelic SNPs found in {vcf_path}")

    out_path.write_text("".join(out_lines))
    print(f"   snp_locations: {out_path.name} | SNPs={kept} | Total_VCF={total}")


def detect_query_id(merged_vcf: Path) -> str:
    """
    In merged files, the last sample is typically the query sample.
    """
    samples = subprocess.check_output(
        ["bcftools", "query", "-l", str(merged_vcf)],
        text=True
    ).strip().splitlines()

    if not samples:
        raise RuntimeError(f"VCF contains no samples in header: {merged_vcf}")

    return samples[-1]


def make_labels_file_for_ref(ref_vcf: Path, out_labels: Path):
    """
    Creates chrN.superpopulation_labels.txt with 2 columns:
      <ID> <SUPERPOP>
    For all samples in ref_vcf (1199).
    """
    ref_ids = subprocess.check_output(
        ["bcftools", "query", "-l", str(ref_vcf)],
        text=True
    ).strip().splitlines()

    if not ref_ids:
        raise RuntimeError(f"Reference VCF contains no samples: {ref_vcf}")

    pop = {}
    with SUPERPOP_MAP.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            pop[parts[0]] = parts[1]

    na = 0
    with out_labels.open("w") as out:
        for sid in ref_ids:
            sp = pop.get(sid, "NA")
            if sp == "NA":
                na += 1
            out.write(f"{sid}\t{sp}\n")

    print(f"   labels: {out_labels.name} | n={len(ref_ids)} | NA={na}")


def main():
    for chrn in CHRS:
        chr_str = str(chrn)
        idir = BASE / f"isec_chr{chr_str}"
        merged = idir / f"chr{chr_str}.panel1199_plus_018.common.merged.vcf.gz"

        if not merged.exists():
            print(f"chr{chr_str}: {merged.name} not found, skipping")
            continue

        print(f"\n=========== chr{chr_str} (018) ===========")
        query_id = detect_query_id(merged)
        print("   Detected query ID:", query_id)

        # ---------- REF ----------
        ref_vcf = idir / f"chr{chr_str}.ref.vcf.gz"
        run([
            "bcftools", "view",
            "-S", str(PANEL_IDS),
            "-Oz", "-o", str(ref_vcf),
            str(merged)
        ])
        run(["tabix", "-p", "vcf", str(ref_vcf)])

        # ---------- QUERY ----------
        query_vcf = idir / f"chr{chr_str}.query.vcf.gz"
        run([
            "bcftools", "view",
            "-s", query_id,
            "-Oz", "-o", str(query_vcf),
            str(merged)
        ])
        run(["tabix", "-p", "vcf", str(query_vcf)])

        # ---------- SNP locations ----------
        snp_loc = idir / f"chr{chr_str}.snp_locations.fixed"
        make_snp_locations_3col(merged, chr_str, snp_loc)

        # ---------- Labels ----------
        labels = idir / f"chr{chr_str}.superpopulation_labels.txt"
        make_labels_file_for_ref(ref_vcf, labels)

        print(f"chr{chr_str} completed: ref | query | snp_locations.fixed | labels")

    print("\nCompleted. Files generated per chromosome:")
    print("chrN.ref.vcf.gz, chrN.query.vcf.gz, chrN.snp_locations.fixed, chrN.superpopulation_labels.txt")


if __name__ == "__main__":
    main()

_______________________________________________________________________________________________________________________________#RUN_RFMIX#

#!/usr/bin/env bash
set -euo pipefail

RFMIX="$HOME/programs/rfmix/rfmix"
THREADS=8

# =========================
# 017
# =========================
BASE_017="$HOME/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/18genfased/017/merged"
OUT_017="${BASE_017}/rfmix_results"

mkdir -p "$OUT_017"

for CHR in {1..22}; do
    echo "===== RFMix 017 chr${CHR} ====="

    mkdir -p "${OUT_017}/chr${CHR}"

    "$RFMIX" \
      -f "${BASE_017}/isec_chr${CHR}/chr${CHR}.query.vcf.gz" \
      -r "${BASE_017}/isec_chr${CHR}/chr${CHR}.ref.vcf.gz" \
      -m "${BASE_017}/isec_chr${CHR}/chr${CHR}.superpopulation_labels.txt" \
      -g "${BASE_017}/isec_chr${CHR}/chr${CHR}.snp_locations" \
      --chromosome="chr${CHR}" \
      -o "${OUT_017}/chr${CHR}/rfmix_chr${CHR}_output" \
      --n-threads="$THREADS"
done

# =========================
# 018
# =========================
BASE_018="$HOME/ivon/outputs/rfmix/merged_018"
OUT_018="${BASE_018}/RFMIX"

mkdir -p "$OUT_018"

for CHR in {1..22}; do
    echo "===== RFMix 018 chr${CHR} ====="

    "$RFMIX" \
      -f "${BASE_018}/chr${CHR}.query.vcf.gz" \
      -r "${BASE_018}/chr${CHR}.ref.vcf.gz" \
      -m "${BASE_018}/superpopulation_labels.txt" \
      -g "${BASE_018}/merged_chr${CHR}.snp_locations" \
      --chromosome="chr${CHR}" \
      -o "${OUT_018}/rfmix_chr${CHR}_output" \
      --n-threads="$THREADS"
done

______________________________________________________________________________________________________________________________#PLOTS_RFMIX#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Paths ===
BASE_DIR = "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
PREFIX = "merged_18_1200_common_qc05_maf05_nomind_pruned"

# IID -> SUPERPOP map (2 columns: IID \t AFR/EUR/EAS/AMR)
MAP_SUPERPOP = "/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map"

# Ks already computed
K_LIST = [1, 2, 3, 4, 5, 6]

# Block order in the plot
BLOCK_ORDER = ["AFR", "EUR", "EAS", "AMR", "OWN", "UNK"]

# === Utilities ===
def short_pb(iid: str) -> str:
    # PB000696_017C... -> "017"
    if iid.startswith("PB000696_"):
        try:
            return iid.split("_")[1][:3]
        except Exception:
            return iid
    return iid

def load_fam_ids(fam_path: str) -> pd.DataFrame:
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None)
    fam.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    return fam[["FID", "IID"]]

def load_map(map_path: str) -> pd.DataFrame:
    m = pd.read_csv(map_path, sep=r"\s+", header=None, names=["IID", "SUPERPOP"])
    return m

def ensure_q_shape(q: pd.DataFrame, k: int):
    if q.shape[1] != k:
        raise ValueError(f"Q file has {q.shape[1]} columns but K={k}. Check the .Q file.")

def plot_one_k(df_ordered: pd.DataFrame, Q: np.ndarray, k: int, out_png: str):
    x = np.arange(df_ordered.shape[0])
    bottom = np.zeros(df_ordered.shape[0])

    plt.figure(figsize=(20, 6))

    for j in range(k):
        plt.bar(x, Q[:, j], bottom=bottom, width=1.0)
        bottom += Q[:, j]

    # Block separators
    blocks = df_ordered["BLOCK"].values
    cuts = np.where(blocks[:-1] != blocks[1:])[0]

    for c in cuts:
        plt.axvline(c + 0.5, linewidth=1)

    plt.xticks([], [])

    # Labels for PB samples
    pb_mask = df_ordered["IID"].str.startswith("PB000696_").fillna(False).values

    for i in np.where(pb_mask)[0]:
        lab = df_ordered.iloc[i]["SHORT"]

        if lab in ("017", "018"):
            plt.text(i, 1.02, lab, ha="center", va="bottom",
                     fontsize=10, fontweight="bold", rotation=90)
        else:
            plt.text(i, 1.02, lab, ha="center", va="bottom",
                     fontsize=8, rotation=90)

    # Block labels
    start = 0

    for b in BLOCK_ORDER:
        idx = np.where(df_ordered["BLOCK"].values == b)[0]

        if len(idx) == 0:
            continue

        mid = (idx[0] + idx[-1]) / 2

        plt.text(
            mid, 1.10,
            f"{b} (n={len(idx)})",
            ha="center",
            va="bottom",
            fontsize=11
        )

        start = idx[-1] + 1

    plt.ylim(0, 1.15)
    plt.ylabel("Proportion (Q)")
    plt.title(
        f"ADMIXTURE (supervised/unsupervised) – {PREFIX} – K={k}\n"
        "Ordered by superpopulation + 18 genomes (017 and 018 highlighted)"
    )

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

# === Main ===
def main():
    os.chdir(BASE_DIR)

    fam_path = f"{PREFIX}.fam"
    fam = load_fam_ids(fam_path)

    mp = load_map(MAP_SUPERPOP)

    # Merge while preserving .fam order
    df = fam.merge(mp, on="IID", how="left")

    # Block classification
    df["BLOCK"] = df["SUPERPOP"].fillna("UNK")

    # Own genomes
    mask_own = df["IID"].str.startswith("PB000696_").fillna(False)
    df.loc[mask_own, "BLOCK"] = "OWN"

    # Short labels for PB samples
    df["SHORT"] = df["IID"].apply(short_pb)

    # Sort by block
    df["BLOCK"] = pd.Categorical(df["BLOCK"], categories=BLOCK_ORDER, ordered=True)
    df_sorted = df.sort_values(["BLOCK", "IID"]).reset_index(drop=True)

    # Reorder Q matrix
    idx_map = pd.Series(np.arange(df.shape[0]), index=df["IID"]).to_dict()
    original_positions = df_sorted["IID"].map(idx_map).values

    print("Counts per block:")
    print(df_sorted["BLOCK"].value_counts(dropna=False))

    for k in K_LIST:
        qfile = f"{PREFIX}.{k}.Q"

        if not os.path.exists(qfile) or os.path.getsize(qfile) == 0:
            print(f"[WARNING] File not found: {qfile} (skipped)")
            continue

        q = pd.read_csv(qfile, sep=r"\s+", header=None)
        ensure_q_shape(q, k)

        Q = q.values.astype(float)
        Q_sorted = Q[original_positions, :]

        out_png = os.path.join(
            BASE_DIR,
            f"{PREFIX}.K{k}.Q.superpop_blocks.png"
        )

        plot_one_k(df_sorted, Q_sorted, k, out_png)

        print(f"[OK] Saved: {out_png}")

if __name__ == "__main__":
    main()









