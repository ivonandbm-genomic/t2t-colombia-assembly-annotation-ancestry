#scripts/
├── 01_filter_deepvariant_snps_autosomes.sh
├── 02_phase_variants_whatshap.sh
├── 03_split_phased_vcfs_by_chromosome.sh
├── 04_build_gnomad_reference_panel.py
├── 05_extract_reference_panel_vcfs.py
├── 06_plink_merge_qc_ldpruning_for_pca_admixture.sh
├── 07_run_pca.sh
├── 08_run_admixture.py
├── 09_plot_pca_admixture.py
├── 10_rfmix_intersect_harmonize_merge_panel_query.sh
├── 11_prepare_rfmix_inputs_query_ref_labels_maps.py
├── 12_run_rfmix_all_samples.sh
└── 13_plot_local_ancestry_rfmix.py#

#1_Filtrar snps_ bialelicos_autosomas_pass_de los genomas#

#!/usr/bin/env python3

"""
Script en Python para procesar múltiples VCFs generados con DeepVariant (PacBio HiFi):
- Filtra solo SNPs bialélicos con filtro PASS
- Aplica filtros de calidad
- Normaliza usando una referencia fasta
- Asigna IDs únicos a cada variante
- Indexa los archivos VCF
"""

import os
import subprocess
from pathlib import Path
from glob import glob

# ──────────────────────────────
# CONFIGURACIÓN DEL USUARIO
# ──────────────────────────────
VCF_DIR = Path("/mnt/diskrare/arlenb/08/small_variants/hg38")
BAM_DIR = Path("/mnt/diskrare/arlenb/08/aligned_reads/hg38")
OUT_DIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil")
REF = "/home/rare/ivon/data/hg38.fa"  # Ruta al genoma de referencia
THREADS = 32


# Crear carpeta de salida si no existe
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Obtener lista de archivos VCF
vcfs = sorted(VCF_DIR.glob("*.vcf.gz"))
print(f"🔍 {len(vcfs)} VCFs encontrados para procesar.")

# Procesamiento por archivo
for vcf in vcfs:
    basename = vcf.stem.replace("_aligned", "")
    print(f"\n▶ Procesando {basename}...")

    # Paso 1: Filtrar PASS, solo SNPs bialélicos
    step1 = OUT_DIR / f"{basename}.pass.snps.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-f", "PASS", "-v", "snps",
        "-Oz", "-o", str(step1), str(vcf)
    ], check=True)

    # Paso 2: Filtro de calidad (QUAL >= 20, DP >= 12, GQ >= 20, no missing)
    step2 = OUT_DIR / f"{basename}.filtered.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-i", 'QUAL >= 20 & FORMAT/DP >= 12 & FORMAT/GQ >= 20 & GT!="mis"',
        "-Oz", "-o", str(step2), str(step1)
    ], check=True)

    # Paso 3: Normalizar usando la referencia
    step3 = OUT_DIR / f"{basename}.norm.vcf.gz"
    subprocess.run([
        "bcftools", "norm",
        "-f", REF,
        "-Oz", "-o", str(step3), str(step2)
    ], check=True)

    # Paso 4: Anotar ID único para cada SNP
    step4 = OUT_DIR / f"{basename}.norm.id.vcf.gz"
    subprocess.run([
        "bcftools", "annotate",
        "--set-id", "+'%CHROM\_%POS\_%REF\_%ALT'",
        "-Oz", "-o", str(step4), str(step3)
    ], check=True)

    # Paso 5: Indexar
    subprocess.run(["bcftools", "index", str(step4)], check=True)
print("\n🎉 Procesamiento completado para todos los VCFs.")

___________________________________________________________________________________________________________________
#2_PANEL REFERENCIA PASO 1  GENERACION DE panel_5superpoblaciones.tsv#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sys

BASE = "/mnt/diskrare/ivonb/refamerindios"

PANEL_TSV = os.path.join(BASE, "panel_5superpoblaciones.tsv")
PANEL_IDS = os.path.join(BASE, "panel_5super.ids")

VCF_BASE = BASE
OUT_BASE = os.path.join(BASE, "vcfs_panel")

BCFTOOLS = "bcftools"          # o ruta absoluta si hace falta
THREADS = 32                   # ⬅️ aquí defines los hilos


def run_cmd(cmd, **kwargs):
    print("  ➤ Ejecutando:", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)


def main():
    os.makedirs(OUT_BASE, exist_ok=True)

    if not os.path.isfile(PANEL_TSV):
        print(f"❌ No encontré {PANEL_TSV}", file=sys.stderr)
        sys.exit(1)

    print("📄 Generando lista de IDs del panel...")
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

    print(f"   -> {len(ids)} IDs escritos en {PANEL_IDS}")

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

        print(f"\n🧬 Chr{chr_str}:")
        print(f"   Entrada: {vcf_in}")
        print(f"   Salida : {out_vcf}")

        if not os.path.isfile(vcf_in):
            print(f"⚠️ No encontré {vcf_in}, lo salto.")
            continue

        # ➜ bcftools view con 32 hilos
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

        # ➜ bcftools index con 32 hilos
        cmd_index = [
            BCFTOOLS,
            "index",
            "--threads", str(THREADS),
            "-t",
            out_vcf,
        ]
        run_cmd(cmd_index)

    print(f"\n✅ VCFs filtrados por panel en: {OUT_BASE}")


if __name__ == "__main__":
    main()

#PANEL DE REFERENCIA PASO 2 SELECCION DE 300 POR SUPERPOBLACION#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Construye un panel de referencia de 5 superpoblaciones
a partir de gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz

Superpoblaciones usadas: AFR, AMR, EUR, EAS, SAS
(AMI se agrupa dentro de AMR).

Salida:
- panel_5superpoblaciones.tsv : tabla con sample, project, population, superpop
- panel_AFR.ids, panel_AMR.ids, ... : listas de IDs por superpoblación
"""

import os
import sys
import json
import argparse
import pandas as pd
from pandas import json_normalize

# ─────────────────────────────────────────────
# Parámetros por defecto (ajusta a tu gusto)
# ─────────────────────────────────────────────
DEFAULT_META = "/mnt/diskrare/ivonb/refamerindios/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz"
DEFAULT_OUTDIR = "/mnt/diskrare/ivonb/refamerindios"
N_PER_SUPERPOP = 300   # máximo individuos por superpoblación
SEED = 12345

# Columnas lógicas que usará el script
SAMPLE_COL   = "s"              # ID de muestra en el meta de gnomAD
PROJECT_COL  = "project"        # de hgdp_tgp_meta
POP_COL      = "population"     # de hgdp_tgp_meta
SUPERPOP_COL = "genetic_region" # de hgdp_tgp_meta (EUR, AFR, AMR, EAS, SAS, AMI, MID, ...)

SUPERPOPS_5 = ["AFR", "AMR", "EUR", "EAS", "SAS"]


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_cols(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(
            f"❌ Faltan columnas en el meta: {missing}\n"
            f"Columnas disponibles: {list(df.columns)}\n"
            f"👉 Revisa el script o el archivo de entrada."
        )


def map_to_5_superpops(genetic_region):
    """
    Mapea genetic_region a 5 grupos: AFR, AMR, EUR, EAS, SAS.
    - AMI se agrupa en AMR (Amerindio/Latino)
    - Otros (MID, OTH, etc.) se devuelven como None (se excluyen)
    """
    if pd.isna(genetic_region):
        return None
    gr = str(genetic_region).upper()
    if gr in SUPERPOPS_5:
        return gr
    if gr == "AMI":   # Indigenous American
        return "AMR"
    return None


def parse_json_field(x):
    """
    Parsea un campo JSON de gnomAD (como hgdp_tgp_meta).
    Maneja NA / nulls de forma segura.
    """
    if pd.isna(x):
        return {}
    s = str(x).strip()
    if s == "" or s.upper() == "NA":
        return {}
    try:
        return json.loads(s)
    except Exception:
        # Si por alguna razón viene algo raro, lo registramos y devolvemos {}
        eprint(f"⚠️ No pude parsear JSON: {s[:80]}...")
        return {}


# ─────────────────────────────────────────────
# Lógica principal
# ─────────────────────────────────────────────
def build_panel(meta_path: str, outdir: str):
    print(f"📂 Leyendo meta: {meta_path}")
    df_raw = pd.read_csv(
        meta_path,
        sep="\t",
        compression="gzip",
        low_memory=False
    )
    print("✅ Columnas crudas del meta:", list(df_raw.columns))

    # Verificación básica
    if "hgdp_tgp_meta" not in df_raw.columns:
        raise SystemExit("❌ No encuentro la columna 'hgdp_tgp_meta' en el meta.")

    if SAMPLE_COL not in df_raw.columns:
        raise SystemExit(f"❌ No encuentro la columna de ID '{SAMPLE_COL}' en el meta.")

    # ─────────────────────────────────────────
    # 1) Desempaquetar hgdp_tgp_meta
    # ─────────────────────────────────────────
    meta_parsed = df_raw["hgdp_tgp_meta"].apply(parse_json_field)
    meta_flat = json_normalize(meta_parsed)

    print("🔎 Subcolumnas en hgdp_tgp_meta:", meta_flat.columns.tolist())

    # Construimos un DataFrame con:
    # - ID de muestra (s)
    # - campos de hgdp_tgp_meta (project, population, genetic_region, etc.)
    cols_from_raw = [SAMPLE_COL]
    # Si quieres añadir flags de calidad:
    for qc_col in ["gnomad_high_quality", "high_quality"]:
        if qc_col in df_raw.columns:
            cols_from_raw.append(qc_col)

    df = pd.concat([df_raw[cols_from_raw], meta_flat], axis=1)

    # Verificar que estén las columnas de interés
    check_cols(df, [SAMPLE_COL, PROJECT_COL, POP_COL, SUPERPOP_COL])

    print("✅ Columnas disponibles para el panel:",
          [SAMPLE_COL, PROJECT_COL, POP_COL, SUPERPOP_COL])

    # ─────────────────────────────────────────
    # 2) Filtros opcionales de calidad / proyecto
    # ─────────────────────────────────────────
    # Filtrar por project: típicamente "1000 Genomes" y/o "HGDP"
    proyectos_aceptados = ["1000 Genomes", "HGDP"]
    df = df[df[PROJECT_COL].isin(proyectos_aceptados)]
    print(f"📉 Tras filtrar por proyecto {proyectos_aceptados}: {len(df)} muestras")

    # Filtrar por calidad, si las columnas existen
    if "gnomad_high_quality" in df.columns:
        antes = len(df)
        df = df[df["gnomad_high_quality"] == True]
        print(f"📉 Tras gnomad_high_quality: {antes} → {len(df)}")

    if "high_quality" in df.columns:
        antes = len(df)
        df = df[df["high_quality"] == True]
        print(f"📉 Tras high_quality: {antes} → {len(df)}")

    # ─────────────────────────────────────────
    # 3) Mapeo a 5 superpoblaciones
    # ─────────────────────────────────────────
    df["superpop_5"] = df[SUPERPOP_COL].apply(map_to_5_superpops)
    antes = len(df)
    df = df[df["superpop_5"].notna()]
    print(f"📉 Tras mapear a 5 superpoblaciones: {antes} → {len(df)}")

    # Nos quedamos solo con las 5 superpops que nos interesan
    df = df[df["superpop_5"].isin(SUPERPOPS_5)]
    print("📊 Recuento por superpoblación (antes de muestrear):")
    print(df["superpop_5"].value_counts())

    # ─────────────────────────────────────────
    # 4) Muestreo equilibrado por superpoblación
    # ─────────────────────────────────────────
    os.makedirs(outdir, exist_ok=True)

    panel_rows = []
    for sp in SUPERPOPS_5:
        sub = df[df["superpop_5"] == sp]
        n_disponibles = len(sub)
        if n_disponibles == 0:
            eprint(f"⚠️ No hay muestras para superpoblación {sp}, se omite.")
            continue
        n_tomar = min(N_PER_SUPERPOP, n_disponibles)
        sub_sel = sub.sample(n_tomar, random_state=SEED)

        # Guardar lista de IDs para cada superpoblación
        out_ids = os.path.join(outdir, f"panel_{sp}.ids")
        sub_sel[SAMPLE_COL].to_csv(out_ids, index=False, header=False)
        print(f"💾 Guardado {n_tomar} IDs para {sp} en {out_ids}")

        # Añadir al panel combinado
        for _, row in sub_sel.iterrows():
            panel_rows.append({
                "sample": row[SAMPLE_COL],
                "project": row[PROJECT_COL],
                "population": row[POP_COL],
                "superpop": row["superpop_5"],
            })

    if not panel_rows:
        raise SystemExit("❌ No se seleccionó ninguna muestra para el panel final.")

    panel_df = pd.DataFrame(panel_rows)
    out_panel = os.path.join(outdir, "panel_5superpoblaciones.tsv")
    panel_df.to_csv(out_panel, sep="\t", index=False)
    print(f"✅ Panel combinado guardado en: {out_panel}")
    print("📊 Recuento final en el panel:")
    print(panel_df["superpop"].value_counts())


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Construye un panel de referencia de 5 superpoblaciones "
                    "a partir del meta HGDP+1KG de gnomAD."
    )
    p.add_argument(
        "--meta",
        default=DEFAULT_META,
        help=f"Ruta al archivo meta (.tsv.bgz) [default: {DEFAULT_META}]"
    )
    p.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help=f"Carpeta de salida [default: {DEFAULT_OUTDIR}]"
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_panel(args.meta, args.outdir)
#ALL CROM#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sys

BASE = "/mnt/diskrare/ivonb/refamerindios"

PANEL_TSV = os.path.join(BASE, "panel_5superpoblaciones.tsv")
PANEL_IDS = os.path.join(BASE, "panel_5super.ids")

VCF_BASE = BASE
OUT_BASE = os.path.join(BASE, "vcfs_panel")

BCFTOOLS = "bcftools"          # o ruta absoluta si hace falta
THREADS = 32                   # ⬅️ aquí defines los hilos


def run_cmd(cmd, **kwargs):
    print("  ➤ Ejecutando:", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)


def main():
    os.makedirs(OUT_BASE, exist_ok=True)

    if not os.path.isfile(PANEL_TSV):
        print(f"❌ No encontré {PANEL_TSV}", file=sys.stderr)
        sys.exit(1)

    print("📄 Generando lista de IDs del panel...")
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

    print(f"   -> {len(ids)} IDs escritos en {PANEL_IDS}")

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

        print(f"\n🧬 Chr{chr_str}:")
        print(f"   Entrada: {vcf_in}")
        print(f"   Salida : {out_vcf}")

        if not os.path.isfile(vcf_in):
            print(f"⚠️ No encontré {vcf_in}, lo salto.")
            continue

        # ➜ bcftools view con 32 hilos
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

        # ➜ bcftools index con 32 hilos
        cmd_index = [
            BCFTOOLS,
            "index",
            "--threads", str(THREADS),
            "-t",
            out_vcf,
        ]
        run_cmd(cmd_index)

    print(f"\n✅ VCFs filtrados por panel en: {OUT_BASE}")


if __name__ == "__main__":
    main()

______________________________________________________________________________
______________________________________________________________________________
#PLINK / PCA / ADMIXTURE
panel + 18 genomas
↓
SNPs comunes
↓
merge PLINK
↓
QC
↓
MAF
↓
LD pruning
↓
PCA / ADMIXTURE#
__________________________________________________________________________________________________________________#PLINK#

#PLINK#

#!/usr/bin/env bash
set -euo pipefail

PLINK="/home/rare/programs/Plink/plink"

WORKDIR="/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
cd "$WORKDIR"

# 1. Intersección de SNPs comunes entre panel y 18 genomas
cut -f2 panel_5super_autosomes_fixid.bim | sort -u > ids.panel.txt
cut -f2 18genomes_auto_fixid.bim        | sort -u > ids.18g.txt
comm -12 ids.panel.txt ids.18g.txt > ids.common.txt

# 2. Extraer SNPs comunes
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

# 3. Merge panel + 18 genomas
$PLINK \
  --bfile panel_5super_autosomes_common \
  --bmerge 18genomes_auto_common.bed 18genomes_auto_common.bim 18genomes_auto_common.fam \
  --make-bed \
  --out merged_18_1200_common_raw

# 4. QC básico
$PLINK \
  --bfile merged_18_1200_common_raw \
  --geno 0.05 \
  --make-bed \
  --out merged_18_1200_common_qc05

# 5. Filtro MAF 0.05
$PLINK \
  --bfile merged_18_1200_common_qc05 \
  --maf 0.05 \
  --make-bed \
  --out merged_18_1200_common_qc05_maf05_nomind

# 6. LD pruning
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

echo "🔹 Ejecutando PCA"

$PLINK \
  --bfile merged_18_1200_common_qc05_maf05_nomind_pruned \
  --pca 20 \
  --out PCA_m05

echo "✅ PCA listo"


___________________________________________________________________________________________________________#PLOT_PCA#
#PLOT_PCA#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# === 0) Rutas y Configuración ===
BASE_DIR = "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
OUT_DIR = Path("/home/rare/ivon/figpaper")
OUT_DIR.mkdir(parents=True, exist_ok=True)

os.chdir(BASE_DIR)

# === 1) Cargar eigenvec del PCA ===
evec = pd.read_csv("PCA_m05.eigenvec", sep=r"\s+", header=None)
evec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, 21)]

# === 2) Cargar mapa de superpoblaciones ===
pop = pd.read_csv(
    "/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map",
    sep=r"\s+",
    header=None,
    names=["IID", "POP"]
)

# Unir PCA + poblaciones
df = evec.merge(pop, on="IID", how="left")
df["POP"] = df["POP"].fillna("UNK")

# === 3) Identificar muestras 017/018 ===
mask_own = df["IID"].astype(str).str.startswith("PB000696_")
df.loc[mask_own, "SHORT"] = (
    df.loc[mask_own, "IID"].str.split("_").str[1].str[:3]
)

mask_label = (mask_own & df["SHORT"].isin(["017", "018"])) | df["IID"].astype(str).isin(["017", "018"])

# === 4) Paleta Maestra Unificada ===
color_map = {
    "AMR": "#1F77B4",  # Azul
    "EUR": "#FF7F0E",  # Naranja
    "AFR": "#2CA02C",  # Verde
    "EAS": "#D62728",  # Rojo
    "SAS": "#9467BD",  # Púrpura
    "UNK": "#D3D3D3",  # Gris
}
df["color"] = df["POP"].map(color_map).fillna("#D3D3D3")

# === 5) Dibujar ===
# Ajustamos el tamaño de la figura para que la resolución de 400 dpi sea efectiva
fig, ax = plt.subplots(figsize=(10, 8))

legend_pops = [p for p in ["AFR", "AMR", "EUR", "EAS", "SAS"] if p in set(df["POP"])]

# Dibujar puntos de referencia
for pop_label in legend_pops:
    sub = df[df["POP"] == pop_label]
    if not sub.empty:
        ax.scatter(
            sub["PC1"], sub["PC2"],
            s=30, # Puntos ligeramente más grandes para 400 dpi
            c=sub["color"],
            label=pop_label,
            alpha=0.7,
            edgecolors='white',
            linewidths=0.2
        )

# Dibujar UNK (referencia tenue)
sub_unk = df[df["POP"] == "UNK"]
if not sub_unk.empty:
    ax.scatter(sub_unk["PC1"], sub_unk["PC2"], s=15, c="#D3D3D3", alpha=0.3, zorder=1)

# Etiquetas para 017 y 018
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

# Ejes y Limpieza Estética
ax.set_xlabel("PC1", fontsize=12, fontweight='bold')
ax.set_ylabel("PC2", fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Leyenda
ax.legend(
    title="Superpopulation",
    title_fontsize=10,
    fontsize=9,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False
)

plt.tight_layout()

# === 6) Guardado a 400 DPI ===
out_name = "PCA_PC1_PC2_400DPI_017_018"
plt.savefig(OUT_DIR / f"{out_name}.png", dpi=400, bbox_inches="tight")
plt.savefig(OUT_DIR / f"{out_name}.pdf", bbox_inches="tight")

plt.close()

print(f"✅ Archivos generados a 400 DPI en: {OUT_DIR}")
________________________________________________________________________________________________________________________#ADMIXURE#
# ADMIXURE#
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
        sys.exit(f"❌ No encuentro ADMIXTURE en: {ADMIXTURE}")

    for ext in ("bed", "bim", "fam"):
        f = WORKDIR / f"{PREFIX}.{ext}"
        if not f.exists() or f.stat().st_size == 0:
            sys.exit(f"❌ Falta o está vacío: {f}")

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
    print(f"🚀 Corriendo ADMIXTURE K={k}")
    print("CMD:", " ".join(cmd))
    print("LOG:", log_file.name)
    print("======================================")

    # Ejecuta y hace "tee" (imprime y guarda)
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
        sys.exit(f"❌ ADMIXTURE falló en K={k} (return code={rc}). Revisa {log_file}")

    return log_file

def extract_cv_error(log_file: Path) -> str | None:
    # Busca línea tipo: "CV error (K=4): 0.38712"
    pat = re.compile(r"CV error\s*\(K=\d+\)\s*:\s*([0-9.eE+-]+)")
    with open(log_file) as f:
        for line in f:
            m = pat.search(line)
            if m:
                return m.group(1)
    return None

def main():
    check_inputs()
    print("📂 WORKDIR:", WORKDIR)
    print("🧬 PREFIX :", PREFIX)
    print("▶ ADMIXTURE:", ADMIXTURE)
    print(f"▶ Params: --cv={CV_FOLDS} -j{THREADS} --seed={SEED}  |  K={K_MIN}..{K_MAX}")

    results = []

    for k in range(K_MIN, K_MAX + 1):
        logf = run_one_k(k)
        cv = extract_cv_error(logf)
        if cv is None:
            print(f"⚠️ No se encontró 'CV error' en {logf.name}")
        else:
            results.append((k, cv))

    # Guardar tabla resumen
    out_tsv = WORKDIR / "cv_error.tsv"
    with open(out_tsv, "w") as out:
        out.write("K\tCV_error\n")
        for k, cv in results:
            out.write(f"{k}\t{cv}\n")

    print("\n✅ Listo.")
    print("✅ Logs: cv_K1.log ... cv_K6.log (en la misma carpeta)")
    print("✅ Resumen:", out_tsv)

    # Mostrar en pantalla
    if results:
        print("\nCV errors:")
        for k, cv in results:
            print(f"K={k}\tCV={cv}")

if __name__ == "__main__":
    main()
__________________________________________________________________________________________________________#PLOT_ADMIXURE#


 #PLOT ADMIXURE#
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
K_MAIN = 4  # barplot principal + pies

# Bloques iguales
BLOCKS_ORDER = ["AFR", "EUR", "EAS", "AMR", "COLM", "COLF"]
BLOCK_WIDTH        = 200
GAP_BETWEEN_BLOCKS = 60

# Project blocks
KEEP_COLM = {"017"}
KEEP_COLF = {"018"}

# Si quieres forzar manualmente los labels (debe tener largo=K_MAIN), descomenta:
# COMP_LABELS_K4 = ["EUR", "AMR", "EAS", "AFR"]
COMP_LABELS_K4 = None  # si None, intentará inferir labels por bloque (AFR/EUR/EAS/AMR)

os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# HELPERS
# =========================
def check_file(path: str):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Falta o está vacío: {path}")

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
    Construye df con IID, SHORT, GROUP y layout de bloques:
    - x_positions (borde izquierdo)
    - widths por fila (para llenar bloques iguales)
    - block_centers, block_ranges, used_blocks (para etiquetas/separadores)
    """
    df = fam.merge(pop, on="IID", how="left")
    df["SHORT"] = df["IID"].apply(short_id)

    # GROUP: por defecto SUPERPOP; si es Project y SHORT coincide -> COLM/COLF
    df["GROUP"] = df["SUPERPOP"]
    df.loc[df["IID"].apply(is_project) & df["SHORT"].isin(KEEP_COLM), "GROUP"] = "COLM"
    df.loc[df["IID"].apply(is_project) & df["SHORT"].isin(KEEP_COLF), "GROUP"] = "COLF"

    # Filtrar
    df = df[df["GROUP"].notna()].copy()
    df = df[df["GROUP"] != "UNK"].copy()
    df = df[df["GROUP"].isin(BLOCKS_ORDER)].copy()

    # Orden por bloques; dentro de bloque por SHORT
    df["BLOCK_ORDER"] = pd.Categorical(df["GROUP"], categories=BLOCKS_ORDER, ordered=True)
    df.sort_values(["BLOCK_ORDER", "SHORT"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Conteo por bloque
    sample_counts = {b: int((df["GROUP"] == b).sum()) for b in BLOCKS_ORDER}
    sample_counts = {b: n for b, n in sample_counts.items() if n > 0}

    # Ancho por bloque para llenar el bloque completo
    bar_widths_block = {b: (BLOCK_WIDTH / sample_counts[b]) for b in sample_counts}

    # Posiciones x por individuo (borde izquierdo)
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
    # etiquetas arriba
    for c, lab in zip(block_centers, used_blocks):
        ax.text(
            c, 1.02, lab,
            ha="center", va="bottom",
            transform=ax.get_xaxis_transform(),
            fontsize=11, fontweight="bold"
        )
    # separadores
    for j in range(1, len(block_ranges)):
        sep_pos = block_ranges[j-1][1] + (GAP_BETWEEN_BLOCKS / 2)
        ax.axvline(sep_pos, color="gray", linestyle="--", linewidth=1)

def infer_component_labels_from_blocks(df: pd.DataFrame, Qk: np.ndarray, ref_blocks=("AFR","EUR","EAS","AMR")):
    """
    Intenta asignar cada componente (columna de Qk) a AFR/EUR/EAS/AMR
    usando el mayor promedio dentro de cada bloque de referencia.
    Devuelve labels del tamaño K (ej: ["AFR","EUR","EAS","AMR"] en algún orden).

    Nota: esto es un heurístico razonable para figuras “supervisadas por bloques”
    cuando tus bloques AFR/EUR/EAS/AMR son referencias.
    """
    K = Qk.shape[1]
    block_means = {}
    for b in ref_blocks:
        idx = df.index[df["GROUP"] == b].to_list()
        if len(idx) == 0:
            continue
        block_means[b] = Qk[idx, :].mean(axis=0)

    # si faltan bloques, cae a C1..CK
    if len(block_means) < 2:
        return [f"C{j+1}" for j in range(K)]

    # greedy assignment: componente -> bloque con mayor mean, evitando repetir bloque si es posible
    labels = [None] * K
    used = set()

    # ordenar componentes por "claridad" (diferencia entre top1 y top2)
    clarity = []
    for j in range(K):
        scores = [(b, block_means[b][j]) for b in block_means]
        scores.sort(key=lambda x: x[1], reverse=True)
        top1 = scores[0][1]
        top2 = scores[1][1] if len(scores) > 1 else 0.0
        clarity.append((j, top1 - top2, scores))
    clarity.sort(key=lambda x: x[1], reverse=True)

    for j, _, scores in clarity:
        # elige el mejor bloque no usado; si todos usados, el mejor bloque
        chosen = None
        for b, _v in scores:
            if b not in used:
                chosen = b
                break
        if chosen is None:
            chosen = scores[0][0]
        labels[j] = chosen
        used.add(chosen)

    # si hay componentes extra (K>4), marca “UNK1..”
    # (porque ADMIXTURE puede separar subcomponentes dentro de un bloque)
    if K > len(ref_blocks):
        # contamos cuántas veces se repite un label
        counts = {}
        for lab in labels:
            counts[lab] = counts.get(lab, 0) + 1
        # renombrar repeticiones: AFR-1, AFR-2...
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
    Paleta consistente para barras y tortas.
    Usa tab10/tab20 sin especificar colores “a mano” pero sí consistente.
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

        # (sin título) pero con etiquetas de bloques arriba (si quieres quitarlas aquí, dímelo)
        add_block_labels_and_separators(ax, block_centers, block_ranges, used_blocks)
        ax.set_xticks([])

    # SIN título global
    plt.tight_layout()
    plt.savefig(out_png, dpi=400, bbox_inches="tight")
    plt.close()
    print("✅ Guardado:", out_png)

def plot_barplot_oneK_K4(df, x_positions, widths, Qk, K, out_png, out_pdf, block_centers, block_ranges, used_blocks, comp_labels):
    """
    Barplot K=4 (principal) SIN título, pero conservando etiquetas de bloques.
    También genera una leyenda de componentes (ancestrías) abajo si comp_labels viene definido.
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

    # ✅ mantener etiquetas de bloques (AFR/EUR/EAS/AMR/COLM/COLF)
    add_block_labels_and_separators(ax, block_centers, block_ranges, used_blocks)

    # Leyenda (sin título)
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
    print("✅ Guardado:", out_png)
    print("✅ Guardado:", out_pdf)

def plot_pies_K4_017_018(df, Qk, K, out_png, comp_labels):
    """
    2 tortas (017=COLM, 018=COLF) SIN título.
    - Colores consistentes con barplot (mismo orden de componentes)
    - Leyenda con labels
    - Porcentajes en cada porción
    """
    df_colm = df[df["GROUP"] == "COLM"]
    df_colf = df[df["GROUP"] == "COLF"]

    if df_colm.shape[0] == 0 or df_colf.shape[0] == 0:
        print("⚠️ No encontré COLM y/o COLF para generar las tortas. Revisa IIDs en .fam y KEEP_COLM/KEEP_COLF.")
        return

    i_colm = int(df_colm.index[0])
    i_colf = int(df_colf.index[0])

    if comp_labels is None or len(comp_labels) != Qk.shape[1]:
        comp_labels = [f"C{j+1}" for j in range(Qk.shape[1])]

    colors = get_component_colors(Qk.shape[1])

    def autopct_fmt(pct):
        # muestra solo si >=1%
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

    # leyenda única abajo (sin título)
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
    print("✅ Guardado:", out_png)

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

    # ---------- 1) Multipanel K=1..6 (SIN título) ----------
    out_multi = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MIN}-K{K_MAX}.equalblocks_COLM_COLF.png")
    plot_multipanel_Ks(
        df, x_positions, widths, fam,
        block_centers, block_ranges, used_blocks,
        out_multi, K_MIN, K_MAX
    )

    # ---------- 2) K=4 barplot principal (SIN título, con labels de bloques) ----------
    q_main = os.path.join(BASE_DIR, f"{PREFIX}.{K_MAIN}.Q")
    check_file(q_main)

    Q = load_Q(q_main)
    Qk_main = align_Q_to_df(Q, fam, df)

    # Labels de componentes para K=4:
    if COMP_LABELS_K4 is not None:
        if len(COMP_LABELS_K4) != Qk_main.shape[1]:
            sys.exit(f"❌ COMP_LABELS_K4 tiene largo {len(COMP_LABELS_K4)} pero K={Qk_main.shape[1]}")
        comp_labels = COMP_LABELS_K4
    else:
        # inferencia automática por bloques AFR/EUR/EAS/AMR
        comp_labels = infer_component_labels_from_blocks(df, Qk_main, ref_blocks=("AFR","EUR","EAS","AMR"))

    out_k4_png = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MAIN}.barplot_equalblocks_COLM_COLF.png")
    out_k4_pdf = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MAIN}.barplot_equalblocks_COLM_COLF.pdf")

    plot_barplot_oneK_K4(
        df, x_positions, widths, Qk_main, K_MAIN,
        out_k4_png, out_k4_pdf,
        block_centers, block_ranges, used_blocks,
        comp_labels=comp_labels
    )

    # ---------- 3) Tortas K=4 (017 y 018) con % + leyenda + colores consistentes ----------
    out_pies = os.path.join(OUT_DIR, f"{PREFIX}.K{K_MAIN}.pies_017_COLM_018_COLF.png")
    plot_pies_K4_017_018(df, Qk_main, K_MAIN, out_pies, comp_labels=comp_labels)

    print("\n✅ Todo listo. Revisa:", OUT_DIR)
    print("   - multipanel K1..K6 (sin título)")
    print("   - barplot K4 (sin título, con labels de bloques)")
    print("   - pies K4 (con % y leyenda, colores consistentes)")

if __name__ == "__main__":
    main()

_________________________________________________________________________________________________________________________
#RFMIX#
_________________________________________________________________________________________________________________________
VCFs faseados por cromosoma
↓
10_rfmix_intersect_harmonize_merge_panel_query.sh
↓
merged panel1199 + muestra
↓
11_prepare_rfmix_inputs_query_ref_labels_maps.py
↓
chrN.query.vcf.gz
chrN.ref.vcf.gz
chrN.superpopulation_labels.txt
chrN.snp_locations
↓
12_run_rfmix_all_samples.sh

___________________________________________________________________________________________________________________________#FASEO#

#Fasear genomas  #
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from pathlib import Path

# ───────────────────────────────
# CONFIGURACIÓN PARA TUS 18 GENOMAS
# ───────────────────────────────
VCF_DIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes")
BAM_DIR = Path("/mnt/diskrare/arlenb/08/aligned_reads/hg38")
OUT_DIR = Path("/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/18genfased")
REF = "/mnt/diskrare/arlenb/reference/hg38.fasta"  # ruta que indicaste
WHATSHAP = "/home/rare/.local/bin/whatshap"

# Crear carpeta de salida
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Log para problemas
log = open("log_faseo_18.txt", "w")

# Tus 18 muestras exactas (de tu ls)
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

print("Iniciando faseo de los 18 genomas con WhatsHap...\n")

for sample in SAMPLES:
    vcf_path = VCF_DIR / f"{sample}.autosomes.vcf.gz"
    bam_path = BAM_DIR / f"{sample}.bam"
    bai_path = BAM_DIR / f"{sample}.bam.bai"

    print(f"🔄 Procesando: {sample}")

    if not vcf_path.exists():
        msg = f"❌ VCF faltante: {vcf_path}"
        print(msg)
        log.write(msg + "\n")
        continue

    if not bam_path.exists():
        msg = f"❌ BAM faltante: {bam_path} (verifica si existe .bam.bai también)"
        print(msg)
        log.write(msg + "\n")
        continue

    if not bai_path.exists():
        print("⚠️  No hay índice .bai → creando uno rápido...")
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
        print("   Faseo OK")

        # Indexar
        subprocess.run(["tabix", "-p", "vcf", str(phased_vcf)], check=True)

        # Opcional: filtrar solo variantes completamente faseadas (ambos alelos con |)
        final_vcf = OUT_DIR / f"{sample}.autosomes.fully_phased.vcf.gz"
        subprocess.run([
            "bcftools", "view",
            "-i", 'GT~"|"',
            "-Oz", "-o", str(final_vcf),
            str(phased_vcf)
        ], check=True)
        subprocess.run(["tabix", "-p", "vcf", str(final_vcf)], check=True)

        print(f"✅ Completado: {sample} → {final_vcf.name}\n")
    except subprocess.CalledProcessError as e:
        msg = f"❌ Error en {sample}: {e}"
        print(msg)
        log.write(msg + "\n")

log.close()
print("🎉 ¡Todo procesado!")
print(f"VCFs faseados en: {OUT_DIR}")
print("Revisa log_faseo_18.txt para cualquier problema.")
___________________________________________________________________________________________________________________________________#SPLIT#

#SPLIT 18 GENOMAS FASEADOS POR CROMOSOMA#

#!/bin/bash
set -euo pipefail

# Directorio con los fully_phased
PHASED_DIR="18genfased"
SPLIT_DIR="18gen_by_chr"

mkdir -p $SPLIT_DIR
cd $PHASED_DIR

echo "Iniciando split por cromosoma de los 18 genomas..."

# Mapeo: nombre largo → código corto (001 a 018)
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
        echo "Advertencia: No reconozco $base → salto"
        continue
    fi
    
    sample_dir="../$SPLIT_DIR/$code"
    mkdir -p "$sample_dir"
    
    echo "Dividiendo genoma $code ($base)..."
    
    # Cromosomas 1-22
    for chr in {1..22}; do
        out="${sample_dir}/chr${chr}.vcf.gz"
        bcftools view -Oz -o "$out" "$full_vcf" $chr && tabix -p vcf "$out" &
    done
    
    # X e Y (si existen en el VCF, no falla si no)
    bcftools view -Oz -o "${sample_dir}/chrX.vcf.gz" "$full_vcf" X 2>/dev/null && tabix -p vcf "${sample_dir}/chrX.vcf.gz" || true &
    bcftools view -Oz -o "${sample_dir}/chrY.vcf.gz" "$full_vcf" Y 2>/dev/null && tabix -p vcf "${sample_dir}/chrY.vcf.gz" || true &
    
    wait  # espera que termine este genoma antes de pasar al siguiente
    echo "Genoma $code dividido y indexado"
done

echo "¡Todos los 18 genomas divididos por cromosoma en $SPLIT_DIR!"
___________________________________________________________________________________________________________________________________________________________#MERGED#
#panel1199 chrN
+
muestras chrN
↓
intersección de sitios comunes
↓
armonización exacta CHROM POS REF ALT
↓
panel final
↓
merge panel + muestra
↓
VCF merged con 1200 muestras#
___________________________________________________________________________________________________________________________________________________________#MERGED 017#


#!/usr/bin/env bash
set -euo pipefail

# =========================
# RUTAS PARA 017
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

  [[ -s "$PANEL" ]] || { echo "🚩 No existe PANEL: $PANEL"; exit 1; }
  [[ -s "$QUERY" ]] || { echo "🚩 No existe QUERY: $QUERY"; exit 1; }

  # Limpiar salidas previas
  rm -f chr${chr}.${SAMPLE}.common.vcf.gz chr${chr}.${SAMPLE}.common.vcf.gz.tbi \
        chr${chr}.panel1199.final.vcf.gz chr${chr}.panel1199.final.vcf.gz.tbi \
        chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz \
        chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz.tbi \
        sites.pos sites${SAMPLE}.4col.sorted panel.header.vcf panel.body.filtered \
        0000.vcf.gz 0000.vcf.gz.tbi 0001.vcf.gz 0001.vcf.gz.tbi README.txt sites.txt

  # 1. Intersección panel + 017
  bcftools isec -n=2 -w1 -O z -p "$odir" "$PANEL" "$QUERY"

  [[ -s sites.txt ]] || { echo "🚩 sites.txt no se generó"; exit 1; }
  nsites=$(wc -l < sites.txt)
  echo "sites.txt: $nsites"

  # 2. Filtrar 017 por posiciones comunes
  cut -f1,2 sites.txt > sites.pos

  bcftools view -T sites.pos -Oz \
    -o chr${chr}.${SAMPLE}.common.vcf.gz \
    "$QUERY"

  tabix -f -p vcf chr${chr}.${SAMPLE}.common.vcf.gz

  nquery=$(bcftools view -H chr${chr}.${SAMPLE}.common.vcf.gz | wc -l)
  echo "${SAMPLE}.common: $nquery"

  [[ "$nquery" -eq "$nsites" ]] || { echo "🚩 ${SAMPLE}.common != sites.txt"; exit 1; }

  # 3. Extraer sitios exactos CHROM POS REF ALT desde 017
  bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' chr${chr}.${SAMPLE}.common.vcf.gz \
    | sort -u > sites${SAMPLE}.4col.sorted

  # 4. Construir panel final con coincidencia exacta de alelos
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

  [[ "$npanel" -eq "$nsites" ]] || { echo "🚩 panel.final != sites.txt"; exit 1; }
  [[ "$nsamp_panel" -eq 1199 ]] || { echo "🚩 panel.final samples != 1199"; exit 1; }

  # 5. Merge panel + 017
  bcftools merge -m none -Oz \
    -o chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz \
    chr${chr}.panel1199.final.vcf.gz \
    chr${chr}.${SAMPLE}.common.vcf.gz

  tabix -f -p vcf chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz

  nmerged=$(bcftools view -H chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz | wc -l)
  nsamp_merged=$(bcftools query -l chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz | wc -l)

  echo "merged: $nmerged | samples(merged)=$nsamp_merged"

  [[ "$nmerged" -eq "$nsites" ]] || { echo "🚩 merged != sites.txt"; exit 1; }
  [[ "$nsamp_merged" -eq 1200 ]] || { echo "🚩 merged samples != 1200"; exit 1; }

  echo "chr${chr} ✅ OK"
done

echo "======== RESUMEN FINAL ${SAMPLE} ========"

for chr in {1..22}; do
  dir="${OUTBASE}/isec_chr${chr}"
  nsites=$(wc -l < "${dir}/sites.txt" 2>/dev/null || echo 0)
  nmerged=$(bcftools view -H "${dir}/chr${chr}.panel1199_plus_${SAMPLE}.common.merged.vcf.gz" 2>/dev/null | wc -l)
  echo "chr${chr} sites=${nsites} merged=${nmerged}"
done





___________________________________________________________________________________________________________________________________________________________#MERGED 018#


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
    """Run a shell command and raise if it fails."""
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
            print(f"🚩 No existe PANEL: {panel}", file=sys.stderr)
            sys.exit(1)
        if not query.exists() or query.stat().st_size == 0:
            print(f"🚩 No existe QUERY: {query}", file=sys.stderr)
            sys.exit(1)

        # limpiar outputs viejos
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

        # 1) isec
        run([
            "bcftools", "isec", "-n=2", "-w1", "-O", "z",
            "-p", str(odir),
            str(panel), str(query)
        ])

        sites_txt = odir / "sites.txt"
        if not sites_txt.exists() or sites_txt.stat().st_size == 0:
            print("🚩 sites.txt no se generó", file=sys.stderr)
            sys.exit(1)

        nsites = count_sites_txt(sites_txt)
        print(f"sites.txt: {nsites}")

        # 2) filtrar 018 por POS
        sites_pos = odir / "sites.pos"
        run(["bash", "-lc", f"cut -f1,2 {sites_txt} > {sites_pos}"])

        vcf_018_common = odir / f"chr{chr_}.018.common.vcf.gz"
        run(["bcftools", "view", "-T", str(sites_pos), "-Oz", "-o", str(vcf_018_common), str(query)])
        run(["tabix", "-f", "-p", "vcf", str(vcf_018_common)])

        n018 = count_vcf_records(vcf_018_common)
        print(f"018.common: {n018}")
        if n018 != nsites:
            print("🚩 018.common != sites.txt", file=sys.stderr)
            sys.exit(1)

        # 3) lista exacta 4col desde 018.common
        sites018_4col = odir / "sites018.4col.sorted"
        run(["bash", "-lc", f"bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_018_common} | sort -u > {sites018_4col}"])

        # 4) construir panel final exacto por alelos (incluye FORMAT=GT)
        panel_header = odir / "panel.header.vcf"
        run(["bash", "-lc", f"bcftools view -h {panel} > {panel_header}"])

        panel_body = odir / "panel.body.filtered"
        # Nota: se añade 'GT' y luego [\t%GT] para mantener el número correcto de columnas
        cmd = (
            f"bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\tGT[\\t%GT]\\n' {panel} | "
            f"awk 'BEGIN{{FS=OFS=\"\\t\"}} "
            f"NR==FNR {{ key[$1\"\\t\"$2\"\\t\"$3\"\\t\"$4]=1; next }} "
            f"{{ k=$1\"\\t\"$2\"\\t\"$3\"\\t\"$4; if(k in key){{ "
            f"for(i=5;i<=NF;i++) printf \"%s%s\", $i, (i==NF?ORS:OFS) "
            f"}} }}' {sites018_4col} - > {panel_body}"
        )
        run(["bash", "-lc", cmd])

        panel_final = odir / f"chr{chr_}.panel1199.final.vcf.gz"
        run(["bash", "-lc", f"cat {panel_header} {panel_body} | bgzip -c > {panel_final}"])
        run(["tabix", "-f", "-p", "vcf", str(panel_final)])

        npanel = count_vcf_records(panel_final)
        nsamp_panel = count_samples(panel_final)
        print(f"panel.final: {npanel} | samples(panel)={nsamp_panel}")

        if npanel != nsites:
            print("🚩 panel.final != sites.txt", file=sys.stderr)
            sys.exit(1)
        if nsamp_panel != 1199:
            print("🚩 panel.final samples != 1199", file=sys.stderr)
            sys.exit(1)

        # 5) merge (panel + 018)
        merged = odir / f"chr{chr_}.panel1199_plus_018.common.merged.vcf.gz"
        run([
            "bcftools", "merge", "-m", "none", "-Oz",
            "-o", str(merged),
            str(panel_final),
            str(vcf_018_common),
        ])
        run(["tabix", "-f", "-p", "vcf", str(merged)])

        nmerged = count_vcf_records(merged)
        nsamp_merged = count_samples(merged)
        print(f"merged: {nmerged} | samples(merged)={nsamp_merged}")

        if nmerged != nsites:
            print("🚩 merged != sites.txt", file=sys.stderr)
            sys.exit(1)
        if nsamp_merged != 1200:
            print("🚩 merged samples != 1200", file=sys.stderr)
            sys.exit(1)

        print(f"chr{chr_} ✅ OK")

    # Resumen final
    print("======== RESUMEN FINAL (018) ========")
    for chr_ in CHRS:
        dir_ = OUTBASE / f"isec_chr{chr_}"
        if not dir_.exists():
            print(f"chr{chr_} 🚩 faltante")
            continue
        sites_txt = dir_ / "sites.txt"
        merged = dir_ / f"chr{chr_}.panel1199_plus_018.common.merged.vcf.gz"
        nsites = count_sites_txt(sites_txt)
        nmerged = count_vcf_records(merged) if merged.exists() else 0
        print(f"chr{chr_} sites={nsites} merged={nmerged}")

if __name__ == "__main__":
    main()



__________________________________________________________________________________________________________________________________#ARCHIVOS_RFMIX#
#PREPARACIÓN ARCHIVOS RFMIX#


_________________________________________________________________________________________________________________________________#018#
#.query_ref_localitation_018# 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from cyvcf2 import VCF

# ================= CONFIG =================
CHRS = range(1, 23)

# Carpeta base de tu muestra 018 (donde están isec_chr1..isec_chr22)
BASE = Path(
    "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/"
    "18genomes/18genfased/018/with_ref_panel1199"
)

# Lista de IDs del panel (1199)
PANEL_IDS = Path(
    "/mnt/diskrare/ivonb/refamerindios/panel_hgdp1kg_1200/panel_1199.ids"
)

# Mapa (IDs -> superpoblación), 2 columnas: ID <tab/space> SUPERPOP
SUPERPOP_MAP = Path("/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map")

# Directorio con genetic_map_chrN_rfmix.txt
MAPDIR = Path("/home/rare/ivon/outputs/rfmix/rfmix_map")
# ==========================================


def run(cmd):
    print("▶", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True)


def load_genetic_map(chr_str: str):
    gmap = MAPDIR / f"genetic_map_chr{chr_str}_rfmix.txt"
    if not gmap.exists():
        raise FileNotFoundError(f"No existe mapa genético: {gmap}")

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
    Salida en 3 columnas (formato que te funcionó):
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
        raise RuntimeError(f"No se encontraron SNPs bialélicos en {vcf_path}")

    out_path.write_text("".join(out_lines))
    print(f"   ✅ snp_locations: {out_path.name} | SNPs={kept} | Total_VCF={total}")


def detect_query_id(merged_vcf: Path) -> str:
    """
    En tus merges, el último sample suele ser tu muestra (como viste en 017 y 018).
    """
    samples = subprocess.check_output(
        ["bcftools", "query", "-l", str(merged_vcf)],
        text=True
    ).strip().splitlines()

    if not samples:
        raise RuntimeError(f"VCF sin samples en header: {merged_vcf}")

    return samples[-1]


def make_labels_file_for_ref(ref_vcf: Path, out_labels: Path):
    """
    Crea chrN.superpopulation_labels.txt con 2 columnas:
      <ID> <SUPERPOP>
    Para TODOS los samples del ref_vcf (1199).
    """
    # Leer orden real de samples del REF (esto es lo más importante)
    ref_ids = subprocess.check_output(
        ["bcftools", "query", "-l", str(ref_vcf)],
        text=True
    ).strip().splitlines()

    if not ref_ids:
        raise RuntimeError(f"REF sin samples: {ref_vcf}")

    # Cargar mapa ID->superpop
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

    # Escribir labels en el mismo orden del ref_vcf
    na = 0
    with out_labels.open("w") as out:
        for sid in ref_ids:
            sp = pop.get(sid, "NA")
            if sp == "NA":
                na += 1
            out.write(f"{sid}\t{sp}\n")

    print(f"   ✅ labels: {out_labels.name} | n={len(ref_ids)} | NA={na}")


def main():
    for chrn in CHRS:
        chr_str = str(chrn)
        idir = BASE / f"isec_chr{chr_str}"
        merged = idir / f"chr{chr_str}.panel1199_plus_018.common.merged.vcf.gz"

        if not merged.exists():
            print(f"⚠ chr{chr_str}: no existe {merged.name}, se omite")
            continue

        print(f"\n=========== chr{chr_str} (018) ===========")
        query_id = detect_query_id(merged)
        print("   Query ID detectado:", query_id)

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

        # ---------- SNP locations (3 columnas) ----------
        snp_loc = idir / f"chr{chr_str}.snp_locations.fixed"
        make_snp_locations_3col(merged, chr_str, snp_loc)

        # ---------- Labels (REF) ----------
        labels = idir / f"chr{chr_str}.superpopulation_labels.txt"
        make_labels_file_for_ref(ref_vcf, labels)

        print(f"✔ chr{chr_str} listo: ref | query | snp_locations.fixed | labels")

    print("\n✅ Listo. Ya tienes por cromosoma: chrN.ref.vcf.gz, chrN.query.vcf.gz, chrN.snp_locations.fixed, chrN.superpopulation_labels.txt")


if __name__ == "__main__":
    main()



_____________________________________________________________________________________________________________________________________________________#017#
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
    print("▶", " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_genetic_map(chr_str):
    gmap = MAPDIR / f"genetic_map_chr{chr_str}_rfmix.txt"
    if not gmap.exists():
        raise FileNotFoundError(gmap)

    df = pd.read_csv(
        gmap, sep=r"\s+", header=None,
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
        raise RuntimeError(f"No SNPs en {vcf_path}")

    out_path.write_text("".join(out))


def main():
    for chrn in CHRS:
        chr_str = str(chrn)
        idir = BASE / f"isec_chr{chr_str}"
        merged = idir / f"chr{chr_str}.panel1199_plus_017.common.merged.vcf.gz"

        if not merged.exists():
            print(f"⚠ chr{chr_str}: no existe merged, se omite")
            continue

        print(f"\n=========== chr{chr_str} ===========")

        # Detectar query real
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

        print(f"✔ chr{chr_str} listo: ref | query | snp_locations")


if __name__ == "__main__":
    main()


_______________________________________________________________________________________________________________________________#RUN_RFMIX#
# RFMIX#
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


#PLOTS#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Rutas (tuyas) ===
BASE_DIR = "/home/rare/ivon/data/vcf_filtrados/newfil/bialelicos/18genomes/plink_individual"
PREFIX = "merged_18_1200_common_qc05_maf05_nomind_pruned"

# Tu mapa IID -> SUPERPOP (2 columnas: IID \t AFR/EUR/EAS/AMR)
MAP_SUPERPOP = "/mnt/diskrare/ivonb/refamerindios/superpop_4groups.map"

# Ks que ya corriste (ajusta si quieres)
K_LIST = [1, 2, 3, 4, 5, 6]

# Orden de bloques en el plot
BLOCK_ORDER = ["AFR", "EUR", "EAS", "AMR", "OWN", "UNK"]

# === Utilidades ===
def short_pb(iid: str) -> str:
    # PB000696_017C... -> "017" (y 018C... -> "018")
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
        raise ValueError(f"Q tiene {q.shape[1]} columnas pero K={k}. Revisa el archivo .Q")

def plot_one_k(df_ordered: pd.DataFrame, Q: np.ndarray, k: int, out_png: str):
    # Barras apiladas
    x = np.arange(df_ordered.shape[0])
    bottom = np.zeros(df_ordered.shape[0])

    plt.figure(figsize=(20, 6))

    for j in range(k):
        plt.bar(x, Q[:, j], bottom=bottom, width=1.0)
        bottom += Q[:, j]

    # Separadores entre bloques
    # Calcula cortes donde cambia el BLOCK
    blocks = df_ordered["BLOCK"].values
    cuts = np.where(blocks[:-1] != blocks[1:])[0]
    for c in cuts:
        plt.axvline(c + 0.5, linewidth=1)

    # Etiquetas en el eje X: solo 18 genomas (SHORT), y resaltar 017/018
    plt.xticks([], [])  # quitamos ticks por defecto

    # Colocar texto solo para tus PB
    pb_mask = df_ordered["IID"].str.startswith("PB000696_").fillna(False).values
    for i in np.where(pb_mask)[0]:
        lab = df_ordered.iloc[i]["SHORT"]
        # resalta 017 y 018
        if lab in ("017", "018"):
            plt.text(i, 1.02, lab, ha="center", va="bottom", fontsize=10, fontweight="bold", rotation=90)
        else:
            plt.text(i, 1.02, lab, ha="center", va="bottom", fontsize=8, rotation=90)

    # Títulos de bloques arriba
    # posición media de cada bloque
    start = 0
    for b in BLOCK_ORDER:
        idx = np.where(df_ordered["BLOCK"].values == b)[0]
        if len(idx) == 0:
            continue
        mid = (idx[0] + idx[-1]) / 2
        plt.text(mid, 1.10, f"{b} (n={len(idx)})", ha="center", va="bottom", fontsize=11)
        start = idx[-1] + 1

    plt.ylim(0, 1.15)
    plt.ylabel("Proporción (Q)")
    plt.title(f"ADMIXTURE (supervised/unsupervised) – {PREFIX} – K={k}\nOrdenado por superpoblación + 18 genomas (017 y 018 resaltados)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

# === Main ===
def main():
    os.chdir(BASE_DIR)

    fam_path = f"{PREFIX}.fam"
    fam = load_fam_ids(fam_path)

    mp = load_map(MAP_SUPERPOP)

    # Merge manteniendo el orden del .fam
    df = fam.merge(mp, on="IID", how="left")

    # Clasificación de bloque:
    df["BLOCK"] = df["SUPERPOP"].fillna("UNK")
    # tus 18
    mask_own = df["IID"].str.startswith("PB000696_").fillna(False)
    df.loc[mask_own, "BLOCK"] = "OWN"

    # Etiqueta corta para PB
    df["SHORT"] = df["IID"].apply(short_pb)

    # Orden por bloque (AFR/EUR/EAS/AMR/OWN/UNK) y dentro de cada bloque por IID
    df["BLOCK"] = pd.Categorical(df["BLOCK"], categories=BLOCK_ORDER, ordered=True)
    df_sorted = df.sort_values(["BLOCK", "IID"]).reset_index(drop=True)

    # Reordenar Q para que coincida con df_sorted:
    # Necesitamos mapear índice original (orden .fam) -> nuevo orden
    idx_map = pd.Series(np.arange(df.shape[0]), index=df["IID"]).to_dict()
    original_positions = df_sorted["IID"].map(idx_map).values

    print("Conteo por bloque:")
    print(df_sorted["BLOCK"].value_counts(dropna=False))

    for k in K_LIST:
        qfile = f"{PREFIX}.{k}.Q"
        if not os.path.exists(qfile) or os.path.getsize(qfile) == 0:
            print(f"[WARN] No existe: {qfile} (lo salto)")
            continue

        q = pd.read_csv(qfile, sep=r"\s+", header=None)
        ensure_q_shape(q, k)

        Q = q.values.astype(float)
        Q_sorted = Q[original_positions, :]

        out_png = os.path.join(BASE_DIR, f"{PREFIX}.K{k}.Q.superpop_blocks.png")
        plot_one_k(df_sorted, Q_sorted, k, out_png)
        print(f"[OK] Guardado: {out_png}")

if __name__ == "__main__":
    main()































