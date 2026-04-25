#1. Filtrar snps_ bialelicos_autosomas_pass_#


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




# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #






























