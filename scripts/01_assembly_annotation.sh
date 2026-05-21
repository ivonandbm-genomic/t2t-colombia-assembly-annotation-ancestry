#libraries
import os, sys
import subprocess
import datetime
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pycircos import Circos
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob 
import csv 
import os
import glob
import shutil
import re



#routes
smrtlink="../programs/smrtlink/smrtlink/smrtcmds/bin"
output="./outputs"
folder_genomes="../data/genomes/bamfiles"
folder_reference="./reference"
reference_t2t="./reference/chm13v22.fasta"
genome_t2t="t2t"
genome_hg38="hg38"
reference_hg38="./reference/hg38.fasta"
bamfiles="./bamfiles"
mnt_diskrare="/mnt/diskrare/arlenb"
genes_t2t="/home/rare/arlen/reference/genes_t2t_renamed.gtf"
reference_17C="./reference/03_1_C01_bc2059_017C_p_ctg.fa"
genome_018C="018C"
genome_017C="017C"
reference_18C="./reference/03_1_D01_bc2060_018C_p_ctg.fa"

def log_output(function_name, output):
    with open("register.txt", "a") as file:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        file.write(f"[{timestamp}] Function: {function_name} | Output: {output}\n")

def index_bamfiles (folder_genomes, smrtlink):
    for x in os.listdir(folder_genomes):
        if x.endswith(".bam"):
            file=os.path.join(folder_genomes,x)
            command=f"{smrtlink}/pbindex {file}"
            subprocess.run(command, shell=True, check=True)


def convert_bam_fastq(folder_genomes, output,smrtlink):
    out=f"{mnt_diskrare}/01"
    if not os.path.exists(out):
        os.makedirs(out)
    for x in os.listdir(folder_genomes):
        if x=="1_A01_bc2044_002P.bam":
            file=os.path.join(folder_genomes,x)
            basename=os.path.splitext(x)[0]
            command=f"{smrtlink}/bam2fastq -o {out}/01_{basename} {file} -j 32"
            subprocess.run(command, shell=True, check=True)
            

def run_all_jellyfish(mnt_diskrare, output, k=21, hash_size="10G", threads=16):
    input_folder=f"{mnt_diskrare}/01"
    gz_files = glob.glob(os.path.join(input_folder, "*.fastq.gz"))
    out=f"{output}/01"
    if not os.path.exists(out):
        os.makedirs(out)
    for gz in gz_files:
        sample = os.path.basename(gz).replace(".fastq.gz", "")
        out_jf = f"{out}/{sample}.jf"
        out_histo = f"{out}/{sample}.histo"
        cmd = (
          f"zcat {gz} | jellyfish count -C -m {k} "
          f"-s {hash_size} -t {threads} -o {out_jf} /dev/fd/0"
        )
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(
          f"jellyfish histo -t {threads} {out_jf} > {out_histo}",
          shell=True, check=True
        )
        print(f"✅ Done {sample}")

#run_all_jellyfish(mnt_diskrare, output, k=31, hash_size="10G", threads=32)

def hifi_assembly(output, mnt_diskrare):
    folder_fastq_files = f"{mnt_diskrare}/01"
    num_cores = os.cpu_count()
    out = f"{output}/02"
    os.makedirs(out, exist_ok=True)

    for x in os.listdir(folder_fastq_files):
        if x.endswith("fastq.gz"):
            file = os.path.join(folder_fastq_files, x)
            x = x.replace("01_1_", "")
            basename = x.split(".")[0]

            # Check if assembly already exists (use .p_ctg.gfa as marker)
            output_gfa = os.path.join(out, f"02_1_{basename}.asm.p_ctg.gfa")
            if os.path.exists(output_gfa):
                print(f"Skipping {basename}, already assembled.")
                continue

            print(f"Assembling: {basename}")
            program = "../programs/02/hifiasm"
            command = f"{program} -o {out}/02_1_{basename}.asm -t {num_cores} {file}"
            #subprocess.run(command, shell=True, check=True)



def convert_gfa_to_fasta(output):
    folder_gfa_files=f"{output}/02"
    out=f"{output}/03"
    if not os.path.exists(out):
        os.makedirs(out)
    for x in os.listdir(folder_gfa_files):
        if x.endswith("p_ctg.gfa"):
            file=os.path.join(folder_gfa_files, x)
            x=x.replace("02_","")
            basename=x.split(".")[0]
            haplotype=x.split(".")[3]
            program="../programs/03/gfatools"
            command=f"{program} gfa2fa  {file} > {out}/03_{basename}_{haplotype}.fasta"
            #print(command)
            subprocess.run(command, shell=True, check=True)
            #log_output("convert_gfa_to_fasta", f"successfully ran: {command}")

#Align the contigs assembly  with hifi reads

def align_contig_assembly(output, mnt_diskrare, smrtlink):
    folder_hifiR=f"{mnt_diskrare}/01"
    folder_assemblies=f"{output}/03"
    out=f"{output}/03_aligned"
    if not os.path.exists(out):
        os.makedirs(out)
    for x in os.listdir(folder_assemblies):
        if x.endswith("P_p_ctg.fasta"):
            fasta_file=os.path.join(folder_assemblies,x)
            basename=os.path.basename(x).replace("_p_ctg.fasta", "").replace("03_", "01_")
            hifi_reads=os.path.join(folder_hifiR,f"{basename}.fastq.gz")
            command=f"{smrtlink}/pbmm2 align {fasta_file} {hifi_reads} {out}/{basename}.bam --sort --preset ccs -j 32"
            subprocess.run(command, check=True, shell=True)
#align_contig_assembly(output, mnt_diskrare, smrtlink)



def decompress_fastq(gz_path, out_path):
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def correct_with_ragtag(output, mnt_diskrare, reference, genome):
    folder_hifiR = f"{mnt_diskrare}/01"
    folder_assemblies = f"{output}/03"
    out = f"{output}/03_ragtag_correct/{genome}"
    
    os.makedirs(out, exist_ok=True)

    for x in os.listdir(folder_assemblies):
        if x.endswith("P_p_ctg.fasta"):
            fasta_file = os.path.join(folder_assemblies, x)
            basename = os.path.basename(x).replace("_p_ctg.fasta", "").replace("03_", "01_")
            hifi_gz = os.path.join(folder_hifiR, f"{basename}.fastq.gz")
            hifi_fastq = os.path.join(folder_hifiR, f"{basename}.fastq")

            # Decompress if necessary
            if os.path.exists(hifi_gz):
                decompress_fastq(hifi_gz, hifi_fastq)
            else:
                raise FileNotFoundError(f"Missing FASTQ.gz file: {hifi_gz}")

            ragtag_out = os.path.join(out, basename)
            os.makedirs(ragtag_out, exist_ok=True)

            # Run ragtag correct
            program = "/home/rare/programs/04/ragtag.py"
            command = f"{program} correct -t 32 -u -F {hifi_fastq} -T corr -o {ragtag_out} {reference} {fasta_file}"
            subprocess.run(command, shell=True, check=True)

            # Clean up the uncompressed FASTQ
            os.remove(hifi_fastq)

#correct_with_ragtag(output, mnt_diskrare, reference_t2t, genome_t2t)
#correct_with_ragtag(output, mnt_diskrare, reference_hg38, genome_hg38)


#Scaffolding

TOOLS = {
    "ragtag": "/home/rare/programs/04/ragtag.py"
}
INPUT_DIR = "/home/rare/arlen/outputs/03"
OUTPUT_ROOT = "/home/rare/arlen/outputs/04/scaffolds"
THREADS = 32

def run(cmd):
    print(">>>", cmd)
    subprocess.run(cmd, shell=True, check=True)

def scaffold_with_reference(reference_name, reference_path):
    pattern = os.path.join(INPUT_DIR, "*.fasta")
    for fasta in glob.glob(pattern):
        prefix = os.path.splitext(os.path.basename(fasta))[0]
        outdir = os.path.join(OUTPUT_ROOT, reference_name, prefix)
        os.makedirs(outdir, exist_ok=True)

        cmd = (
            f"{TOOLS['ragtag']} scaffold "
            f"{reference_path} {fasta} "
            f"-o {outdir} -t {THREADS}"
        )
        run(cmd)

        src = os.path.join(outdir, "ragtag.scaffold.fasta")
        dst = os.path.join(outdir, f"{prefix}.fa")
        if os.path.exists(src):
            os.rename(src, dst)
        else:
            print(f"Warning: scaffold missing for {prefix} vs {reference_name}")

def scaffold_t2t():
    scaffold_with_reference("t2t", "/home/rare/arlen/reference/chm13v22.fasta")

def scaffold_hg38():
    scaffold_with_reference("hg38", "/home/rare/arlen/reference/hg38.fasta")

# Immediately run both functions, without any conditional guard
#scaffold_t2t()
#scaffold_hg38()


def clean_fasta_headers(root_dir):
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.endswith(".fa"):
                fasta_path = os.path.join(dirpath, fname)
                temp_path = fasta_path + ".tmp"

                with open(fasta_path, 'r') as fin, open(temp_path, 'w') as fout:
                    for line in fin:
                        if line.startswith(">"):
                            # Keep only the part before the first underscore
                            prefix = line[1:].split("_", 1)[0]
                            new_header = f">{prefix}\n"
                            fout.write(new_header)
                        else:
                            fout.write(line)

                os.replace(temp_path, fasta_path)
                print(f"✅ Cleaned: {fasta_path}")

# Example usage
#clean_fasta_headers("./outputs/04/scaffolds")



#Genome evaluation

def quast_contig(reference,output, genome):
    folder_contigs=os.path.join(output, "03")
    out_1=f"{output}/05/quast_contigs_no_reference/{genome}"
    if not os.path.exists(out_1):
        os.makedirs(out_1)
    out_2=f"{output}/05/quast_contigs_reference/{genome}"
    if not os.path.exists(out_2):
        os.makedirs(out_2)    
    for x in os.listdir(folder_contigs):
        if x.endswith(".fasta"):
            file=os.path.join(folder_contigs,x)
            basename=os.path.basename(x).replace("03_", "")
            program="../programs/05/quast-lg.py"
            command_1=f"{program} {file}  -o {out_1}/05_{basename} --large -t 32"
            subprocess.run(command_1, shell=True, check=True)
            command_2=f"{program} {file} -r {reference}  -o {out_2}/05_{basename} -t 32 --large  --eukaryote --circos"
            subprocess.run(command_2, shell=True, check=True)

#quast_contig(reference_t2t, output, genome_t2t)
#quast_contig(reference_hg38, output, genome_hg38)

#Quality control 

def quast_scaffolds(reference, output, genome,genes, mnt_diskrare):
    folder_scaffolds=os.path.join(output,"04", "scaffolds", genome)
    out=f"{output}/05/quast_scaffolds/{genome}"
    if not os.path.exists(out):
        os.makedirs(out)
    for root, dirs,  files in os.walk(folder_scaffolds):
        for x in files:
            if  x.endswith(".fasta"):
                file=os.path.join(root,x)
                basename=re.sub(r'(_hap1|_hap2|_p_ctg)?\.fasta$', '', os.path.basename(x).replace("04_", ""))
                program="../programs/05/quast-lg.py"
                command=f"{program} {file} -r {reference} --pacbio {mnt_diskrare}/01/01_1_{basename}.fastq.gz -g {genes} -o {out}/05_{basename} --upper-bound-min-con 1 -L --eukaryote  --large --circos  -t 32" 
                #print(command)
                subprocess.run(command, shell=True, check=True)
        
 
#quast_scaffolds(reference_t2t, output, genome_t2t, genes_t2t, mnt_diskrare)

def busco(output,mnt_diskrare, genome):
    folder_scaffolds=os.path.join(mnt_diskrare,"04", "scaffolds", genome)
    out=f"{output}/busco/{genome}"
    os.makedirs(out, exist_ok=True)
    for root, dirs,  files in os.walk(folder_scaffolds):
        for x in files:
            if x.endswith("_hap2.fasta"):
                file=os.path.join(root,x)
                basename=x.replace(".fasta","")
                program="../programs/06/busco/bin/busco"
                command=[f"{program}", 
                         "-i", f"{file}",
                         "-o", f"{out}/{basename}",
                         "-m", "genome",
                         "-l", "primates_odb12",
                         "-c", "24",
                         "-f", 
                         "--config", "/home/rare/programs/06/busco/config/config.ini"]
                #print(command)
                subprocess.run(command, check=True)

#busco(output, mnt_diskrare, genome_t2t)

def busco_plot(output, genome):
    folder_busco_results=f"{output}/busco/{genome}"
    all_results=f"{output}/busco/all_results_pctg_Control"
    os.makedirs(all_results, exist_ok=True)
    for root, dirs, files in os.walk(folder_busco_results):
        for x in files:
            if x.endswith("C_p_ctg.json"):
                file=os.path.join(root,x)
                command_1=["cp", f"{file}",
                         f"{all_results}"]
                subprocess.run(command_1, check=True)
                program="../programs/06/busco/bin/busco"
                command_2=["busco", "--plot", f"{all_results}"]
                subprocess.run(command_2,check=True)
#busco_plot(output,genome_t2t)



# merqury + meryl
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import shutil
import subprocess as sp
from pathlib import Path
from typing import List

# ──────────────────────────────────────────
# CONFIG: ajusta sólo lo necesario
# ──────────────────────────────────────────
ROOT = Path("/mnt/diskrare/arlenb")
INPUT_DIRS = [
    ROOT / "04/scaffolds/t2t/03_1_C01_bc2059_017C_p_ctg",
    ROOT / "04/scaffolds/t2t/03_1_D01_bc2060_018C_p_ctg",
    ROOT / "04/scaffolds/hg38/03_1_C01_bc2059_017C_p_ctg",
    ROOT / "04/scaffolds/hg38/03_1_D01_bc2060_018C_p_ctg",
]
DATA_ROOT = Path("/home/rare/ivon/data/merqury")

# merqury + meryl
MERQURY = Path(os.environ.get("MERQURY", "/home/rare/programs/merqury-1.3"))
MERYL_BIN = shutil.which("meryl") or "/home/rare/programs/meryl-src/build/bin/meryl"
SAMTOOLS = shutil.which("samtools") or "samtools"  # en PATH

# Parámetros de cómputo (ajusta MEM_GB según tu RAM; tienes ~187 GB → 128–140 va bien)
K = 21
THREADS = 32
MEM_GB = 128     # límite de RAM para meryl (GB aprox)
DRY_RUN = False  # pon True si quieres ver comandos sin ejecutarlos

# Lecturas por muestra: AJUSTA rutas/patrones si cambian
READS_BY_SAMPLE = {
    "017C": "/mnt/diskrare/arlenb/01/01_1_C01_bc2059_017C.fastq.gz",
    "018C": "/mnt/diskrare/arlenb/01/01_1_D01_bc2060_018C.fastq.gz",
}

# ──────────────────────────────────────────
# utilidades
# ──────────────────────────────────────────
def run(cmd: List[str], cwd: Path | None = None):
    print("→", " ".join(map(str, cmd)), f"(cwd={cwd})" if cwd else "")
    if DRY_RUN:
        return
    sp.run(cmd, check=True, cwd=cwd)

def find_fasta(d: Path) -> Path | None:
    cands = []
    for pat in ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"):
        cands.extend(d.glob(pat))
    pri = [x for x in cands if "p_ctg" in x.name]
    if pri:
        return pri[0]
    return cands[0] if cands else None

def sample_from_path(p: Path) -> str | None:
    s = str(p)
    if "017C" in s: return "017C"
    if "018C" in s: return "018C"
    m = re.search(r"(\d{3}C)", s)
    return m.group(1) if m else None

def context_from_path(p: Path) -> str:
    s = str(p)
    if "/t2t/" in s: return "t2t"
    if "/hg38/" in s: return "hg38"
    return "asm"

def ensure_tools():
    if not Path(MERYL_BIN).exists():
        raise SystemExit(f"❌ No encuentro meryl en {MERYL_BIN} (agrega al PATH o ajusta MERYL_BIN)")
    if not (MERQURY / "merqury.sh").exists():
        raise SystemExit(f"❌ No encuentro merqury.sh en {MERQURY}")
    if shutil.which(SAMTOOLS) is None and not shutil.which("samtools"):
        raise SystemExit("❌ No encuentro samtools en el PATH (necesario para faidx)")
    print(f"✓ meryl: {MERYL_BIN}")
    print(f"✓ merqury.sh: {MERQURY / 'merqury.sh'}")
    print(f"✓ samtools en PATH")

def ensure_reads(sample: str) -> List[str]:
    pattern = READS_BY_SAMPLE.get(sample)
    if not pattern:
        raise SystemExit(f"❌ No definiste patrón de lecturas para {sample} en READS_BY_SAMPLE")
    hits = glob.glob(pattern)
    if not hits:
        raise SystemExit(f"❌ No encontré lecturas con patrón: {pattern}")
    return hits

def check_disk_space(path: Path, need_gb: int = 200):
    try:
        usage = shutil.disk_usage(path)
        free_gb = usage.free / (1024**3)
        if free_gb < need_gb:
            print(f"⚠️  Espacio libre bajo en {path.resolve()}: {free_gb:.1f} GB (sug. ≥ {need_gb} GB)")
    except Exception as e:
        print(f"ℹ️  No pude medir espacio en {path}: {e}")

def is_meryl_db_ok(db: Path) -> bool:
    if not db.exists(): return False
    expected = ["histogram", "info", "README"]
    return any((db / x).exists() for x in expected)

def ensure_faidx(fa: Path):
    """Crea .fai si no existe."""
    fai = Path(str(fa) + ".fai")
    if not fai.exists():
        run([SAMTOOLS, "faidx", str(fa)])

# ──────────────────────────────────────────
# pipeline por ensamblaje
# ──────────────────────────────────────────
def process_dir(d: Path):
    asm = find_fasta(d)
    if not asm:
        print(f"⚠ No FASTA en {d}, salto.")
        return

    sample = sample_from_path(d)
    ctx = context_from_path(d)
    if not sample:
        print(f"⚠ No pude deducir muestra (017C/018C) desde {d}, salto.")
        return

    reads = ensure_reads(sample)

    # carpeta de trabajo por muestra/contexto
    work = DATA_ROOT / "merqury" / f"{sample}_{ctx}"
    work.mkdir(parents=True, exist_ok=True)
    (work / "logs").mkdir(exist_ok=True)
    check_disk_space(work, need_gb=200)

    # prefijo relativo (IMPORTANTÍSIMO para merqury.sh)
    outprefix = f"{sample}__{ctx}"

    # 1) meryl count (limita memoria y hilos)
    meryl_db = work / f"{sample}.meryl"
    if is_meryl_db_ok(meryl_db):
        print(f"✓ meryl DB existe: {meryl_db}")
    else:
        if meryl_db.exists() and not DRY_RUN:
            shutil.rmtree(meryl_db, ignore_errors=True)
        cmd = [
            MERYL_BIN, "count",
            f"k={K}",
            f"memory={MEM_GB}",
            f"threads={THREADS}",
            "output", str(meryl_db),
        ] + reads
        run(cmd, cwd=work)

   # 2) Ensure assembly index (prevents .fasta.fai errors)
ensure_faidx(asm)

    # 3 Merqury (NO pasar -k ni -m aquí)
    qv_file = work / f"{outprefix}.qv"
    if qv_file.exists():
        print(f"✓ QV ya generado: {qv_file}")
    else:
        cmd = [
            str(MERQURY / "merqury.sh"),
            "-t", str(THREADS),
            str(meryl_db),
            str(asm),
            outprefix,
        ]
        run(cmd, cwd=work)


# ──────────────────────────────────────────
# main
# ──────────────────────────────────────────
def main():
    ensure_tools()
    print(f"Parámetros: k={K}, threads={THREADS}, mem≈{MEM_GB} GB")
    for d in INPUT_DIRS:
        print(f"\n=== Procesando {d} ===")
        try:
            process_dir(d)
        except sp.CalledProcessError as e:
            print(f"❌ Falló un comando en {d} (returncode={e.returncode}). Revisa logs en {DATA_ROOT}/merqury/*/logs/")
        except SystemExit as e:
            print(e)

if __name__ == "__main__":
    main()




#Repeat annotation

import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_repeat_masker_single(file, output, genome, threads_per_job=6):
    out = f"{output}/10/{genome}"
    basename = os.path.basename(file).replace(".fa", "")
    outdir = os.path.join(out, basename)

    # Skip if output folder already exists
    if os.path.exists(outdir):
        print(f"Skipping {file}: folder already exists at {outdir}")
        return

    os.makedirs(outdir, exist_ok=True)

    program = "../programs/repeats/RepeatMasker/RepeatMasker"
    cmd = [
        program,
        "-species", "human",
        "-dir", outdir,
        "-pa", str(threads_per_job),
        "-engine", "rmblast",
        "-gff",
        "-xsmall",
        file
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"Finished: {file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {file}: {e}")

def Repeat_Masker_parallel(output, genome, max_workers=5, threads_per_job=6):
    folder_assemblies = f"{output}/04/scaffolds/{genome}"
    input_files = []

    for root, dirs, files in os.walk(folder_assemblies):
        for file in files:
            if file.endswith("p_ctg.fa"):
                input_files.append(os.path.join(root, file))

    print(f"Total input files found: {len(input_files)}")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(run_repeat_masker_single, file, output, genome, threads_per_job)
            for file in input_files
        ]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Unexpected error: {e}")


#Repeat_Masker_parallel(output, genome_t2t, max_workers=5, threads_per_job=6)



#Genome annotation with Liftover  

import os
import subprocess
from Bio import SeqIO

def extract_chromosomes(fasta_path, filtered_path):
    """Keep only chr1–chr22, chrX, chrY from FASTA and write to new file"""
    allowed = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    with open(filtered_path, "w") as out_handle:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in allowed:
                SeqIO.write(record, out_handle, "fasta")

def liftover(reference, output, genome, gene_annotation):
    out = f"{output}/13/{genome}"
    out_2 = f"{output}/13/intermediate_files/{genome}"
    os.makedirs(out, exist_ok=True)
    os.makedirs(out_2, exist_ok=True)

    folder_assemblies = f"{output}/04/scaffolds/{genome}"

    for root, dirs, files in os.walk(folder_assemblies):
        for file in files:
            if file.endswith(".fasta"):
                original_fasta = os.path.join(root, file)

                if not os.path.isfile(original_fasta):
                    print(f"[✗] File not found: {original_fasta} — skipping.")
                    continue

                basename = os.path.basename(file).replace(".fa", "").replace(".fasta", "")
                filtered_fasta = os.path.join(root, f"{basename}.chroms.fa")

                try:
                    print(f"[→] Filtering chromosomes from: {original_fasta}")
                    extract_chromosomes(original_fasta, filtered_fasta)
                except Exception as e:
                    print(f"[✗] Failed to extract chromosomes from {original_fasta}: {e}")
                    continue

                out_gtf = os.path.join(out, f"{basename}.gtf")

                cmd = [
                    "liftoff",
                    filtered_fasta,
                    reference,
                    "-g", gene_annotation,
                    "-o", out_gtf,
                    "-p", "8",
                    "-infer_genes",
                    "-dir", out_2
                ]

                try:
                    print(f"[✓] Running Liftoff on: {basename} (chromosomes only)")
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"[✗] Liftoff failed on {basename}: {e}")
                    continue

                # Clean up
                try:
                    os.remove(filtered_fasta)
                except OSError as e:
                    print(f"[!] Warning: could not delete temp file {filtered_fasta}: {e}")

# Example call:
#liftover(reference_t2t, output, genome_t2t, genes_t2t)

def gene_annotation(mnt_diskrare, reference_t2t, output, genome_t2t, genes_t2t):
    out=f"{output}/liftoff/{genome_t2t}"
    os.makedirs(out, exist_ok=True)
    out2=f"{output}/liftoff/{genome_t2t}/intermediate_files"
    os.makedirs(out2, exist_ok=True)
    folder_assemblies=f"{mnt_diskrare}/04/scaffolds/{genome_t2t}"
    for root, dirs, files  in os.walk(folder_assemblies):
        for file in files:
            if file.endswith("p_ctg.fasta.wrap.fa"):
                assembly=os.path.join(root,file)
                basename=os.path.basename(assembly).replace(".fasta.wrap.fa", "")
                command=f"liftoff {assembly} {reference_t2t} -g {genes_t2t} -o {out}/{basename}.gtf -dir {out2}/{basename} -u {out}/{basename}_unlifted.gtf -p 24  -a 0.9 -s 0.9"
                subprocess.run(command, check=True, shell=True)
#gene_annotation(mnt_diskrare, reference_t2t, output, genome_t2t, genes_t2t)


#Structural Variant Calling using T2T and hg38 as references

def aligment(folder_genomes, reference, output, smrtlink, genome):
    out=f"{output}/08/aligned_reads/{genome}"
    if not os.path.exists(out):
        os.makedirs(out)
    for x in os.listdir(folder_genomes):
        if x.endswith(".bam"):
            file = os.path.join(folder_genomes, x)
            basename = os.path.splitext(x)[0]
            output_file = f"{out}/08_{basename}_aligned.bam"
            if not os.path.exists(output_file):
                command = (
                    f"{smrtlink}/pbmm2 align {reference} {file} "
                    f"--sort --preset HIFI --log-level INFO -j 32 > {output_file}"
                )
                subprocess.run(command, shell=True, check=True)
            else:
                print(f"Skipping alignment for {x} as {output_file} already exists.")
#aligment(folder_genomes, reference_17C, output, smrtlink, genome_017C)
#aligment(folder_genomes, reference_18C, output, smrtlink, genome_018C)

def discover_signatures(output, mnt_diskrare, smrtlink, genome):
    aligned_reads = f"{mnt_diskrare}/08/aligned_reads/{genome}"
    out = f"{output}/08/structural_variants/signatures/{genome}"

    os.makedirs(out, exist_ok=True)

    for x in os.listdir(aligned_reads):
        if x.endswith(".bam"):
            file = os.path.join(aligned_reads, x)
            basename = "_".join(x.split("_")[0:5])
            out_file = f"{out}/{basename}.svsig.gz"

            if os.path.exists(out_file):
                print(f"Skipping {basename}, signature already exists.")
                continue

            command = (
                f"{smrtlink}/pbsv discover {file} {out_file} "
                f"--tandem-repeats /home/rare/arlen/reference/human_T2T_CHM13v2.trf.bed"
            )
            print(f"Running: {command}")
            subprocess.run(command, shell=True, check=True)
#discover_signatures(output, mnt_diskrare, smrtlink, genome_t2t)


def structural_variants_call(output, reference, smrtlink, genome):
    signatures = f"{output}/08/structural_variants/signatures/{genome}"
    out = f"{output}/08/structural_variants/vcfs/{genome}"

    os.makedirs(out, exist_ok=True)

    for x in os.listdir(signatures):
        if x.endswith(".bam.svsig.gz"):
            file = os.path.join(signatures, x)
            basename = os.path.basename(file).replace(".bam.svsig.gz", "")
            out_file = f"{out}/{basename}.vcf"

            if os.path.exists(out_file):
                print(f"Skipping {basename}, VCF already exists.")
                continue

            command = f"{smrtlink}/pbsv call {reference} {file} {out_file}"
            print(f"Running: {command}")
            subprocess.run(command, shell=True, check=True)
            
#structural_variants_call(output, reference_t2t, smrtlink, genome_t2t)

def filter_SVs(output, genome):
    SVs_folder = f"{output}/08/structural_variants/vcfs/{genome}"
    out = f"{output}/08/structural_variants/vcfs_filtered/{genome}"

    os.makedirs(out, exist_ok=True)

    for x in os.listdir(SVs_folder):
        if x.endswith("vcf"):
            file = os.path.join(SVs_folder, x)
            basename = os.path.basename(file).replace(".vcf", "")
            filtered_file = f"{out}/{basename}.vcf"

            if os.path.exists(filtered_file):
                print(f"Skipping {basename}, filtered VCF already exists.")
                continue

            program = "../programs/svpack/svpack"
            command = f"{program} filter --pass-only --min-svlen 50 {file} > {filtered_file}"
            print(f"Running: {command}")
            subprocess.run(command, shell=True, check=True)

#filter_SVs(output,genome_t2t) PASS >50pb

##Small variants

def small_variants(output, reference, genome, force=False):
    aligned_reads = os.path.abspath(f"{output}/08/aligned_reads/{genome}")
    out = os.path.abspath(f"{output}/08/small_variants/{genome}")
    os.makedirs(out, exist_ok=True)  # Ensure output directory exists

    reference_name = os.path.basename(reference)
    dir_reference = os.path.abspath(os.path.dirname(reference))
    num_cores = os.cpu_count()

    bam_files = [f for f in os.listdir(aligned_reads) if f.endswith(".bam")]

    for file in bam_files:
        basename = os.path.splitext(file)[0]
        vcf_path = os.path.join(out, f"{basename}.vcf.gz")

        if force or not os.path.exists(vcf_path):
            print(f"Calling variants for {file}...")

            BIN_VERSION = "1.8.0"
            command = [
                "docker", "run", "--rm",
                "-v", f"{aligned_reads}:/input",
                "-v", f"{dir_reference}:/reference",
                "-v", f"{out}:/output",
                f"google/deepvariant:{BIN_VERSION}",
                "/opt/deepvariant/bin/run_deepvariant",
                f"--model_type=PACBIO",
                f"--ref=/reference/{reference_name}",
                f"--reads=/input/{file}",
                f"--output_vcf=/output/{basename}.vcf.gz",
                f"--output_gvcf=/output/{basename}.g.vcf.gz",
                f"--num_shards={num_cores}"
            ]
            subprocess.run(command, check=True)
        else:
            print(f"Skipping {file}, VCF already exists.")

#small_variants(output, reference_t2t, genome_t2t, force=False)
#small_variants(output, reference_hg38, genome_hg38, force=False)

def filter_SNVs(output, genome):
    folder_SNVs=f"{output}/08/small_variants/{genome}"
    out=f"{output}/08/SNVs_filtered/{genome}"
    if not os.path.exists(out):
        os.makedirs(out)
    for x in os.listdir(folder_SNVs):
        if x.endswith(".vcf.gz"):
            file=os.path.join(folder_SNVs, x)
            basename=os.path.basename(file)
            command=f"bcftools view -f 'PASS' {file}  -o {out}/{basename}"
            subprocess.run(command, check=True, shell=True)
            











