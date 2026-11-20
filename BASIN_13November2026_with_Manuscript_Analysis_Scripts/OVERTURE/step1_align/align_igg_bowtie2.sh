#!/usr/bin/env bash
# file: trim_end2end_dm6.sbatch
#
# Purpose: Trim adapters stringently with cutadapt, then align ALL samples end-to-end with Bowtie2 (dm6).
#          Outputs sorted & indexed BAMs (+ flagstat/idxstats). No TLEN splitting here.
#
# Example:
#   sbatch trim_end2end_dm6.sbatch \
#     --proj /work/users/h/s/hsalzler/251208_Make_Universal_Peak_Set_New_IgG_Norm_Strategy/bam_files \
#     --fastq_dir /work/users/h/s/hsalzler/251208_Make_Universal_Peak_Set_New_IgG_Norm_Strategy \
#     --bt2_index_prefix /proj/materalb/users/harmony/Bowtie2_Index_6.38/bowtie2Index/dm6 \
#     --cores 8 --minlen 25 --qcut 20
#

#SBATCH -J trim_end2end_dm6
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o trim_end2end_dm6.%j.out
#SBATCH -e trim_end2end_dm6.%j.err

set -euo pipefail

# ---------- Args ----------
PROJ="BAMS_Trimmed_Aligned_Same_CnR_CnT"
FASTQ_DIR="."
BT2_INDEX_PREFIX="/proj/materalb/users/harmony/Bowtie2_Index_6.38/bowtie2Index/dm6"
CORES="${SLURM_CPUS_PER_TASK:-8}"
MINLEN=25          # cutadapt -m
QCUT=20            # cutadapt -q
ADAPTERS_FILE=""   # optional FASTA with adapters/barcodes to augment defaults

while [[ $# -gt 0 ]]; do
  case "$1" in
    --proj) PROJ="$2"; shift 2;;
    --fastq_dir) FASTQ_DIR="$2"; shift 2;;
    --bt2_index_prefix) BT2_INDEX_PREFIX="$2"; shift 2;;
    --cores) CORES="$2"; shift 2;;
    --minlen) MINLEN="$2"; shift 2;;
    --qcut) QCUT="$2"; shift 2;;
    --adapters_file) ADAPTERS_FILE="$2"; shift 2;;
    -h|--help) sed -n '1,200p' "$0"; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; exit 1;;
  esac
done

[[ -n "$PROJ" && -n "$FASTQ_DIR" && -n "$BT2_INDEX_PREFIX" ]] || {
  echo "[ERROR] --proj, --fastq_dir, --bt2_index_prefix are required." >&2; exit 1; }

# ---------- Modules (Longleaf) ----------
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load cutadapt      2>/dev/null || module load cutadapt/4        2>/dev/null || true
  module load bowtie2       2>/dev/null || module load bowtie2/2.5       2>/dev/null || true
  module load samtools      2>/dev/null || module load samtools/1.17     2>/dev/null || true
fi
command -v cutadapt >/dev/null || { echo "[ERROR] cutadapt not found"; exit 1; }
command -v bowtie2  >/dev/null || { echo "[ERROR] bowtie2 not found"; exit 1; }
command -v samtools >/dev/null || { echo "[ERROR] samtools not found"; exit 1; }

# ---------- Dirs ----------
mkdir -p "${PROJ}/alignment/sam/bowtie2_summary" \
         "${PROJ}/alignment/bam" \
         "${PROJ}/logs" \
         "${PROJ}/trimmed"

# ---------- Adapters (defaults) ----------
TRUSEQ_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
NEXTERA_R1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
NEXTERA_R2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
ME="CTGTCTCTTATACACATCT"   # mosaic end (anywhere)

log(){ echo "[$(date '+%F %T')] $*"; }

is_gz(){ [[ "$1" == *.gz ]]; }

# Derive mate path from R1 common patterns
mate_of_r1(){
  local r1="$1" r2=""
  if   [[ "$r1" =~ _R1_([0-9]{3})\.fastq(\.gz)?$ ]]; then r2="${r1/_R1_/_R2_}"
  elif [[ "$r1" =~ _R1\.fastq(\.gz)?$             ]]; then r2="${r1/_R1./_R2.}"
  elif [[ "$r1" =~ _1\.fastq(\.gz)?$              ]]; then r2="${r1/_1./_2.}"
  else r2="${r1/_R1/_R2}"; r2="${r2/_1/_2}"; fi
  printf "%s" "$r2"
}

# Base sample name from R1 filename
sample_name_from_r1(){
  local base; base="$(basename "$1")"
  base="${base%_R1_001.fastq.gz}"; base="${base%_R1_001.fastq}"
  base="${base%_R1.fastq.gz}";     base="${base%_R1.fastq}"
  base="${base%_1.fastq.gz}";      base="${base%_1.fastq}"
  printf "%s" "$base"
}

# ---------- Collect R1s ----------
shopt -s nullglob
mapfile -t R1S < <(
  find "$FASTQ_DIR" -maxdepth 1 -type f \( \
      -name "*_R1_*.fastq" -o -name "*_R1.fastq" -o -name "*_1.fastq" -o \
      -name "*_R1_*.fastq.gz" -o -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \
    \) | grep -Ev '_I1|_I2' | LC_ALL=C sort
)
(( ${#R1S[@]} > 0 )) || { echo "[ERROR] No R1 FASTQs found in: $FASTQ_DIR"; exit 1; }

log "Project: $PROJ"
log "FASTQ_DIR: $FASTQ_DIR"
log "Bowtie2 index prefix: $BT2_INDEX_PREFIX"
log "Cores: $CORES"
log "Found R1 files: ${#R1S[@]}"

# ---------- Main ----------
for R1 in "${R1S[@]}"; do
  R2="$(mate_of_r1 "$R1")"
  if [[ ! -f "$R2" ]]; then
    log "[WARN] Mate R2 not found for: $R1 (expected: $R2). Skipping."
    continue
  fi

  sample="$(sample_name_from_r1 "$R1")"
  log "Sample=$sample"

  # ---- Trim (stringent, safe if already trimmed) ----
  T1="${PROJ}/trimmed/${sample}_R1.trim.fastq.gz"
  T2="${PROJ}/trimmed/${sample}_R2.trim.fastq.gz"
  if [[ ! -s "$T1" || ! -s "$T2" ]]; then
    log "Cutadapt: q=${QCUT}, m=${MINLEN}, O=7, e=0.1, no-indels"
    if [[ -n "$ADAPTERS_FILE" ]]; then
      cutadapt -j "$CORES" -q "$QCUT" -m "$MINLEN" -O 7 -e 0.10 --no-indels \
        -a file:"$ADAPTERS_FILE" -A file:"$ADAPTERS_FILE" \
        -b file:"$ADAPTERS_FILE" -B file:"$ADAPTERS_FILE" \
        -o "$T1" -p "$T2" "$R1" "$R2" > "${PROJ}/logs/${sample}.cutadapt.log" 2>&1
    else
      cutadapt -j "$CORES" -q "$QCUT" -m "$MINLEN" -O 7 -e 0.10 --no-indels \
        -a "$TRUSEQ_R1"  -A "$TRUSEQ_R2" \
        -a "$NEXTERA_R1" -A "$NEXTERA_R2" \
        -b "$ME"         -B "$ME" \
        -o "$T1" -p "$T2" "$R1" "$R2" > "${PROJ}/logs/${sample}.cutadapt.log" 2>&1
    fi
  else
    log "Trimmed FASTQs exist; reusing"
  fi

  # ---- Align (uniform end-to-end for all) ----
  SAM_OUT="${PROJ}/alignment/sam/${sample}_bowtie2.sam"
  LOG_OUT="${PROJ}/alignment/sam/bowtie2_summary/${sample}_bowtie2.txt"
  if [[ ! -s "$SAM_OUT" ]]; then
    log "Bowtie2 (end-to-end, very-sensitive, no-mixed, no-discordant, -I 10 -X 700)"
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
            -I 10 -X 700 -p "$CORES" -x "$BT2_INDEX_PREFIX" \
            -1 "$T1" -2 "$T2" -S "$SAM_OUT" &> "$LOG_OUT"
  else
    log "SAM exists; skipping align"
  fi

  # ---- BAM sort/index + QC ----
  BAM="${PROJ}/alignment/bam/${sample}.sorted.bam"
  if [[ ! -s "$BAM" ]]; then
    log "SAM→sorted BAM"
    samtools view -@ "$CORES" -bS "$SAM_OUT" | samtools sort -@ "$CORES" -o "$BAM"
    rm -f "${BAM}.bai" "${BAM}.csi"
    samtools index -@ "$CORES" "$BAM"
    samtools flagstat "$BAM" > "${PROJ}/logs/${sample}.flagstat.txt"
    samtools idxstats "$BAM" > "${PROJ}/logs/${sample}.idxstats.txt"
  else
    log "BAM exists; skipping convert/sort"
  fi

  log "Done: $sample"
done

log "All samples complete."
