#!/bin/bash
#SBATCH --job-name=coda
#SBATCH --output=coda_%j.out
#SBATCH --error=coda_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

set -euo pipefail
shopt -s nullglob

# Report where a failure happened. With `set -e`, a non-zero command (e.g. a
# grep that matches nothing, a missing tool) otherwise aborts the job with no
# message in the .out and an empty .err. This trap prints the failing line.
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]##*/} aborted at line ${LINENO} (exit ${rc})." >&2' ERR

module load bedtools
module load deeptools
module load r

########################################
# CONFIG
########################################

# Path to the Pipeline 2 run directory to QC.
# Defaults to the current working directory if not specified as argument 1.
RUN_DIR="${1:-$(pwd)}"

# DOMAIN_BED: BED file of regions of expected signal enrichment for the
# antibody being assessed (e.g., broad domains for H3K27me3, or a peak file
# for a transcription factor). Used to compute domain/interpeak signal ratios
# for QC. Does not need to be perfect — it serves as a reference signal class.
#
# IMPORTANT: If this peak/domain BED was generated empirically from coverage tracks, those
# tracks should use the same canonical geometry as CODA's analytical inputs:
# --centerReads OFF and --extendReads OFF. This keeps peak/domain boundaries and
# derived interpeak regions matched to the signal representation CODA evaluates.
#
# No example domain BED is distributed with BASIN. Supply a BED appropriate
# for your antibody and experiment; set the path below or export DOMAIN_BED
# before running.
#
# If DOMAIN_BED is not an absolute path, the script will also look for it
# inside the bigwigs_baseline/ directory of the run folder.
DOMAIN_BED="${DOMAIN_BED:-/path/to/your/domain_or_peak_file.bed}"

# ANCHOR_BED: The BASIN anchor regions (BARs) file used in Pipeline 2.
# Three universal Drosophila dm6 BAR files are provided in the repo under:
#   HARMONY/templates/anchor_files/dm6/
#
# CRITICAL: this MUST be the SAME file HARMONY used for scaling (its
# universal_bed). CODA's anchor metrics -- including the post-BASIN anchor
# convergence check -- are only meaningful when computed over the same BAR set
# the scale factors came from. Check the HARMONY run's resolved config if unsure.
#
# Choose the appropriate file for your experiment (and match it in HARMONY):
#   ...chr.bed                — full file including all chromosomes (recommended
#                                  default; matches the HARMONY config template)
#   ...NO_HISTONE_GENE_PEAKS.bed — all chromosomes, histone gene peaks removed
#   ...AUTOSOMES_ONLY.bed     — autosomes only; preferable when sex ratio is
#                                  unknown or mixed, BUT only if HARMONY also used it
#
# NOTE: edit the path by replacing the WHOLE quoted value, or `export ANCHOR_BED`
# before running. If you edit only the inner /path/to text, take care not to
# strand the closing `}` of the ${ANCHOR_BED:-...} default onto the filename.
#
# CODA currently supports Drosophila dm6 primary chromosomes only.
# OVERTURE and HARMONY can handle organism-specific BAR files; CODA cannot yet
# evaluate other organisms.
ANCHOR_BED="${ANCHOR_BED:-/path/to/HARMONY/templates/anchor_files/dm6/IgG_consensus_fdr0p005_relaxed_90pct_filter_away_lt50bp.chr.bed}"

# REPEAT_GTF: GTF file of repeat/transposable element annotations for your
# organism, used to exclude repeat regions when defining interpeak intervals.
#
# For Drosophila dm6, two versions are available from the Hammell lab
# TEtranscripts resource folder:
#   https://www.dropbox.com/scl/fo/jdpgn6fl8ngd3th3zebap/ACdZkShDC1au-OckIipI5kM/TEtranscripts/TE_GTF?rlkey=41oz6ppggy82uha5i3yo1rnlx&dl=0
#
#   dm6_BDGP_rmsk_TE.gtf.gz  — BDGP repeat annotations (recommended)
#   dm6_rmsk_TE.gtf.gz        — UCSC RepeatMasker annotations
#
# Both files use non-prefixed chromosome names (2L, 3R, X); CODA normalizes
# these to chr-prefixed names automatically. Both gzipped (.gtf.gz) and
# uncompressed (.gtf) files are accepted.
#
# Files for many other organisms are available at the same Dropbox location.
REPEAT_GTF="${REPEAT_GTF:-/path/to/dm6_BDGP_rmsk_TE.gtf.gz}"

# GENOME_SIZES: Chromosome sizes file for your organism (two-column TSV:
# chromosome name, size in bp). Must use the same chromosome naming convention
# as your BAM files (e.g., chr-prefixed for dm6).
GENOME_SIZES="${GENOME_SIZES:-/path/to/genome.chrom.sizes.txt}"

# PRIMARY_DROP: The primary BASIN drop threshold to use for post-normalization
# QC metrics. Should match one of the drop_fractions values used in Pipeline 2.
# Recommended default: 0.20
PRIMARY_DROP="${PRIMARY_DROP:-0.20}"

# PREFERRED_DROPS: Comma-separated list of drop thresholds to evaluate in QC
# figures. Should be a subset of the drop_fractions used in Pipeline 2.
PREFERRED_DROPS="${PREFERRED_DROPS:-0.10,0.20,0.35,0.50}"

MIN_INTERPEAK_BP="${MIN_INTERPEAK_BP:-100}"
RUN_LABEL="${RUN_LABEL:-postrun_qc}"

########################################
# OPTIONAL MANUAL QC EXPERIMENT GROUPS
########################################
# Use this ONLY when multiple experiments or fractions were run through BASIN
# at the same time and should be QC'd separately. All QC metrics — including
# z-scores, genotype consistency, and experiment-level recommendations — are
# computed within groups, not across the full sample sheet. Samples from
# different experiments, timepoints, or fractions (e.g., supernatant vs pellet)
# should always be placed in separate groups.
#
# Format:
#   group_name:sample1,sample2,sample3;group_name2:sample4,sample5
#   Multiple groups separated by semicolons.
#   Sample names must match sample_id entries in samples.resolved.tsv exactly.
#
# If empty, all samples are placed in one QC group called "all_samples".
# This is appropriate only when the entire sample sheet is a single coherent
# experiment.
#
# If non-empty, ALL samples in the manifest must be listed here or the script
# will fail with an informative error listing the missing sample(s).
#
# Example (replace with your actual sample IDs):
#   MANUAL_QC_GROUPS="exp1:sampleA_rep1,sampleA_rep2,sampleB_rep1;exp2:sampleC_rep1,sampleC_rep2"

MANUAL_QC_GROUPS="${MANUAL_QC_GROUPS:-}"
########################################
# QC THRESHOLDS / DECISION PARAMETERS
########################################
# Defaults are shown in comments. Edit here or export variables before running.
BASELINE_WEIRD_Z="${BASELINE_WEIRD_Z:-2.5}"                    # default 2.5
LOW_ANCHOR_FRAC_MIN="${LOW_ANCHOR_FRAC_MIN:-0.15}"             # default 0.15
LOW_ANCHOR_COUNT_MIN="${LOW_ANCHOR_COUNT_MIN:-50}"             # default 50
INTERPEAK_DIVERSION_Z="${INTERPEAK_DIVERSION_Z:-2}"            # default 2
BASELINE_RATIO_OUTLIER_Z="${BASELINE_RATIO_OUTLIER_Z:-2}"       # default 2
INTERPEAK_DIVERSION_SIGNALS_MIN="${INTERPEAK_DIVERSION_SIGNALS_MIN:-2}" # default 2
INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX="${INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX:-0.25}" # default 0.25; protect low-background groups
INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX="${INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX:-0.30}"   # default 0.30; protect only if domain signal is tight
CLASSWIDE_MIN_SAMPLES="${CLASSWIDE_MIN_SAMPLES:-2}"            # default 2
CLASSWIDE_MIN_FLAGGED_SAMPLES="${CLASSWIDE_MIN_FLAGGED_SAMPLES:-2}" # default 2
CLASSWIDE_MIN_FLAGGED_FRAC="${CLASSWIDE_MIN_FLAGGED_FRAC:-0.5}" # default 0.5
GENOTYPE_SPREAD_RATIO_WARN="${GENOTYPE_SPREAD_RATIO_WARN:-1.5}" # default 1.5
GENOTYPE_SPREAD_CV_WARN="${GENOTYPE_SPREAD_CV_WARN:-0.25}"     # default 0.25
GENOTYPE_LOW_INTERPEAK_PROTECT_MAX="${GENOTYPE_LOW_INTERPEAK_PROTECT_MAX:-0.25}" # default 0.25; protects groups with very low interpeak
GENOTYPE_DOMAIN_RANGE_PROTECT_MAX="${GENOTYPE_DOMAIN_RANGE_PROTECT_MAX:-0.30}"  # default 0.30; only protect if absolute domain signal is tight
PRIMARY_OUTLIER_Z="${PRIMARY_OUTLIER_Z:-2}"                    # default 2
SEX_CHROM_SKEW_Z="${SEX_CHROM_SKEW_Z:-2}"                      # default 2
# Post-BASIN anchor convergence check (on the SCALED bigWigs): per-sample robust
# z for an anchor outlier, and per-group coefficient of variation for whole-group
# non-convergence.
POST_BASIN_ANCHOR_CV_WARN="${POST_BASIN_ANCHOR_CV_WARN:-0.15}" # default 0.15
DROP10_FAIL_FRAC="${DROP10_FAIL_FRAC:-0.25}"                   # default 0.25
DROP50_FAIL_FRAC="${DROP50_FAIL_FRAC:-0.50}"                   # default 0.50
FAIL_INCONSISTENT_GENOTYPE_GROUPS="${FAIL_INCONSISTENT_GENOTYPE_GROUPS:-2}" # default 2
FAIL_MAX_DOMAIN_RATIO_SPREAD="${FAIL_MAX_DOMAIN_RATIO_SPREAD:-2.25}" # default 2.25
FAIL_FRACTION_FLAGGED="${FAIL_FRACTION_FLAGGED:-0.3333333333}" # default 1/3
MANUAL_WARNING_SAMPLE_MIN="${MANUAL_WARNING_SAMPLE_MIN:-3}"    # default 3
WARNING_SAMPLE_MAX="${WARNING_SAMPLE_MAX:-2}"                  # default 2

export BASELINE_WEIRD_Z
export LOW_ANCHOR_FRAC_MIN
export LOW_ANCHOR_COUNT_MIN
export INTERPEAK_DIVERSION_Z
export INTERPEAK_DIVERSION_SIGNALS_MIN
export INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX
export INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX
export CLASSWIDE_MIN_SAMPLES
export CLASSWIDE_MIN_FLAGGED_SAMPLES
export CLASSWIDE_MIN_FLAGGED_FRAC
export GENOTYPE_SPREAD_RATIO_WARN
export GENOTYPE_SPREAD_CV_WARN
export GENOTYPE_LOW_INTERPEAK_PROTECT_MAX
export GENOTYPE_DOMAIN_RANGE_PROTECT_MAX
export PRIMARY_OUTLIER_Z
export SEX_CHROM_SKEW_Z
export POST_BASIN_ANCHOR_CV_WARN
export DROP10_FAIL_FRAC
export DROP50_FAIL_FRAC
export FAIL_INCONSISTENT_GENOTYPE_GROUPS
export FAIL_MAX_DOMAIN_RATIO_SPREAD
export FAIL_FRACTION_FLAGGED
export MANUAL_WARNING_SAMPLE_MIN
export WARNING_SAMPLE_MAX
########################################
# PATHS
########################################

BASELINE_DIR="${RUN_DIR}/bigwigs_baseline"
SCALED_DIR="${RUN_DIR}/bigwigs_scaled"
MANIFEST_TSV="${RUN_DIR}/run_manifest/samples.resolved.tsv"
COVERAGE_CONTRACT_TSV="${RUN_DIR}/run_manifest/coverage_contract.tsv"
SCALES_TSV="${RUN_DIR}/scales/scales_long.tsv"
RETAINED_BARS_TSV="${RUN_DIR}/scales/retained_bars.tsv"
SIGNAL_SUMMARY_TSV="${RUN_DIR}/signal_summary/universal_peaks.rawCounts.tsv"
OUTDIR="${RUN_DIR}/qc_${RUN_LABEL}"

echo "=============================="
echo "RUN_DIR = ${RUN_DIR}"
echo "DOMAIN_BED (initial) = ${DOMAIN_BED}"
echo "BASELINE_DIR = ${BASELINE_DIR}"
echo "SCALED_DIR = ${SCALED_DIR}"
echo "=============================="

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/region_sets"
mkdir -p "${OUTDIR}/per_set"
mkdir -p "${OUTDIR}/manifests"

if [[ ! -f "${MANIFEST_TSV}" ]]; then
  echo "Missing manifest: ${MANIFEST_TSV}"
  exit 1
fi

if [[ ! -f "${SCALES_TSV}" ]]; then
  echo "Missing scales table: ${SCALES_TSV}"
  exit 1
fi

if [[ ! -f "${SIGNAL_SUMMARY_TSV}" ]]; then
  echo "Missing canonical CPM signal summary anchor matrix: ${SIGNAL_SUMMARY_TSV}"
  exit 1
fi

if [[ -f "${COVERAGE_CONTRACT_TSV}" ]]; then
  contract_get() {
    awk -F'\t' -v k="$1" '$1==k {gsub(/\r$/, "", $2); print $2; exit}' "${COVERAGE_CONTRACT_TSV}"
  }

  [[ "$(contract_get canonical_bin_size)" == "1" ]] || { echo "ERROR: CODA requires canonical_bin_size=1." >&2; exit 2; }
  [[ "$(contract_get canonical_center_reads)" == "false" ]] || { echo "ERROR: CODA requires canonical_center_reads=false." >&2; exit 2; }
  [[ "$(contract_get canonical_extend_reads)" == "false" ]] || { echo "ERROR: CODA requires canonical_extend_reads=false." >&2; exit 2; }
  [[ "$(contract_get canonical_qc_norm)" == "CPM" ]] || { echo "ERROR: CODA requires canonical_qc_norm=CPM." >&2; exit 2; }

  echo "Verified HARMONY coverage contract; min_mapq=$(contract_get min_mapq)"
else
  echo "WARNING: legacy HARMONY run has no coverage_contract.tsv; proceeding for backward compatibility." >&2
  echo "         Manually confirm canonical files were 1-bp, uncentered, unextended CPM/NONE." >&2
fi

if [[ -f "${RETAINED_BARS_TSV}" ]]; then
  echo "Using exact HARMONY retained BAR identities: ${RETAINED_BARS_TSV}"
else
  echo "WARNING: retained_bars.tsv not found (legacy run); CODA will reconstruct retained BARs with the corrected HARMONY trim rule." >&2
fi

if [[ ! -f "${DOMAIN_BED}" ]]; then
  if [[ -f "${BASELINE_DIR}/${DOMAIN_BED}" ]]; then
    DOMAIN_BED="${BASELINE_DIR}/${DOMAIN_BED}"
  else
    echo "Missing DOMAIN_BED: ${DOMAIN_BED}"
    exit 1
  fi
fi

echo "FINAL DOMAIN_BED = ${DOMAIN_BED}"

for f in "${ANCHOR_BED}" "${REPEAT_GTF}" "${GENOME_SIZES}"; do
  [[ -f "$f" ]] || { echo "Missing required file: $f"; exit 1; }
done

########################################
# STEP 1: BUILD CLEAN REGION SETS ONCE
########################################

# Normalize a whitespace-delimited file (BED or chrom.sizes) to dm6 primary
# chromosomes. Accepts BOTH chr-prefixed (chr2L) and bare (2L) names and emits
# chr-prefixed output to match the bigWigs (the same convention the repeat GTF
# block below uses). Unlike a bare `grep` under `set -e -o pipefail`, this cannot
# kill the run silently on a zero match: if nothing survives it exits with a
# clear message naming the offending file. ncol=3 for BED, ncol=2 for chrom.sizes.
normalize_dm6_chr() {
  local infile="$1" outfile="$2" label="$3" ncol="$4"
  awk -v ncol="${ncol}" 'BEGIN{OFS="\t"}
    /^[[:space:]]*$/ {next}
    /^#/ {next}
    {
      chr=$1
      if (chr !~ /^chr/) chr="chr" chr
      if (chr=="chr2L"||chr=="chr2R"||chr=="chr3L"||chr=="chr3R"||chr=="chr4"||chr=="chrX"||chr=="chrY") {
        if (ncol==2) print chr,$2
        else         print chr,$2,$3
      }
    }' "${infile}" > "${outfile}"
  if [[ ! -s "${outfile}" ]]; then
    {
      echo "ERROR: no dm6 primary contigs (2L/2R/3L/3R/4/X/Y, chr-prefixed or bare) found in ${label}:"
      echo "         ${infile}"
      echo "       CODA expects dm6 primary chromosomes; an empty result usually means the wrong"
      echo "       file, the wrong genome, or a naming convention CODA does not recognize."
      echo "       First 3 lines of that file:"
      head -3 "${infile}" 2>/dev/null | sed 's/^/         /'
    } >&2
    exit 1
  fi
}

normalize_dm6_chr "${DOMAIN_BED}"   "${OUTDIR}/region_sets/domains.clean.bed"     "DOMAIN_BED"   3
normalize_dm6_chr "${ANCHOR_BED}"   "${OUTDIR}/region_sets/anchors.clean.bed"     "ANCHOR_BED"   3
normalize_dm6_chr "${GENOME_SIZES}" "${OUTDIR}/region_sets/genome.filtered.sizes" "GENOME_SIZES" 2

# Decompress GTF if gzipped, then convert to BED handling both
# chr-prefixed (BDGP: chr2L) and non-prefixed (UCSC: 2L) chromosome names
if [[ "${REPEAT_GTF}" == *.gz ]]; then
  REPEAT_GTF_CMD="zcat ${REPEAT_GTF}"
else
  REPEAT_GTF_CMD="cat ${REPEAT_GTF}"
fi

eval "${REPEAT_GTF_CMD}" | \
awk 'BEGIN{OFS="\t"}
    /^[[:space:]]*$/ {next}
    /^#/ {next}
    {
        chr=$1
        start=$4 - 1
        end=$5
        # Normalize to chr-prefixed names
        # Handles both BDGP (chr2L) and UCSC (2L) GTF conventions
        if (chr !~ /^chr/) {
            if (chr=="2L" || chr=="2R" || chr=="3L" || chr=="3R" || \
                chr=="4"  || chr=="X"  || chr=="Y") {
                chr = "chr" chr
            } else {
                next
            }
        } else {
            if (chr!="chr2L" && chr!="chr2R" && chr!="chr3L" && \
                chr!="chr3R" && chr!="chr4"  && chr!="chrX"  && chr!="chrY") {
                next
            }
        }
        if (start < 0) start = 0
        if (end > start) print chr, start, end
    }' \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> "${OUTDIR}/region_sets/repeats.merged.bed"

cat \
  "${OUTDIR}/region_sets/domains.clean.bed" \
  "${OUTDIR}/region_sets/anchors.clean.bed" \
  "${OUTDIR}/region_sets/repeats.merged.bed" \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> "${OUTDIR}/region_sets/excluded_for_interpeak.merged.bed"

bedtools complement \
  -i "${OUTDIR}/region_sets/excluded_for_interpeak.merged.bed" \
  -g "${OUTDIR}/region_sets/genome.filtered.sizes" \
  > "${OUTDIR}/region_sets/interpeak.raw.bed"

awk -v minbp="${MIN_INTERPEAK_BP}" 'BEGIN{OFS="\t"} ($3-$2) >= minbp {print $1,$2,$3}' \
  "${OUTDIR}/region_sets/interpeak.raw.bed" \
  > "${OUTDIR}/region_sets/interpeak.nonanchor.nonrepeat.min${MIN_INTERPEAK_BP}bp.bed"

########################################
# STEP 2: DISCOVER CANONICAL BIGWIG SETS
########################################

# CODA baseline QC is defined on HARMONY's canonical CPM track only. Optional
# centered/extended derivatives live under bigwigs_downstream/ and are invisible
# to CODA.
find "${BASELINE_DIR}" -maxdepth 1   \( -type f -o -type l \) -name "*.CPM.bw" | sort   > "${OUTDIR}/manifests/baseline_bw_list.txt"

find "${SCALED_DIR}" -maxdepth 1 -mindepth 1 -type d -name "drop*" | sort   > "${OUTDIR}/manifests/drop_dirs.txt"

if [[ ! -s "${OUTDIR}/manifests/baseline_bw_list.txt" ]]; then
  echo "No canonical CPM QC bigWigs (*.CPM.bw) found in ${BASELINE_DIR}." >&2
  echo "CODA no longer accepts RPGC/NONE as interchangeable baseline units." >&2
  exit 1
fi

########################################
# STEP 3: SUMMARIZE ONE BIGWIG SET
########################################

summarize_set() {
  local set_label="$1"
  local bw_list_file="$2"
  local set_outdir="${OUTDIR}/per_set/${set_label}"
  local nproc="${SLURM_CPUS_PER_TASK:-1}"

  mkdir -p "${set_outdir}"

  mapfile -t BWFILES < "${bw_list_file}"
  if [[ ${#BWFILES[@]} -eq 0 ]]; then
    echo "No bigWigs for set ${set_label}; skipping"
    return 0
  fi

  multiBigwigSummary BED-file \
    --BED "${OUTDIR}/region_sets/domains.clean.bed" \
    --bwfiles "${BWFILES[@]}" \
    --outFileName "${set_outdir}/domains.npz" \
    --outRawCounts "${set_outdir}/domains.rawCounts.tsv" \
    --numberOfProcessors "${nproc}"

  multiBigwigSummary BED-file \
    --BED "${OUTDIR}/region_sets/anchors.clean.bed" \
    --bwfiles "${BWFILES[@]}" \
    --outFileName "${set_outdir}/anchors.npz" \
    --outRawCounts "${set_outdir}/anchors.rawCounts.tsv" \
    --numberOfProcessors "${nproc}"

  multiBigwigSummary BED-file \
    --BED "${OUTDIR}/region_sets/interpeak.nonanchor.nonrepeat.min${MIN_INTERPEAK_BP}bp.bed" \
    --bwfiles "${BWFILES[@]}" \
    --outFileName "${set_outdir}/interpeak.npz" \
    --outRawCounts "${set_outdir}/interpeak.rawCounts.tsv" \
    --numberOfProcessors "${nproc}"
}

summarize_set "baseline" "${OUTDIR}/manifests/baseline_bw_list.txt"

if [[ -s "${OUTDIR}/manifests/drop_dirs.txt" ]]; then
  while IFS= read -r dropdir; do
    [[ -n "${dropdir}" ]] || continue
    dropname="$(basename "${dropdir}")"

    find "${dropdir}" -maxdepth 1 \( -type f -o -type l \) \
      -name "*.bw" \
      ! -name "*.gt[0-9]*.bw" \
      ! -name "*.lt[0-9]*.bw" | sort \
      > "${OUTDIR}/manifests/${dropname}_bw_list.txt"

    summarize_set "${dropname}" "${OUTDIR}/manifests/${dropname}_bw_list.txt"
  done < "${OUTDIR}/manifests/drop_dirs.txt"
else
  echo "No drop directories found under ${SCALED_DIR}"
fi

########################################
# STEP 4: WRITE R SCRIPT
########################################

cat > "${OUTDIR}/run_basin_qc.R" <<'EOF'
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
manifest_tsv <- args[2]
scales_tsv <- args[3]
signal_summary_tsv <- args[4]
primary_drop <- as.numeric(args[5])
preferred_drops <- as.numeric(strsplit(args[6], ",")[[1]])
manual_qc_groups <- args[7]
retained_bars_tsv <- args[8]

# Secondary outputs (all per-metric tables and figures) are written into a
# secondary_outputs/ subdirectory. The two primary deliverables --
# sample_qc_summary.tsv and experiment_qc_summary.tsv -- remain at the top of
# the qc_[label]/ directory.
sec_dir <- file.path(outdir, "secondary_outputs")
dir.create(sec_dir, showWarnings = FALSE, recursive = TRUE)

# ========================================================
# QC THRESHOLDS / DECISION PARAMETERS
# ========================================================
# These read the shell variables above. Defaults are shown in comments.
env_num <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || val == "") return(default)
  as.numeric(val)
}
env_int <- function(name, default) as.integer(round(env_num(name, default)))
BASELINE_WEIRD_Z <- env_num("BASELINE_WEIRD_Z", 2.5)                         # default 2.5
LOW_ANCHOR_FRAC_MIN <- env_num("LOW_ANCHOR_FRAC_MIN", 0.15)                  # default 0.15
LOW_ANCHOR_COUNT_MIN <- env_num("LOW_ANCHOR_COUNT_MIN", 50)                  # default 50
INTERPEAK_DIVERSION_Z <- env_num("INTERPEAK_DIVERSION_Z", 2)                 # default 2
BASELINE_RATIO_OUTLIER_Z <- env_num("BASELINE_RATIO_OUTLIER_Z", 2)            # default 2
INTERPEAK_DIVERSION_SIGNALS_MIN <- env_int("INTERPEAK_DIVERSION_SIGNALS_MIN", 2) # default 2
INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX <- env_num("INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX", 0.25) # default 0.25
INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX <- env_num("INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX", 0.30) # default 0.30
CLASSWIDE_MIN_SAMPLES <- env_int("CLASSWIDE_MIN_SAMPLES", 2)                 # default 2
CLASSWIDE_MIN_FLAGGED_SAMPLES <- env_int("CLASSWIDE_MIN_FLAGGED_SAMPLES", 2) # default 2
CLASSWIDE_MIN_FLAGGED_FRAC <- env_num("CLASSWIDE_MIN_FLAGGED_FRAC", 0.5)     # default 0.5
GENOTYPE_SPREAD_RATIO_WARN <- env_num("GENOTYPE_SPREAD_RATIO_WARN", 1.5)     # default 1.5
GENOTYPE_SPREAD_CV_WARN <- env_num("GENOTYPE_SPREAD_CV_WARN", 0.25)          # default 0.25
GENOTYPE_LOW_INTERPEAK_PROTECT_MAX <- env_num("GENOTYPE_LOW_INTERPEAK_PROTECT_MAX", 0.25) # default 0.25
GENOTYPE_DOMAIN_RANGE_PROTECT_MAX <- env_num("GENOTYPE_DOMAIN_RANGE_PROTECT_MAX", 0.30)   # default 0.30
PRIMARY_OUTLIER_Z <- env_num("PRIMARY_OUTLIER_Z", 2)                         # default 2
SEX_CHROM_SKEW_Z <- env_num("SEX_CHROM_SKEW_Z", 2)                           # default 2
# Post-BASIN pipeline integrity check (computed on the SCALED bigWigs). After
# correct BASIN scaling, every sample's trimmed anchor median should be ~1.0
# by construction. POST_BASIN_ANCHOR_CV_WARN flags an experiment group whose
# scaled anchor medians fail to converge, indicating a mechanical problem
# (wrong anchor file, scale factor not applied) rather than biology.
POST_BASIN_ANCHOR_CV_WARN <- env_num("POST_BASIN_ANCHOR_CV_WARN", 0.15)      # default 0.15
DROP10_FAIL_FRAC <- env_num("DROP10_FAIL_FRAC", 0.25)                        # default 0.25
DROP50_FAIL_FRAC <- env_num("DROP50_FAIL_FRAC", 0.50)                        # default 0.50
FAIL_INCONSISTENT_GENOTYPE_GROUPS <- env_int("FAIL_INCONSISTENT_GENOTYPE_GROUPS", 2) # default 2
FAIL_MAX_DOMAIN_RATIO_SPREAD <- env_num("FAIL_MAX_DOMAIN_RATIO_SPREAD", 2.25) # default 2.25
FAIL_FRACTION_FLAGGED <- env_num("FAIL_FRACTION_FLAGGED", 1/3)               # default 1/3
MANUAL_WARNING_SAMPLE_MIN <- env_int("MANUAL_WARNING_SAMPLE_MIN", 3)         # default 3
WARNING_SAMPLE_MAX <- env_int("WARNING_SAMPLE_MAX", 2)                       # default 2

normalize_sample_name <- function(x) {
  x %>%
    basename() %>%
    str_replace_all("^['\"]+|['\"]+$", "") %>%
    str_trim() %>%
    # Scaled bigWigs can come from either duplicate-removed or duplicate-kept runs.
    # Strip the BASIN/drop/dup-handling suffix back to the manifest sample_id.
    # Examples handled:
    #   sample.scaled.drop0.20.dedup.bw
    #   sample.scaled.drop0.20.keepdups.bw
    #   sample.scaled.drop0.20.keep_dups.bw
    #   sample.scaled.drop0.20.keep-dup.bw
    str_replace("\\.scaled\\.drop[0-9.]+\\.(dedup|keepdups|keep_dups|keep-?dup|dups?kept)\\.bw$", "") %>%
    # Baseline bigWigs may also carry duplicate-handling suffixes before .CPM.bw.
    str_replace("\\.(dedup|keepdups|keep_dups|keep-?dup|dups?kept)\\.CPM\\.bw$", "") %>%
    str_replace("\\.CPM\\.bw$", "") %>%
    # Scale factors are summarized from the raw (NONE) baseline bigWigs, whose
    # columns in scales_long.tsv are {sample_id}.NONE.bw. Strip that too so the
    # scale table joins to the scaled-bigwig sample ids.
    str_replace("\\.(dedup|keepdups|keep_dups|keep-?dup|dups?kept)\\.NONE\\.bw$", "") %>%
    str_replace("\\.NONE\\.bw$", "")
}

drop_all_zero_rows <- function(df) {
  signal_cols <- setdiff(colnames(df), c("chr", "start", "end"))
  df %>%
    filter(rowSums(across(all_of(signal_cols), ~ . != 0)) > 0)
}

has_exact_retained_bars <- !is.na(retained_bars_tsv) && retained_bars_tsv != "" && file.exists(retained_bars_tsv)
retained_bars_exact <- NULL
if (has_exact_retained_bars) {
  retained_bars_exact <- read_tsv(retained_bars_tsv, show_col_types = FALSE) %>%
    transmute(
      sample = normalize_sample_name(sample_colname),
      threshold = as.numeric(drop_fraction),
      chr = as.character(chr),
      chr = ifelse(str_starts(chr, "chr"), chr, paste0("chr", chr)),
      start = as.character(start),
      end = as.character(end),
      region_id = paste(chr, start, end, sep = ":"),
      raw_signal = suppressWarnings(as.numeric(raw_signal))
    )
}

safe_mad <- function(x) {
  out <- mad(x, na.rm = TRUE)
  ifelse(is.na(out) | out == 0, NA_real_, out)
}

parse_manual_qc_groups <- function(x) {
  if (is.na(x) || x == "") {
    return(tibble(sample = character(), qc_group = character()))
  }

  group_entries <- strsplit(x, ";", fixed = TRUE)[[1]]
  group_entries <- trimws(group_entries)
  group_entries <- group_entries[group_entries != ""]

  parsed <- purrr::map_dfr(group_entries, function(entry) {
    parts <- strsplit(entry, ":", fixed = TRUE)[[1]]

    if (length(parts) != 2) {
      stop(paste0("Bad MANUAL_QC_GROUPS entry: ", entry))
    }

    group_name <- trimws(parts[1])
    if (group_name == "") {
      stop(paste0("Blank QC group name in entry: ", entry))
    }

    sample_names <- strsplit(parts[2], ",", fixed = TRUE)[[1]]
    sample_names <- trimws(sample_names)
    sample_names <- sample_names[sample_names != ""]

    if (length(sample_names) == 0) {
      stop(paste0("QC group has no samples listed: ", entry))
    }

    tibble(
      sample = sample_names,
      qc_group = group_name
    )
  })

  dup_samples <- parsed %>%
    count(sample) %>%
    filter(n > 1)

  if (nrow(dup_samples) > 0) {
    stop(
      paste0(
        "Samples assigned to more than one QC group: ",
        paste(dup_samples$sample, collapse = ", ")
      )
    )
  }

  parsed
}

read_mbs_long <- function(file, region_class, set_label, drop_zero_rows = FALSE) {
  df <- read_tsv(file, show_col_types = FALSE)
  colnames(df)[1:3] <- c("chr", "start", "end")

  if (drop_zero_rows) {
    df <- drop_all_zero_rows(df)
  }

  sample_cols <- setdiff(colnames(df), c("chr", "start", "end"))

  df %>%
    mutate(region_id = paste(chr, start, end, sep = ":")) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample_file", values_to = "signal") %>%
    mutate(
      sample = normalize_sample_name(sample_file),
      region_class = region_class,
      set_label = set_label
    )
}

compute_retained_anchor_chrom_fractions <- function(file, set_label, drop_zero_rows = TRUE) {
  df <- read_tsv(file, show_col_types = FALSE)
  colnames(df)[1:3] <- c("chr", "start", "end")

  if (drop_zero_rows) {
    df <- drop_all_zero_rows(df)
  }

  drop_fraction <- ifelse(set_label == "baseline", 0, as.numeric(str_extract(set_label, "[0-9.]+")))
  sample_cols <- setdiff(colnames(df), c("chr", "start", "end"))

  long_df <- df %>%
    mutate(anchor_row_id = row_number()) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample_file", values_to = "signal") %>%
    mutate(
      sample = normalize_sample_name(sample_file),
      chrom_group = case_when(
        chr %in% c("chrX", "X") ~ "X",
        chr %in% c("chrY", "Y") ~ "Y",
        chr %in% c("chr2L", "2L", "chr2R", "2R",
                   "chr3L", "3L", "chr3R", "3R",
                   "chr4", "4") ~ "autosomes",
        TRUE ~ "other"
      )
    ) %>%
    filter(chrom_group != "other")

  retained_df <- long_df %>%
    group_by(sample) %>%
    arrange(desc(signal), .by_group = TRUE) %>%
    mutate(
      n_anchor_rows = n(),
      n_trim = round(n_anchor_rows * drop_fraction),
      rank_desc_signal = row_number(),
      retained = rank_desc_signal > n_trim
    ) %>%
    ungroup() %>%
    filter(retained)

  chrom_counts <- retained_df %>%
    group_by(sample, chrom_group) %>%
    summarise(n_retained_anchor_rows = n(), .groups = "drop") %>%
    pivot_wider(names_from = chrom_group, values_from = n_retained_anchor_rows, values_fill = 0)

  if (!"autosomes" %in% colnames(chrom_counts)) chrom_counts$autosomes <- 0
  if (!"X" %in% colnames(chrom_counts)) chrom_counts$X <- 0
  if (!"Y" %in% colnames(chrom_counts)) chrom_counts$Y <- 0

  chrom_counts %>%
    mutate(
      total_retained_anchor_rows = autosomes + X + Y,
      frac_retained_anchor_rows_autosomes = ifelse(total_retained_anchor_rows > 0, autosomes / total_retained_anchor_rows, NA_real_),
      frac_retained_anchor_rows_X = ifelse(total_retained_anchor_rows > 0, X / total_retained_anchor_rows, NA_real_),
      frac_retained_anchor_rows_Y = ifelse(total_retained_anchor_rows > 0, Y / total_retained_anchor_rows, NA_real_),
      sex_chrom_retained_anchor_fraction = frac_retained_anchor_rows_X + frac_retained_anchor_rows_Y,
      set_label = set_label,
      threshold = ifelse(set_label == "baseline", NA_real_, drop_fraction)
    )
}

compute_retained_anchor_chrom_fractions_exact <- function(set_label) {
  drop_fraction <- as.numeric(str_extract(set_label, "[0-9.]+"))
  retained_bars_exact %>%
    filter(abs(threshold - drop_fraction) < 1e-8) %>%
    mutate(
      chrom_group = case_when(
        chr %in% c("chrX", "X") ~ "X",
        chr %in% c("chrY", "Y") ~ "Y",
        chr %in% c("chr2L", "2L", "chr2R", "2R", "chr3L", "3L", "chr3R", "3R", "chr4", "4") ~ "autosomes",
        TRUE ~ "other"
      )
    ) %>%
    filter(chrom_group != "other") %>%
    count(sample, chrom_group, name = "n_retained_anchor_rows") %>%
    complete(sample, chrom_group = c("autosomes", "X", "Y"), fill = list(n_retained_anchor_rows = 0L)) %>%
    pivot_wider(names_from = chrom_group, values_from = n_retained_anchor_rows, values_fill = 0) %>%
    mutate(
      total_retained_anchor_rows = autosomes + X + Y,
      frac_retained_anchor_rows_autosomes = ifelse(total_retained_anchor_rows > 0, autosomes / total_retained_anchor_rows, NA_real_),
      frac_retained_anchor_rows_X = ifelse(total_retained_anchor_rows > 0, X / total_retained_anchor_rows, NA_real_),
      frac_retained_anchor_rows_Y = ifelse(total_retained_anchor_rows > 0, Y / total_retained_anchor_rows, NA_real_),
      sex_chrom_retained_anchor_fraction = frac_retained_anchor_rows_X + frac_retained_anchor_rows_Y,
      set_label = set_label,
      threshold = drop_fraction
    )
}

read_one_set <- function(set_dir, set_label, signal_summary_tsv) {
  domain_long <- read_mbs_long(file.path(set_dir, "domains.rawCounts.tsv"), "domain", set_label, drop_zero_rows = FALSE)
  interpeak_long <- read_mbs_long(file.path(set_dir, "interpeak.rawCounts.tsv"), "interpeak", set_label, drop_zero_rows = FALSE)

  if (set_label == "baseline") {
    anchor_long <- read_mbs_long(signal_summary_tsv, "anchor", set_label, drop_zero_rows = TRUE)
  } else {
    anchor_long <- read_mbs_long(file.path(set_dir, "anchors.rawCounts.tsv"), "anchor", set_label, drop_zero_rows = TRUE)
  }

  bind_rows(domain_long, anchor_long, interpeak_long)
}

set_dirs <- list.dirs(file.path(outdir, "per_set"), recursive = FALSE, full.names = TRUE)
set_labels <- basename(set_dirs)

all_long <- map2_dfr(set_dirs, set_labels, ~ read_one_set(.x, .y, signal_summary_tsv))

retained_anchor_chrom_fractions <- map2_dfr(
  set_dirs,
  set_labels,
  function(set_dir, set_label) {
    if (set_label == "baseline") {
      compute_retained_anchor_chrom_fractions(signal_summary_tsv, set_label, drop_zero_rows = TRUE)
    } else if (has_exact_retained_bars) {
      compute_retained_anchor_chrom_fractions_exact(set_label)
    } else {
      compute_retained_anchor_chrom_fractions(file.path(set_dir, "anchors.rawCounts.tsv"), set_label, drop_zero_rows = TRUE)
    }
  }
)

write_tsv(retained_anchor_chrom_fractions,
          file.path(sec_dir, "retained_anchor_chromosome_fraction_by_threshold.tsv"))

sample_stats <- all_long %>%
  group_by(set_label, sample, region_class) %>%
  reframe(
    n_regions = n(),
    n_nonzero_rows = sum(signal != 0, na.rm = TRUE),
    frac_nonzero_rows = mean(signal != 0, na.rm = TRUE),
    frac_zero = mean(signal == 0, na.rm = TRUE),
    mean_signal = mean(signal, na.rm = TRUE),
    median_signal = median(signal, na.rm = TRUE),
    p90_signal = as.numeric(quantile(signal, 0.90, na.rm = TRUE))
  )

write_tsv(sample_stats, file.path(sec_dir, "sample_region_class_stats_all_sets.tsv"))

baseline_region_class_stats <- sample_stats %>%
  filter(set_label == "baseline")

write_tsv(baseline_region_class_stats,
          file.path(sec_dir, "baseline_region_class_stats.tsv"))

wide_median <- sample_stats %>%
  select(set_label, sample, region_class, median_signal) %>%
  pivot_wider(names_from = region_class, values_from = median_signal, names_prefix = "median_")

wide_mean <- sample_stats %>%
  select(set_label, sample, region_class, mean_signal) %>%
  pivot_wider(names_from = region_class, values_from = mean_signal, names_prefix = "mean_")

wide_support <- sample_stats %>%
  filter(region_class == "anchor") %>%
  select(set_label, sample, n_nonzero_anchor_rows = n_nonzero_rows, frac_nonzero_anchor_rows = frac_nonzero_rows)

wide_retained_chrom <- retained_anchor_chrom_fractions %>%
  select(set_label, sample, threshold,
         total_retained_anchor_rows,
         frac_retained_anchor_rows_autosomes,
         frac_retained_anchor_rows_X,
         frac_retained_anchor_rows_Y,
         sex_chrom_retained_anchor_fraction)

wide <- wide_median %>%
  left_join(wide_mean, by = c("set_label", "sample")) %>%
  left_join(wide_support, by = c("set_label", "sample")) %>%
  mutate(
    threshold = case_when(
      set_label == "baseline" ~ NA_real_,
      TRUE ~ as.numeric(str_extract(set_label, "[0-9.]+"))
    ),
    set_type = ifelse(set_label == "baseline", "baseline", "scaled")
  ) %>%
  left_join(wide_retained_chrom, by = c("set_label", "sample", "threshold")) %>%
  mutate(
    domain_interpeak = median_domain / ifelse(median_interpeak == 0, NA, median_interpeak),
    anchor_interpeak = median_anchor / ifelse(median_interpeak == 0, NA, median_interpeak),
    anchor_domain = median_anchor / ifelse(median_domain == 0, NA, median_domain),
    mean_anchor_interpeak = mean_anchor / ifelse(mean_interpeak == 0, NA, mean_interpeak)
  )

baseline_ratio_table <- wide %>%
  filter(set_type == "baseline") %>%
  select(
    sample,
    threshold,
    median_domain,
    median_anchor,
    median_interpeak,
    mean_anchor,
    mean_interpeak,
    n_nonzero_anchor_rows,
    frac_nonzero_anchor_rows,
    domain_interpeak,
    anchor_interpeak,
    anchor_domain,
    mean_anchor_interpeak
  ) %>%
  arrange(sample)

write_tsv(baseline_ratio_table,
          file.path(sec_dir, "baseline_ratio_table.tsv"))

manual_qc_group_df <- parse_manual_qc_groups(manual_qc_groups)

manifest <- read_tsv(manifest_tsv, show_col_types = FALSE) %>%
  mutate(sample = sample_id) %>%
  distinct(sample, .keep_all = TRUE)

if (nrow(manual_qc_group_df) == 0) {
  manifest <- manifest %>%
    mutate(qc_group = "all_samples")
} else {
  extra_grouped_samples <- setdiff(manual_qc_group_df$sample, manifest$sample)
  if (length(extra_grouped_samples) > 0) {
    stop(
      paste0(
        "sample(s) listed in MANUAL_QC_GROUPS but absent from manifest: ",
        paste(extra_grouped_samples, collapse = ", ")
      )
    )
  }

  manifest <- manifest %>%
    left_join(manual_qc_group_df, by = "sample")

  missing_samples <- manifest %>%
    filter(is.na(qc_group) | qc_group == "") %>%
    pull(sample)

  if (length(missing_samples) > 0) {
    message(
      paste0(
        "Ignoring sample(s) not listed in MANUAL_QC_GROUPS: ",
        paste(missing_samples, collapse = ", ")
      )
    )

    manifest <- manifest %>%
      filter(!is.na(qc_group) & qc_group != "")
  }
}

# Preserve any original manifest metadata columns, including fraction if present.
# IMPORTANT: fraction is NOT used to define QC groups.
# QC grouping is explicit: MANUAL_QC_GROUPS, or all_samples if MANUAL_QC_GROUPS is empty.
if (!"fraction" %in% colnames(manifest)) {
  manifest$fraction <- NA_character_
}
manifest <- manifest %>%
  mutate(
    experiment_group = paste(antibody, qc_group, sep = "__"),
    genotype_group = paste(antibody, qc_group, genotype, sep = "__")
  )

scales <- read_tsv(scales_tsv, show_col_types = FALSE) %>%
  transmute(
    sample = normalize_sample_name(sample_colname),
    threshold = as.numeric(drop_fraction),
    scale_factor
  ) %>%
  distinct(sample, threshold, .keep_all = TRUE)

########################################################
# DROP-FRACTION STABILITY (advisory)
# ------------------------------------------------------
# Helps answer "which drop fraction is legitimate?" WITHOUT relying on
# post-scaling anchor convergence (which is partly circular, since scaling
# forces the anchor median toward a reference at every drop).
#
#   (1) Cross-replicate agreement of DEPTH-NORMALIZED (CPM) medBAR at each drop.
#       For replicates of the same genotype the background estimate should agree;
#       the drop that minimizes the within-genotype coefficient of variation (CV)
#       of CPM medBAR gives the most reproducible background. Non-circular
#       because medBAR is measured before BASIN scaling. CPM medBAR (not raw) is
#       used so that sequencing-depth differences between replicates do not
#       inflate the CV -- this is the same principle behind all the baseline QC:
#       depth must be factored out before comparing samples.
#   (2) Elbow: per-sample relative change in the scale factor between consecutive
#       drops. Where the curve flattens, additional trimming is removing genuine
#       background rather than antibody-enriched anchors.
#
# Both are ADVISORY. The authoritative choice still requires a positive control
# (e.g. a known titration) to confirm that genotype differences survive while
# replicate noise collapses. Outputs land in secondary_outputs/.
########################################################

# Compute depth-normalized (CPM) medBAR per sample per drop from the CPM signal
# summary. This replicates HARMONY's trim logic (sort BARs descending, drop the
# top fraction, take the median) but on CPM values, so the result is comparable
# across samples regardless of sequencing depth.
ds_cpm_anchor <- read_tsv(signal_summary_tsv, show_col_types = FALSE)
colnames(ds_cpm_anchor)[1:3] <- c("chr", "start", "end")
ds_cpm_anchor <- drop_all_zero_rows(ds_cpm_anchor)
ds_cpm_sample_cols <- setdiff(colnames(ds_cpm_anchor), c("chr", "start", "end"))
ds_cpm_long <- ds_cpm_anchor %>%
  pivot_longer(cols = all_of(ds_cpm_sample_cols), names_to = "sample_file", values_to = "cpm_signal") %>%
  mutate(sample = normalize_sample_name(sample_file)) %>%
  filter(!is.na(cpm_signal)) %>%
  group_by(sample) %>%
  filter(any(cpm_signal != 0)) %>%
  ungroup()

# Scale factors and drop fractions from scales_long.tsv (for the elbow signal).
ds_sf <- read_tsv(scales_tsv, show_col_types = FALSE) %>%
  transmute(
    sample = normalize_sample_name(sample_colname),
    threshold = as.numeric(drop_fraction),
    scale_factor = suppressWarnings(as.numeric(scale_factor))
  ) %>%
  distinct(sample, threshold, .keep_all = TRUE)

# For each sample x drop, trim the top fraction of BARs by CPM signal and take
# the median of the kept set. This is medBAR in depth-normalized units.
ds_drop_fracs <- sort(unique(ds_sf$threshold))
ds_cpm_medbar <- map_dfr(ds_drop_fracs, function(d) {
  ds_cpm_long %>%
    group_by(sample) %>%
    arrange(desc(cpm_signal), .by_group = TRUE) %>%
    mutate(
      n_bars = n(),
      rank_desc = row_number(),
      n_trim = round(n_bars * d)
    ) %>%
    filter(rank_desc > n_trim) %>%
    summarise(med_bar = median(cpm_signal, na.rm = TRUE), .groups = "drop") %>%
    mutate(threshold = d)
})

ds_scales <- ds_sf %>%
  left_join(ds_cpm_medbar, by = c("sample", "threshold")) %>%
  left_join(
    manifest %>% select(sample, qc_group, experiment_group, genotype_group, genotype),
    by = "sample"
  )

# (1) Non-circular replicate agreement: CV of depth-normalized (CPM) medBAR across replicates of the
#     same genotype at each drop (n >= 3 required for a meaningful CV).
ds_by_genotype <- ds_scales %>%
  filter(!is.na(med_bar)) %>%
  group_by(experiment_group, genotype_group, genotype, threshold) %>%
  summarise(
    n_reps = n(),
    medbar_mean = mean(med_bar, na.rm = TRUE),
    medbar_sd = sd(med_bar, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    medbar_cv = ifelse(n_reps >= 3 & medbar_mean != 0, medbar_sd / medbar_mean, NA_real_)
  ) %>%
  arrange(experiment_group, genotype_group, threshold)
write_tsv(ds_by_genotype, file.path(sec_dir, "drop_stability_by_genotype.tsv"))

# (2) Elbow: per-sample relative step in scale factor between consecutive drops.
ds_elbow <- ds_scales %>%
  filter(!is.na(scale_factor)) %>%
  arrange(sample, threshold) %>%
  group_by(sample) %>%
  mutate(
    prev_sf = lag(scale_factor),
    sf_rel_step = ifelse(!is.na(prev_sf) & prev_sf != 0,
                         abs(scale_factor - prev_sf) / prev_sf, NA_real_)
  ) %>%
  ungroup()

# (3) Per-experiment, per-drop summary combining both signals.
ds_by_experiment <- ds_by_genotype %>%
  group_by(experiment_group, threshold) %>%
  summarise(
    n_genotype_groups_evaluated = sum(!is.na(medbar_cv)),
    mean_within_genotype_medbar_cv = ifelse(any(!is.na(medbar_cv)),
                                            mean(medbar_cv, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  left_join(
    ds_elbow %>%
      group_by(experiment_group, threshold) %>%
      summarise(
        n_samples = n(),
        median_scalefactor_rel_step = ifelse(any(!is.na(sf_rel_step)),
                                             median(sf_rel_step, na.rm = TRUE), NA_real_),
        .groups = "drop"
      ),
    by = c("experiment_group", "threshold")
  ) %>%
  arrange(experiment_group, threshold)
write_tsv(ds_by_experiment, file.path(sec_dir, "drop_stability_by_experiment.tsv"))

# (4) Advisory recommendation: the drop minimizing cross-replicate medBAR CV.
ds_recommendation <- ds_by_experiment %>%
  filter(!is.na(mean_within_genotype_medbar_cv)) %>%
  group_by(experiment_group) %>%
  arrange(mean_within_genotype_medbar_cv, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    experiment_group,
    suggested_drop_by_medbar_cv = threshold,
    medbar_cv_at_suggested = mean_within_genotype_medbar_cv,
    median_scalefactor_rel_step_at_suggested = median_scalefactor_rel_step,
    primary_drop = primary_drop,
    note = paste0(
      "Advisory only: drop minimizing cross-replicate CV of depth-normalized (CPM) medBAR ",
      "(most reproducible background estimate). Confirm against a positive ",
      "control before locking quantitative cross-sample claims."
    )
  )
write_tsv(ds_recommendation, file.path(sec_dir, "drop_recommendation.tsv"))

########################################################
# REPLICATE COUNT CHECK
# Compute per-sample replicate counts within each
# genotype group and experiment group. Used to suppress
# z-score metrics and issue warnings when n < 3.
#
# n=1: z-scores undefined (MAD=0, division by zero).
#      Flags all z-score metrics as NA.
#      Issues manual_inspection_required at experiment level.
# n=2: z-scores always +/-0.67 regardless of true signal
#      difference. Flags all z-score metrics as NA.
#      Issues warning at experiment level.
# n>=3: Normal operation. No flag.
########################################################

replicate_counts <- manifest %>%
  group_by(genotype_group, experiment_group) %>%
  summarise(
    n_reps_in_genotype_group = n(),
    .groups = "drop"
  )

experiment_rep_counts <- manifest %>%
  group_by(experiment_group) %>%
  summarise(
    n_samples_in_experiment = n(),
    .groups = "drop"
  )

manifest <- manifest %>%
  left_join(replicate_counts, by = c("genotype_group", "experiment_group")) %>%
  left_join(experiment_rep_counts, by = "experiment_group") %>%
  mutate(
    low_replicate_count = n_reps_in_genotype_group < 3,
    very_low_replicate_count = n_reps_in_genotype_group < 2,
    low_replicate_note = case_when(
      n_reps_in_genotype_group == 1 ~
        paste0("only 1 sample in genotype group '",
               genotype, "' — z-score and consistency metrics are undefined and have been suppressed"),
      n_reps_in_genotype_group == 2 ~
        paste0("only 2 samples in genotype group '",
               genotype, "' — z-score metrics always +/-0.67 regardless of signal; suppressed as uninformative"),
      TRUE ~ NA_character_
    ),
    experiment_low_replicate_note = case_when(
      n_samples_in_experiment == 1 ~
        "experiment has only 1 sample — CODA's statistical framework requires at least 3 replicates; most metrics suppressed",
      n_samples_in_experiment == 2 ~
        "experiment has only 2 samples — z-score metrics are uninformative with n=2 and have been suppressed",
      TRUE ~ NA_character_
    )
  )

qc <- wide %>%
  distinct(set_label, sample, .keep_all = TRUE) %>%
  left_join(manifest, by = "sample") %>%
  left_join(scales, by = c("sample", "threshold"))

scaled_metrics <- qc %>%
  filter(set_type == "scaled")

cat("\nSamples in wide but not manifest:\n")
print(setdiff(unique(wide$sample), unique(manifest$sample)))

cat("\nSamples in manifest but not wide:\n")
print(setdiff(unique(manifest$sample), unique(wide$sample)))

cat("\nSamples in scaled metrics but not scales table:\n")
print(setdiff(unique(na.omit(scaled_metrics$sample)), unique(scales$sample)))

baseline_metrics <- qc %>%
  filter(set_type == "baseline") %>%
  group_by(experiment_group) %>%
  mutate(
    baseline_interpeak_z = as.numeric(scale(median_interpeak)),
    baseline_anchor_interpeak_z = as.numeric(scale(mean_anchor_interpeak))
  ) %>%
  ungroup() %>%
  select(
    sample, qc_group, experiment_group, genotype_group, antibody, genotype, replicate,
    n_reps_in_genotype_group, low_replicate_count, very_low_replicate_count,
    low_replicate_note, experiment_low_replicate_note,
    baseline_median_domain = median_domain,
    baseline_median_anchor = median_anchor,
    baseline_median_interpeak = median_interpeak,
    baseline_mean_anchor = mean_anchor,
    baseline_mean_interpeak = mean_interpeak,
    baseline_domain_interpeak = domain_interpeak,
    baseline_anchor_interpeak = anchor_interpeak,
    baseline_anchor_domain = anchor_domain,
    baseline_mean_anchor_interpeak = mean_anchor_interpeak,
    baseline_interpeak_z,
    baseline_anchor_interpeak_z,
    baseline_n_nonzero_anchor_rows = n_nonzero_anchor_rows,
    baseline_frac_nonzero_anchor_rows = frac_nonzero_anchor_rows
  ) %>%
  mutate(
    # Suppress z-scores when n < 3 — they are undefined (n=1) or always
    # +/-0.67 regardless of true signal difference (n=2)
    baseline_interpeak_z = ifelse(low_replicate_count, NA_real_, baseline_interpeak_z),
    baseline_anchor_interpeak_z = ifelse(low_replicate_count, NA_real_, baseline_anchor_interpeak_z),
    baseline_weird = case_when(
      low_replicate_count ~ NA,
      TRUE ~ abs(baseline_interpeak_z) > BASELINE_WEIRD_Z | abs(baseline_anchor_interpeak_z) > BASELINE_WEIRD_Z
    ),
    low_anchor_support = baseline_frac_nonzero_anchor_rows < LOW_ANCHOR_FRAC_MIN | baseline_n_nonzero_anchor_rows < LOW_ANCHOR_COUNT_MIN
  )

########################################################
# NEW: INTERPEAK-DIVERSION QC
########################################################

interpeak_diversion_metrics <- baseline_metrics %>%
  group_by(experiment_group) %>%
  mutate(
    exp_n = n(),

    exp_median_interpeak = median(baseline_median_interpeak, na.rm = TRUE),
    exp_domain_range = diff(range(baseline_median_domain, na.rm = TRUE)),

    med_interpeak = median(baseline_median_interpeak, na.rm = TRUE),
    mad_interpeak = safe_mad(baseline_median_interpeak),

    med_domain_interpeak = median(baseline_domain_interpeak, na.rm = TRUE),
    mad_domain_interpeak = safe_mad(baseline_domain_interpeak),

    med_anchor_domain = median(baseline_anchor_domain, na.rm = TRUE),
    mad_anchor_domain = safe_mad(baseline_anchor_domain),

    # Suppress z-scores when any sample in the genotype group has n < 3
    interpeak_high_robust_z = ifelse(
      low_replicate_count | is.na(mad_interpeak),
      NA_real_,
      (baseline_median_interpeak - med_interpeak) / mad_interpeak
    ),

    domain_interpeak_low_robust_z = ifelse(
      low_replicate_count | is.na(mad_domain_interpeak),
      NA_real_,
      (baseline_domain_interpeak - med_domain_interpeak) / mad_domain_interpeak
    ),

    anchor_domain_high_robust_z = ifelse(
      low_replicate_count | is.na(mad_anchor_domain),
      NA_real_,
      (baseline_anchor_domain - med_anchor_domain) / mad_anchor_domain
    )
  ) %>%
  ungroup() %>%
  mutate(
    interpeak_high_flag = case_when(
      low_replicate_count ~ NA,
      !is.na(interpeak_high_robust_z) & interpeak_high_robust_z > INTERPEAK_DIVERSION_Z ~ TRUE,
      TRUE ~ FALSE
    ),

    domain_interpeak_low_flag = case_when(
      low_replicate_count ~ NA,
      !is.na(domain_interpeak_low_robust_z) & domain_interpeak_low_robust_z < -INTERPEAK_DIVERSION_Z ~ TRUE,
      TRUE ~ FALSE
    ),

    baseline_domain_interpeak_outlier = case_when(
      low_replicate_count ~ NA,
      !is.na(domain_interpeak_low_robust_z) &
        abs(domain_interpeak_low_robust_z) > BASELINE_RATIO_OUTLIER_Z ~ TRUE,
      TRUE ~ FALSE
    ),

    anchor_domain_high_flag = case_when(
      low_replicate_count ~ NA,
      !is.na(anchor_domain_high_robust_z) & anchor_domain_high_robust_z > INTERPEAK_DIVERSION_Z ~ TRUE,
      TRUE ~ FALSE
    ),

    n_interpeak_diversion_signals = case_when(
      low_replicate_count ~ NA_integer_,
      TRUE ~ as.integer(interpeak_high_flag) +
             as.integer(domain_interpeak_low_flag) +
             as.integer(anchor_domain_high_flag)
    ),

    interpeak_diversion_flag_raw = case_when(
      low_replicate_count ~ NA,
      TRUE ~ n_interpeak_diversion_signals >= INTERPEAK_DIVERSION_SIGNALS_MIN
    ),

    interpeak_diversion_protected = case_when(
      low_replicate_count ~ NA,
      TRUE ~ interpeak_diversion_flag_raw &
             !is.na(exp_median_interpeak) &
             !is.na(exp_domain_range) &
             exp_median_interpeak <= INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX &
             exp_domain_range <= INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX
    ),

    interpeak_diversion_flag = case_when(
      low_replicate_count ~ NA,
      TRUE ~ interpeak_diversion_flag_raw & !interpeak_diversion_protected
    )
  ) %>%
  select(
    sample,
    experiment_group,
    low_replicate_count,
    exp_median_interpeak,
    exp_domain_range,
    interpeak_high_robust_z,
    domain_interpeak_low_robust_z,
    anchor_domain_high_robust_z,
    interpeak_high_flag,
    domain_interpeak_low_flag,
    baseline_domain_interpeak_outlier,
    anchor_domain_high_flag,
    n_interpeak_diversion_signals,
    interpeak_diversion_flag_raw,
    interpeak_diversion_protected,
    interpeak_diversion_flag
  )

write_tsv(
  interpeak_diversion_metrics,
  file.path(sec_dir, "interpeak_diversion_metrics.tsv")
)

genotype_group_diversion_summary <- baseline_metrics %>%
  left_join(
    interpeak_diversion_metrics %>%
      select(sample, experiment_group, interpeak_diversion_flag),
    by = c("sample", "experiment_group")
  ) %>%
  group_by(experiment_group, qc_group, genotype_group, antibody, genotype) %>%
  summarise(
    n_samples_in_genotype_group = n(),
    n_interpeak_diversion_samples = sum(interpeak_diversion_flag, na.rm = TRUE),
    frac_interpeak_diversion_samples = n_interpeak_diversion_samples / n_samples_in_genotype_group,
    genotype_median_interpeak = median(baseline_median_interpeak, na.rm = TRUE),
    genotype_median_domain_interpeak = median(baseline_domain_interpeak, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    baseline_metrics %>%
      group_by(experiment_group) %>%
      summarise(
        experiment_median_interpeak = median(baseline_median_interpeak, na.rm = TRUE),
        experiment_median_domain_interpeak = median(baseline_domain_interpeak, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "experiment_group"
  ) %>%
  mutate(
    classwide_interpeak_high = genotype_median_interpeak > experiment_median_interpeak,
    classwide_domain_interpeak_low = genotype_median_domain_interpeak < experiment_median_domain_interpeak,
    classwide_basin_assumption_violation = (
      n_samples_in_genotype_group >= CLASSWIDE_MIN_SAMPLES &
      n_interpeak_diversion_samples >= CLASSWIDE_MIN_FLAGGED_SAMPLES &
      frac_interpeak_diversion_samples >= CLASSWIDE_MIN_FLAGGED_FRAC &
      classwide_interpeak_high &
      classwide_domain_interpeak_low
    )
  )

write_tsv(
  genotype_group_diversion_summary,
  file.path(sec_dir, "genotype_group_diversion_summary.tsv")
)

experiment_classwide_diversion_summary <- genotype_group_diversion_summary %>%
  group_by(experiment_group) %>%
  summarise(
    n_classwide_diversion_genotype_groups = sum(classwide_basin_assumption_violation, na.rm = TRUE),
    classwide_diversion_genotype_groups = paste(
      genotype_group[classwide_basin_assumption_violation],
      collapse = ";"
    ),
    any_classwide_basin_assumption_violation = any(classwide_basin_assumption_violation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    classwide_diversion_genotype_groups = ifelse(
      classwide_diversion_genotype_groups == "",
      NA_character_,
      classwide_diversion_genotype_groups
    )
  )

genotype_group_inconsistency_summary <- baseline_metrics %>%
  group_by(genotype_group, experiment_group, qc_group, antibody, genotype) %>%
  summarise(
    n_samples_in_genotype_group = n(),
    baseline_domain_interpeak_range = diff(range(baseline_domain_interpeak, na.rm = TRUE)),
    baseline_anchor_interpeak_range = diff(range(baseline_mean_anchor_interpeak, na.rm = TRUE)),
    baseline_domain_interpeak_mad = mad(baseline_domain_interpeak, na.rm = TRUE),
    baseline_anchor_interpeak_mad = mad(baseline_mean_anchor_interpeak, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(genotype_group_inconsistency_summary,
          file.path(sec_dir, "genotype_group_inconsistency_summary.tsv"))

########################################################
# NEW: GENOTYPE SPREAD QC
########################################################
genotype_spread_qc <- baseline_metrics %>%
  group_by(genotype_group, experiment_group, qc_group, antibody, genotype) %>%
  summarise(
    n_samples_in_genotype_group = n(),
    any_low_replicate_count = any(low_replicate_count, na.rm = TRUE),

    baseline_median_domain_median = median(baseline_median_domain, na.rm = TRUE),
    baseline_median_domain_min = min(baseline_median_domain, na.rm = TRUE),
    baseline_median_domain_max = max(baseline_median_domain, na.rm = TRUE),
    baseline_median_domain_range = diff(range(baseline_median_domain, na.rm = TRUE)),

    baseline_median_interpeak_median = median(baseline_median_interpeak, na.rm = TRUE),
    baseline_median_interpeak_min = min(baseline_median_interpeak, na.rm = TRUE),
    baseline_median_interpeak_max = max(baseline_median_interpeak, na.rm = TRUE),
    baseline_median_interpeak_range = diff(range(baseline_median_interpeak, na.rm = TRUE)),

    baseline_domain_interpeak_median = median(baseline_domain_interpeak, na.rm = TRUE),
    baseline_anchor_interpeak_median = median(baseline_mean_anchor_interpeak, na.rm = TRUE),

    baseline_domain_interpeak_min = min(baseline_domain_interpeak, na.rm = TRUE),
    baseline_domain_interpeak_max = max(baseline_domain_interpeak, na.rm = TRUE),
    baseline_anchor_interpeak_min = min(baseline_mean_anchor_interpeak, na.rm = TRUE),
    baseline_anchor_interpeak_max = max(baseline_mean_anchor_interpeak, na.rm = TRUE),

    baseline_domain_interpeak_range = diff(range(baseline_domain_interpeak, na.rm = TRUE)),
    baseline_anchor_interpeak_range = diff(range(baseline_mean_anchor_interpeak, na.rm = TRUE)),

    baseline_domain_interpeak_ratio_spread = ifelse(
      any_low_replicate_count | isTRUE(all.equal(min(baseline_domain_interpeak, na.rm = TRUE), 0)),
      NA_real_,
      max(baseline_domain_interpeak, na.rm = TRUE) / min(baseline_domain_interpeak, na.rm = TRUE)
    ),
    baseline_anchor_interpeak_ratio_spread = ifelse(
      any_low_replicate_count | isTRUE(all.equal(min(baseline_mean_anchor_interpeak, na.rm = TRUE), 0)),
      NA_real_,
      max(baseline_mean_anchor_interpeak, na.rm = TRUE) / min(baseline_mean_anchor_interpeak, na.rm = TRUE)
    ),

    baseline_domain_interpeak_cv = ifelse(
      any_low_replicate_count | isTRUE(all.equal(mean(baseline_domain_interpeak, na.rm = TRUE), 0)),
      NA_real_,
      sd(baseline_domain_interpeak, na.rm = TRUE) / abs(mean(baseline_domain_interpeak, na.rm = TRUE))
    ),
    baseline_anchor_interpeak_cv = ifelse(
      any_low_replicate_count | isTRUE(all.equal(mean(baseline_mean_anchor_interpeak, na.rm = TRUE), 0)),
      NA_real_,
      sd(baseline_mean_anchor_interpeak, na.rm = TRUE) / abs(mean(baseline_mean_anchor_interpeak, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>%
  mutate(
    domain_ratio_high_spread_flag = case_when(
      any_low_replicate_count ~ NA,
      n_samples_in_genotype_group >= 2 &
        !is.na(baseline_domain_interpeak_ratio_spread) &
        baseline_domain_interpeak_ratio_spread >= GENOTYPE_SPREAD_RATIO_WARN ~ TRUE,
      TRUE ~ FALSE
    ),

    anchor_ratio_high_spread_flag = case_when(
      any_low_replicate_count ~ NA,
      n_samples_in_genotype_group >= 2 &
        !is.na(baseline_anchor_interpeak_ratio_spread) &
        baseline_anchor_interpeak_ratio_spread >= GENOTYPE_SPREAD_RATIO_WARN ~ TRUE,
      TRUE ~ FALSE
    ),

    high_spread_flag = case_when(
      any_low_replicate_count ~ NA,
      TRUE ~ domain_ratio_high_spread_flag | anchor_ratio_high_spread_flag
    ),

    high_cv_flag = case_when(
      any_low_replicate_count ~ NA,
      n_samples_in_genotype_group >= 2 & (
        (!is.na(baseline_domain_interpeak_cv) &
           baseline_domain_interpeak_cv >= GENOTYPE_SPREAD_CV_WARN) |
        (!is.na(baseline_anchor_interpeak_cv) &
           baseline_anchor_interpeak_cv >= GENOTYPE_SPREAD_CV_WARN)
      ) ~ TRUE,
      TRUE ~ FALSE
    ),

    low_interpeak_protected = case_when(
      any_low_replicate_count ~ NA,
      high_spread_flag &
        !is.na(baseline_median_interpeak_median) &
        !is.na(baseline_median_domain_range) &
        baseline_median_interpeak_median <= GENOTYPE_LOW_INTERPEAK_PROTECT_MAX &
        baseline_median_domain_range <= GENOTYPE_DOMAIN_RANGE_PROTECT_MAX ~ TRUE,
      TRUE ~ FALSE
    ),

    genotype_inconsistent_flag = case_when(
      any_low_replicate_count ~ NA,
      TRUE ~ high_spread_flag & !low_interpeak_protected
    )
  )

write_tsv(genotype_spread_qc,
          file.path(sec_dir, "genotype_spread_qc.tsv"))

experiment_genotype_spread_summary <- genotype_spread_qc %>%
  group_by(experiment_group, qc_group, antibody) %>%
  summarise(
    any_genotype_high_spread = any(high_spread_flag, na.rm = TRUE),
    any_genotype_inconsistent = any(genotype_inconsistent_flag, na.rm = TRUE),
    n_low_interpeak_protected_groups = sum(low_interpeak_protected, na.rm = TRUE),
    low_interpeak_protected_groups = paste(genotype_group[low_interpeak_protected], collapse = ";"),
    n_genotype_groups_high_spread = sum(high_spread_flag, na.rm = TRUE),
    n_genotype_groups_inconsistent = sum(genotype_inconsistent_flag, na.rm = TRUE),
    max_ratio_spread_domain_interpeak = suppressWarnings(max(baseline_domain_interpeak_ratio_spread, na.rm = TRUE)),
    max_ratio_spread_anchor_interpeak = suppressWarnings(max(baseline_anchor_interpeak_ratio_spread, na.rm = TRUE)),
    inconsistent_genotype_groups = paste(genotype_group[genotype_inconsistent_flag], collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    max_ratio_spread_domain_interpeak = ifelse(is.infinite(max_ratio_spread_domain_interpeak), NA_real_, max_ratio_spread_domain_interpeak),
    max_ratio_spread_anchor_interpeak = ifelse(is.infinite(max_ratio_spread_anchor_interpeak), NA_real_, max_ratio_spread_anchor_interpeak),
    inconsistent_genotype_groups = ifelse(inconsistent_genotype_groups == "", NA_character_, inconsistent_genotype_groups),
    low_interpeak_protected_groups = ifelse(low_interpeak_protected_groups == "", NA_character_, low_interpeak_protected_groups),
  )

sample_threshold_summary <- scaled_metrics %>%
  group_by(sample, qc_group, experiment_group, genotype_group, antibody, genotype, replicate) %>%
  summarise(
    n_thresholds_present = sum(!is.na(threshold)),
    thresholds_present = paste(sort(unique(threshold)), collapse = ","),
    domain_interpeak_range = diff(range(domain_interpeak, na.rm = TRUE)),
    anchor_interpeak_range = diff(range(mean_anchor_interpeak, na.rm = TRUE)),
    .groups = "drop"
  )

# Trimmed anchor medians at the primary drop. Current HARMONY writes the exact
# retained BAR identities used for scaling. Legacy runs fall back to a corrected
# reconstruction using global all-zero filtering and round(N * drop_fraction).
primary_anchor_all <- all_long %>%
  filter(
    region_class == "anchor",
    set_label != "baseline",
    !is.na(signal),
    abs(suppressWarnings(as.numeric(str_extract(set_label, "[0-9.]+"))) - primary_drop) < 1e-8
  )

if (has_exact_retained_bars) {
  primary_retained_ids <- retained_bars_exact %>%
    filter(abs(threshold - primary_drop) < 1e-8) %>%
    select(sample, region_id) %>%
    distinct()
  primary_anchor_trimmed_by_sample <- primary_anchor_all %>%
    inner_join(primary_retained_ids, by = c("sample", "region_id")) %>%
    group_by(sample) %>%
    summarise(primary_anchor_trimmed = median(signal, na.rm = TRUE), .groups = "drop")
} else {
  primary_anchor_trimmed_by_sample <- primary_anchor_all %>%
    group_by(sample) %>%
    arrange(desc(signal), .by_group = TRUE) %>%
    mutate(
      n_anchor_rows = n(),
      rank_desc_signal = row_number(),
      n_trim = round(n_anchor_rows * primary_drop)
    ) %>%
    filter(rank_desc_signal > n_trim) %>%
    summarise(primary_anchor_trimmed = median(signal, na.rm = TRUE), .groups = "drop")
}

primary_scaled <- scaled_metrics %>%
  filter(abs(threshold - primary_drop) < 1e-8) %>%
  left_join(primary_anchor_trimmed_by_sample, by = "sample") %>%
  group_by(experiment_group) %>%
  mutate(
    group_n = n(),
    primary_domain_z = ifelse(group_n >= 3, as.numeric(scale(median_domain)), NA_real_),
    primary_interpeak_z = ifelse(group_n >= 3, as.numeric(scale(median_interpeak)), NA_real_)
  ) %>%
  ungroup() %>%
  group_by(experiment_group) %>%
  mutate(
    chrom_group_n = n(),
    primary_sex_chrom_retained_anchor_fraction_z = ifelse(
      chrom_group_n >= 3,
      as.numeric(scale(sex_chrom_retained_anchor_fraction)),
      NA_real_
    ),
    # Post-BASIN anchor convergence: after correct scaling, every sample's
    # TRIMMED scaled anchor median (the kept set at the primary drop, with the
    # top primary_drop fraction of anchors removed -- the same anchors the scale
    # factor was built from) should converge within the experiment group at ~1.0.
    # We measure each sample's deviation with a robust (MAD-based) z and the
    # group's overall convergence with a coefficient of variation.
    primary_anchor_group_median = median(primary_anchor_trimmed, na.rm = TRUE),

    primary_anchor_group_cv = ifelse(
      mean(primary_anchor_trimmed, na.rm = TRUE) != 0,
      sd(primary_anchor_trimmed, na.rm = TRUE) / mean(primary_anchor_trimmed, na.rm = TRUE),
      NA_real_
    ),

  ) %>%
  ungroup() %>%
  select(
    sample,
    qc_group,
    primary_threshold = threshold,
    primary_scale_factor = scale_factor,
    primary_domain = median_domain,
    primary_anchor = median_anchor,
    primary_anchor_trimmed,
    primary_interpeak = median_interpeak,
    primary_mean_anchor = mean_anchor,
    primary_mean_interpeak = mean_interpeak,
    primary_domain_interpeak = domain_interpeak,
    primary_anchor_interpeak = mean_anchor_interpeak,
    primary_anchor_domain = anchor_domain,
    primary_domain_z,
    primary_interpeak_z,
    primary_anchor_group_median,
    primary_anchor_group_cv,
    total_retained_anchor_rows,
    frac_retained_anchor_rows_autosomes,
    frac_retained_anchor_rows_X,
    frac_retained_anchor_rows_Y,
    sex_chrom_retained_anchor_fraction,
    primary_sex_chrom_retained_anchor_fraction_z
  ) %>%
  mutate(
    post_basin_group_outlier = ifelse(
      is.na(primary_domain_z) | is.na(primary_interpeak_z),
      FALSE,
      abs(primary_domain_z) > PRIMARY_OUTLIER_Z | abs(primary_interpeak_z) > PRIMARY_OUTLIER_Z
    ),
    sex_chrom_anchor_skew = ifelse(
      is.na(primary_sex_chrom_retained_anchor_fraction_z),
      FALSE,
      abs(primary_sex_chrom_retained_anchor_fraction_z) > SEX_CHROM_SKEW_Z
    ),
    post_basin_anchor_nonconvergence = ifelse(
      is.na(primary_anchor_group_cv),
      FALSE,
      primary_anchor_group_cv > POST_BASIN_ANCHOR_CV_WARN
    )
  )

# Post-BASIN pipeline-integrity table (secondary output). One row per sample at
# the primary drop threshold. Because the canonical scale factor is defined as
# 1 / retained-BAR median and applied to the same canonical representation, the
# retained scaled BAR median should be ~1.0 by construction. This is therefore a
# mechanical reconciliation check (wrong BAR file, scale factor not applied,
# file mismatch), NOT an independent biological validation metric.
post_basin_anchor_check <- primary_scaled %>%
  select(
    sample, qc_group,
    primary_threshold,
    primary_anchor_trimmed,
    primary_anchor_group_median,
    primary_anchor_group_cv,
    post_basin_anchor_nonconvergence
  ) %>%
  arrange(qc_group, sample)
write_tsv(post_basin_anchor_check, file.path(sec_dir, "post_basin_anchor_check.tsv"))

########################################################
# DROPOUT QC SUMMARIES
########################################################

thresholds_run <- sort(unique(na.omit(scales$threshold)))

sample_threshold_grid <- manifest %>%
  select(sample, qc_group, experiment_group, genotype_group, antibody, genotype, replicate) %>%
  distinct() %>%
  crossing(threshold = thresholds_run)

sample_dropout_summary <- sample_threshold_grid %>%
  left_join(scales, by = c("sample", "threshold")) %>%
  mutate(
    survived_threshold = !is.na(scale_factor),
    threshold_label = paste0("drop", formatC(threshold, format = "f", digits = 2))
  ) %>%
  group_by(sample, qc_group, experiment_group, genotype_group, antibody, genotype, replicate) %>%
  summarise(
    n_thresholds_run = n_distinct(threshold),
    thresholds_run = paste(sort(unique(threshold)), collapse = ","),
    n_thresholds_survived = n_distinct(threshold[survived_threshold]),
    n_thresholds_dropped = n_distinct(threshold[!survived_threshold]),
    dropped_any_threshold = any(!survived_threshold),
    dropped_drop10 = any(abs(threshold - 0.10) < 1e-8 & !survived_threshold),
    dropped_drop50 = any(abs(threshold - 0.50) < 1e-8 & !survived_threshold),
    thresholds_dropped = paste(sort(unique(threshold[!survived_threshold])), collapse = ","),
    dropped_threshold_labels = paste(sort(unique(threshold_label[!survived_threshold])), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    thresholds_dropped = ifelse(thresholds_dropped == "", NA_character_, thresholds_dropped),
    dropped_threshold_labels = ifelse(dropped_threshold_labels == "", NA_character_, dropped_threshold_labels)
  )

write_tsv(sample_dropout_summary,
          file.path(sec_dir, "sample_dropout_summary.tsv"))

experiment_threshold_dropout_summary <- sample_threshold_grid %>%
  left_join(scales, by = c("sample", "threshold")) %>%
  mutate(
    survived_threshold = !is.na(scale_factor)
  ) %>%
  group_by(experiment_group, qc_group, antibody, threshold) %>%
  summarise(
    n_samples_expected = n_distinct(sample),
    n_samples_survived = n_distinct(sample[survived_threshold]),
    n_samples_dropped = n_distinct(sample[!survived_threshold]),
    frac_samples_dropped = n_samples_dropped / n_samples_expected,
    dropped_samples = paste(sort(unique(sample[!survived_threshold])), collapse = ","),
    any_dropout_this_threshold = n_samples_dropped > 0,
    hard_fail_this_threshold = case_when(
      abs(threshold - 0.10) < 1e-8 ~ frac_samples_dropped > DROP10_FAIL_FRAC,
      abs(threshold - 0.50) < 1e-8 ~ frac_samples_dropped > DROP50_FAIL_FRAC,
      TRUE ~ FALSE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    dropped_samples = ifelse(dropped_samples == "", NA_character_, dropped_samples)
  ) %>%
  distinct() %>%
  arrange(experiment_group, threshold)

write_tsv(experiment_threshold_dropout_summary,
          file.path(sec_dir, "experiment_threshold_dropout_summary.tsv"))

experiment_dropout_summary <- experiment_threshold_dropout_summary %>%
  group_by(experiment_group, qc_group, antibody) %>%
  summarise(
    any_dropout_any_threshold = any(any_dropout_this_threshold, na.rm = TRUE),
    thresholds_with_any_dropout = paste(sort(unique(threshold[any_dropout_this_threshold])), collapse = ","),
    dropout_threshold_labels = paste(
      paste0(
        "drop",
        formatC(threshold[any_dropout_this_threshold], format = "f", digits = 2),
        "_n",
        n_samples_dropped[any_dropout_this_threshold],
        "of",
        n_samples_expected[any_dropout_this_threshold]
      ),
      collapse = ";"
    ),
    drop10_frac_dropped = ifelse(any(abs(threshold - 0.10) < 1e-8),
                                 frac_samples_dropped[which.min(abs(threshold - 0.10))],
                                 NA_real_),
    drop50_frac_dropped = ifelse(any(abs(threshold - 0.50) < 1e-8),
                                 frac_samples_dropped[which.min(abs(threshold - 0.50))],
                                 NA_real_),
    hard_fail_drop10 = any(abs(threshold - 0.10) < 1e-8 & frac_samples_dropped > DROP10_FAIL_FRAC, na.rm = TRUE),
    hard_fail_drop50 = any(abs(threshold - 0.50) < 1e-8 & frac_samples_dropped > DROP50_FAIL_FRAC, na.rm = TRUE),
    experiment_dropout_hard_fail = hard_fail_drop10 | hard_fail_drop50,
    experiment_manual_inspection_required = any_dropout_any_threshold,
    .groups = "drop"
  ) %>%
  mutate(
    thresholds_with_any_dropout = ifelse(thresholds_with_any_dropout == "", NA_character_, thresholds_with_any_dropout),
    dropout_threshold_labels = ifelse(dropout_threshold_labels == "", NA_character_, dropout_threshold_labels)
  )

########################################################
# SAMPLE QC SUMMARY
########################################################

sample_qc_summary <- manifest %>%
  mutate(replicate = as.character(replicate)) %>%
  left_join(
    baseline_metrics %>%
      mutate(replicate = as.character(replicate)) %>%
      select(-replicate, -n_reps_in_genotype_group, -low_replicate_count,
             -very_low_replicate_count, -low_replicate_note,
             -experiment_low_replicate_note),
    by = c("sample", "qc_group", "experiment_group", "genotype_group", "antibody", "genotype")
  ) %>%
  left_join(
    interpeak_diversion_metrics %>% select(-low_replicate_count),
    by = c("sample", "experiment_group")
  ) %>%
  left_join(
    sample_threshold_summary %>% mutate(replicate = as.character(replicate)) %>% select(-replicate),
    by = c("sample", "qc_group", "experiment_group", "genotype_group", "antibody", "genotype")
  ) %>%
  left_join(
    sample_dropout_summary %>% mutate(replicate = as.character(replicate)) %>% select(-replicate),
    by = c("sample", "qc_group", "experiment_group", "genotype_group", "antibody", "genotype")
  ) %>%
  left_join(primary_scaled, by = c("sample", "qc_group")) %>%
  mutate(
    # Derive low_replicate_count unconditionally from n_reps_in_genotype_group
    # which is always present from the manifest. This handles the case where
    # the baseline_metrics join dropped the column due to type mismatch.
    low_replicate_count = !is.na(n_reps_in_genotype_group) & n_reps_in_genotype_group < 3,
    baseline_weird = ifelse(is.na(baseline_weird), FALSE, baseline_weird),
    low_anchor_support = ifelse(is.na(low_anchor_support), FALSE, low_anchor_support),
    interpeak_diversion_flag_for_summary = ifelse(
      low_replicate_count, FALSE,
      ifelse(is.na(interpeak_diversion_flag), FALSE, interpeak_diversion_flag)
    ),
    baseline_domain_interpeak_outlier = ifelse(
      low_replicate_count, FALSE,
      ifelse(is.na(baseline_domain_interpeak_outlier), FALSE, baseline_domain_interpeak_outlier)
    ),
    post_basin_group_outlier = ifelse(is.na(post_basin_group_outlier), FALSE, post_basin_group_outlier),
    sex_chrom_anchor_skew = ifelse(is.na(sex_chrom_anchor_skew), FALSE, sex_chrom_anchor_skew),
    post_basin_anchor_nonconvergence = ifelse(is.na(post_basin_anchor_nonconvergence), FALSE, post_basin_anchor_nonconvergence),
    dropped_any_threshold = ifelse(is.na(dropped_any_threshold), FALSE, dropped_any_threshold),
    dropped_drop10 = ifelse(is.na(dropped_drop10), FALSE, dropped_drop10),
    dropped_drop50 = ifelse(is.na(dropped_drop50), FALSE, dropped_drop50)
  ) %>%
  rowwise() %>%
  mutate(
    flag_components = paste(
      c(
        if (!is.na(low_replicate_note)) paste0("low_replicate_count(", n_reps_in_genotype_group, ")") else NULL,
        if (is.na(primary_scale_factor)) "missing_primary_drop" else NULL,
        if (dropped_any_threshold) "dropped_at_some_threshold" else NULL,
        if (baseline_weird) "baseline_weird" else NULL,
        if (!low_replicate_count && baseline_domain_interpeak_outlier) "baseline_domain_interpeak_outlier" else NULL,
        if (low_anchor_support) "low_anchor_support" else NULL,
        if (!low_replicate_count && interpeak_diversion_flag_for_summary) "interpeak_diversion" else NULL,
        if (post_basin_group_outlier) "post_basin_group_outlier" else NULL,
        if (sex_chrom_anchor_skew) "sex_chrom_anchor_skew" else NULL
      ),
      collapse = "|"
    )
  ) %>%
  ungroup() %>%
  mutate(
    flag_components = ifelse(flag_components == "", "ok", flag_components),
    sample_flag = case_when(
      is.na(primary_scale_factor) ~ "missing_primary_drop",
      dropped_any_threshold & low_anchor_support ~ "suspect",
      baseline_weird & low_anchor_support ~ "suspect",
      interpeak_diversion_flag_for_summary & low_anchor_support ~ "suspect",
      low_replicate_count | dropped_any_threshold | baseline_weird | low_anchor_support |
        interpeak_diversion_flag_for_summary | baseline_domain_interpeak_outlier |
        post_basin_group_outlier | sex_chrom_anchor_skew ~ "warning",
      TRUE ~ "ok"
    ),
    sample_note = case_when(
      !is.na(low_replicate_note) & sample_flag == "warning" ~
        paste0("Warning: ", low_replicate_note, "."),
      sample_flag == "missing_primary_drop" ~
        "Primary drop missing for this sample.",
      sample_flag == "suspect" & dropped_any_threshold & low_anchor_support ~
        "Sample drops out at one or more BASIN thresholds and also has poor baseline anchor support.",
      sample_flag == "warning" & baseline_domain_interpeak_outlier ~
        "Baseline domain/interpeak ratio is unusually high or low relative to the experiment group; inspect normalization behavior.",
      sample_flag == "suspect" & baseline_weird & low_anchor_support ~
        "Baseline weirdness coincides with poor anchor support.",
      sample_flag == "suspect" & interpeak_diversion_flag_for_summary & low_anchor_support ~
        "Sample shows an interpeak-diversion pattern together with poor anchor support.",
      sample_flag == "warning" & dropped_any_threshold ~
        "Sample drops out at one or more BASIN thresholds; inspect anchor support and scale computation.",
      sample_flag == "warning" & interpeak_diversion_flag_for_summary ~
        "Sample shows an interpeak-diversion pattern: interpeak high, domain/interpeak low, and/or anchor/domain high.",
      sample_flag == "warning" & baseline_weird ~
        "Baseline region-class behavior is unusual for this experiment group.",
      sample_flag == "warning" & low_anchor_support ~
        "Anchor support is low: few nonzero informative anchor rows.",
      sample_flag == "warning" & post_basin_group_outlier ~
        "Sample is >2 SD from the experiment-group mean at the primary BASIN drop.",
      sample_flag == "warning" & sex_chrom_anchor_skew ~
        "Retained anchor composition shows sex-chromosome skew; consider rerunning with sex chromosomes removed from the anchor file.",
      sample_flag == "ok" ~
        "No major issue detected by current QC rules.",
      TRUE ~ "QC note unavailable."
    )
  ) %>%
  arrange(experiment_group, sample)

write_tsv(sample_qc_summary, file.path(outdir, "sample_qc_summary.tsv"))

########################################################
# EXPERIMENT QC SUMMARY
########################################################
# Decision policy (recommended next action):
#   pass: proceed
#   pass with caveat: proceed while noting the stated limitation
#   manual review: examine the identified issue before proceeding
#   failed QC: do not proceed without remediation
#
# The underlying decision thresholds and ordering are unchanged.

experiment_qc_summary <- sample_qc_summary %>%
  left_join(
    genotype_group_inconsistency_summary %>%
      group_by(experiment_group) %>%
      summarise(
        max_within_genotype_domain_interpeak_range = max(baseline_domain_interpeak_range, na.rm = TRUE),
        max_within_genotype_anchor_interpeak_range = max(baseline_anchor_interpeak_range, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "experiment_group"
  ) %>%
  left_join(
    experiment_dropout_summary,
    by = c("experiment_group", "qc_group", "antibody")
  ) %>%
  left_join(
    experiment_classwide_diversion_summary,
    by = "experiment_group"
  ) %>%
  left_join(
    experiment_genotype_spread_summary,
    by = c("experiment_group", "qc_group", "antibody")
  ) %>%
  group_by(experiment_group, qc_group, antibody) %>%
  summarise(
    antibody = first(antibody),
    qc_group = first(qc_group),
    fractions_present = paste(sort(unique(na.omit(fraction))), collapse = ";"),
    n_samples = n(),
    any_low_replicate_count = any(low_replicate_count, na.rm = TRUE),
    min_reps_in_any_genotype_group = min(n_reps_in_genotype_group, na.rm = TRUE),
    experiment_low_replicate_note = first(na.omit(experiment_low_replicate_note)),
    recommended_drop = ifelse(any(abs(primary_threshold - primary_drop) < 1e-8, na.rm = TRUE), primary_drop, NA_real_),
    n_warning_samples = sum(sample_flag == "warning" & !low_replicate_count, na.rm = TRUE),
    n_low_rep_warning_samples = sum(sample_flag == "warning" & low_replicate_count, na.rm = TRUE),
    n_suspect_samples = sum(sample_flag == "suspect", na.rm = TRUE),
    n_interpeak_diversion_samples = sum(interpeak_diversion_flag_for_summary, na.rm = TRUE),
    frac_interpeak_diversion_samples = n_interpeak_diversion_samples / n_samples,
    any_post_basin_anchor_nonconvergence = any(post_basin_anchor_nonconvergence, na.rm = TRUE),
    max_post_basin_anchor_group_cv = suppressWarnings(max(primary_anchor_group_cv, na.rm = TRUE)),
    fraction_flagged = (n_warning_samples + n_suspect_samples) / n_samples,
    # Count how broadly non-low-replicate warning/suspect flags are distributed.
    # This distinguishes experiment-wide instability from a coherent pattern
    # confined to one genotype group.
    n_flagged_genotype_groups = n_distinct(
      genotype_group[
        sample_flag %in% c("warning", "suspect") &
        !low_replicate_count
      ]
    ),
    flagged_genotype_groups = paste(
      sort(unique(genotype_group[
        sample_flag %in% c("warning", "suspect") &
        !low_replicate_count
      ])),
      collapse = ";"
    ),
    threshold_coverage = paste(sort(unique(na.omit(thresholds_present))), collapse = ";"),
    max_within_genotype_domain_interpeak_range = first(max_within_genotype_domain_interpeak_range),
    max_within_genotype_anchor_interpeak_range = first(max_within_genotype_anchor_interpeak_range),
    any_genotype_high_spread = first(any_genotype_high_spread),
    any_genotype_inconsistent = first(any_genotype_inconsistent),
    n_genotype_groups_high_spread = first(n_genotype_groups_high_spread),
    n_genotype_groups_inconsistent = first(n_genotype_groups_inconsistent),
    max_ratio_spread_domain_interpeak = first(max_ratio_spread_domain_interpeak),
    max_ratio_spread_anchor_interpeak = first(max_ratio_spread_anchor_interpeak),
    inconsistent_genotype_groups = first(inconsistent_genotype_groups),
    any_dropout_any_threshold = first(any_dropout_any_threshold),
    thresholds_with_any_dropout = first(thresholds_with_any_dropout),
    dropout_threshold_labels = first(dropout_threshold_labels),
    drop10_frac_dropped = first(drop10_frac_dropped),
    drop50_frac_dropped = first(drop50_frac_dropped),
    hard_fail_drop10 = first(hard_fail_drop10),
    hard_fail_drop50 = first(hard_fail_drop50),
    experiment_dropout_hard_fail = first(experiment_dropout_hard_fail),
    experiment_manual_inspection_required = first(experiment_manual_inspection_required),
    any_classwide_basin_assumption_violation = first(any_classwide_basin_assumption_violation),
    n_classwide_diversion_genotype_groups = first(n_classwide_diversion_genotype_groups),
    classwide_diversion_genotype_groups = first(classwide_diversion_genotype_groups),
    .groups = "drop"
  ) %>%
  mutate(
    max_post_basin_anchor_group_cv = ifelse(is.infinite(max_post_basin_anchor_group_cv), NA_real_, max_post_basin_anchor_group_cv)
  ) %>%
  mutate(
    experiment_recommendation = case_when(
      experiment_dropout_hard_fail ~ "failed QC",
      n_genotype_groups_inconsistent >= FAIL_INCONSISTENT_GENOTYPE_GROUPS ~ "failed QC",
      !is.na(max_ratio_spread_domain_interpeak) &
        max_ratio_spread_domain_interpeak >= FAIL_MAX_DOMAIN_RATIO_SPREAD &
        any_genotype_inconsistent ~ "failed QC",
      fraction_flagged >= FAIL_FRACTION_FLAGGED &
        n_flagged_genotype_groups >= 2 ~ "failed QC",

      # A high flagged fraction confined to one genotype is not evidence of
      # broad experiment-wide instability. Escalate it for genotype-specific
      # review rather than automatically hard-failing the whole experiment.
      fraction_flagged >= FAIL_FRACTION_FLAGGED &
        n_flagged_genotype_groups == 1 ~ "manual review",
      n_suspect_samples > 0 ~ "manual review",
      n_warning_samples >= MANUAL_WARNING_SAMPLE_MIN ~ "manual review",
      n_genotype_groups_inconsistent == 1 ~ "manual review",
      any_classwide_basin_assumption_violation ~ "manual review",
      n_interpeak_diversion_samples >= 2 ~ "manual review",
      any_post_basin_anchor_nonconvergence ~ "manual review",
      experiment_manual_inspection_required ~ "manual review",

      # Low replicate count (n=1 or n=2) is a warning only — it does not escalate
      # to manual inspection required. The per-sample low_replicate_note and the
      # experiment_low_replicate_note carry the relevant detail.
      any_low_replicate_count ~ "pass with caveat",

      n_warning_samples >= 1 & n_warning_samples <= WARNING_SAMPLE_MAX &
        n_suspect_samples == 0 &
        !any_genotype_inconsistent &
        (
          !any_dropout_any_threshold |
          (
            !is.na(dropout_threshold_labels) &
            is.na(drop10_frac_dropped) == FALSE &
            drop10_frac_dropped == 0 &
            !is.na(drop50_frac_dropped) &
            drop50_frac_dropped <= (1 / n_samples)
          )
        ) ~ "pass with caveat",

      TRUE ~ "pass"
    ),

    decision_trigger = case_when(
      experiment_dropout_hard_fail & hard_fail_drop10 & hard_fail_drop50 ~
        "dropout failed-QC threshold reached at drop0.10 and drop0.50",
      experiment_dropout_hard_fail & hard_fail_drop10 ~
        "dropout failed-QC threshold reached at drop0.10",
      experiment_dropout_hard_fail & hard_fail_drop50 ~
        "dropout failed-QC threshold reached at drop0.50",
      n_genotype_groups_inconsistent >= FAIL_INCONSISTENT_GENOTYPE_GROUPS ~
        ">=2 inconsistent genotype groups",
      !is.na(max_ratio_spread_domain_interpeak) &
        max_ratio_spread_domain_interpeak >= FAIL_MAX_DOMAIN_RATIO_SPREAD &
        any_genotype_inconsistent ~
        "max domain/interpeak ratio spread >= 2.25 with genotype inconsistency",
      fraction_flagged >= FAIL_FRACTION_FLAGGED &
        n_flagged_genotype_groups >= 2 ~
        ">=1/3 of samples flagged across >=2 genotype groups",
      fraction_flagged >= FAIL_FRACTION_FLAGGED &
        n_flagged_genotype_groups == 1 ~
        ">=1/3 of samples flagged, concentrated in one genotype group",

      n_suspect_samples > 0 ~
        ">=1 suspect sample",
      n_warning_samples >= MANUAL_WARNING_SAMPLE_MIN ~
        ">=3 warning samples",
      n_genotype_groups_inconsistent == 1 ~
        "exactly 1 inconsistent genotype group",
      any_classwide_basin_assumption_violation ~
        "class-wide BASIN assumption violation",
      n_interpeak_diversion_samples >= 2 ~
        ">=2 interpeak-diversion samples",
      any_post_basin_anchor_nonconvergence ~
        "pipeline integrity: scaled anchor signal failed to converge (post-BASIN sanity check -- likely a mechanical issue, not biology)",
      experiment_manual_inspection_required ~
        "dropout at one or more thresholds",
      any_low_replicate_count & min_reps_in_any_genotype_group == 1 ~
        "n=1 in one or more genotype groups — z-score metrics undefined",
      any_low_replicate_count & min_reps_in_any_genotype_group == 2 ~
        "n=2 in one or more genotype groups — z-score metrics suppressed",

      n_warning_samples >= 1 & n_warning_samples <= WARNING_SAMPLE_MAX ~
        "1-2 warning samples",

      TRUE ~
        "no major experiment-level QC trigger"
    ),

experiment_note = case_when(
  experiment_dropout_hard_fail &
    hard_fail_drop10 &
    hard_fail_drop50 ~
    paste0(
      "Failed QC: too many samples dropped out at both drop0.10 and drop0.50. ",
      "Do not proceed without remediation."
    ),

  experiment_dropout_hard_fail &
    hard_fail_drop10 ~
    paste0(
      "Failed QC: more than 25% of samples dropped out at drop0.10. ",
      "Do not proceed without remediation."
    ),

  experiment_dropout_hard_fail &
    hard_fail_drop50 ~
    paste0(
      "Failed QC: more than 50% of samples dropped out at drop0.50. ",
      "Do not proceed without remediation."
    ),

  n_genotype_groups_inconsistent >=
    FAIL_INCONSISTENT_GENOTYPE_GROUPS ~
    paste0(
      "Failed QC: multiple genotype groups show within-genotype instability, ",
      "suggesting an unreliable experiment. Do not proceed without remediation."
    ),

  !is.na(max_ratio_spread_domain_interpeak) &
    max_ratio_spread_domain_interpeak >=
      FAIL_MAX_DOMAIN_RATIO_SPREAD &
    any_genotype_inconsistent ~
    paste0(
      "Failed QC: domain/interpeak spread is extreme within at least one ",
      "genotype group and replicates appear inconsistent. ",
      "Do not proceed without remediation."
    ),

  fraction_flagged >= FAIL_FRACTION_FLAGGED &
    n_flagged_genotype_groups >= 2 ~
    paste0(
      "Failed QC: at least one-third of samples are flagged across multiple ",
      "genotype groups, indicating broad experiment-level instability. ",
      "Do not proceed without remediation."
    ),

  fraction_flagged >= FAIL_FRACTION_FLAGGED &
    n_flagged_genotype_groups == 1 ~
    paste0(
      "Manual review: at least one-third of samples are flagged, ",
      "but the pattern is concentrated in one genotype group (",
      flagged_genotype_groups,
      ") rather than distributed across the experiment."
    ),

  n_suspect_samples > 0 ~
    paste0(
      "Manual review: one or more samples are suspect. ",
      "Examine the identified issue before proceeding."
    ),

  n_warning_samples >= MANUAL_WARNING_SAMPLE_MIN ~
    paste0(
      "Manual review: multiple warning-level samples are present. ",
      "Examine the identified issues before proceeding."
    ),

  n_genotype_groups_inconsistent == 1 ~
    paste0(
      "Manual review: one genotype group shows within-genotype instability. ",
      "Examine this group before proceeding."
    ),

  any_classwide_basin_assumption_violation ~
    paste0(
      "Manual review: possible class-wide BASIN assumption violation; ",
      "between-group comparisons may be unreliable. ",
      "Examine this issue before proceeding."
    ),

  n_interpeak_diversion_samples >= 2 ~
    paste0(
      "Manual review: multiple samples show an interpeak-diversion pattern. ",
      "Examine the identified samples before proceeding."
    ),

  any_post_basin_anchor_nonconvergence ~
    paste0(
      "Manual review: post-BASIN anchor convergence exceeded the configured ",
      "CV threshold, indicating that scaled anchor medians were less consistent ",
      "across samples than expected. This may reflect a normalization, ",
      "file-matching, or configuration problem. Inspect ",
      "post_basin_anchor_check.tsv and confirm that HARMONY and CODA used the ",
      "same anchor BED and matching scaled bigWigs before proceeding."
    ),

  experiment_manual_inspection_required ~
    paste0(
      "Manual review: at least one sample dropped out at a run threshold. ",
      "Examine the dropout before proceeding."
    ),

  any_low_replicate_count &
    min_reps_in_any_genotype_group == 1 ~
    paste0(
      "Pass with caveat: one or more genotype groups have only 1 sample. ",
      "Z-score and consistency metrics are undefined and have been suppressed. ",
      "Review the sample QC summary; CODA's statistical framework requires ",
      "at least 3 replicates per group."
    ),

  any_low_replicate_count &
    min_reps_in_any_genotype_group == 2 ~
    paste0(
      "Pass with caveat: one or more genotype groups have only 2 samples. ",
      "Z-score metrics are always +/-0.67 with n=2 regardless of true signal, ",
      "and have been suppressed as uninformative. ",
      "Interpret remaining QC metrics with caution."
    ),

  n_warning_samples >= 1 &
    n_warning_samples <= WARNING_SAMPLE_MAX ~
    paste0(
      "Pass with caveat: one or two samples are flagged, but the experiment ",
      "is otherwise coherent. Proceed while noting the stated limitation."
    ),

  TRUE ~
    "Pass: no major experiment-level QC trigger was detected. Proceed."
  )
)

write_tsv(
  experiment_qc_summary,
  file.path(outdir, "experiment_qc_summary.tsv")
)

sample_qc_by_threshold <- scaled_metrics %>%
  select(
    sample, qc_group, experiment_group, genotype_group, antibody, genotype, replicate, fraction,
    threshold, scale_factor,
    median_domain, median_anchor, median_interpeak,
    mean_anchor, mean_interpeak,
    n_nonzero_anchor_rows, frac_nonzero_anchor_rows,
    total_retained_anchor_rows,
    frac_retained_anchor_rows_autosomes,
    frac_retained_anchor_rows_X,
    frac_retained_anchor_rows_Y,
    sex_chrom_retained_anchor_fraction,
    domain_interpeak, mean_anchor_interpeak, anchor_domain
  ) %>%
  arrange(experiment_group, sample, threshold)

write_tsv(sample_qc_by_threshold, file.path(sec_dir, "sample_qc_by_threshold.tsv"))

plot_df <- sample_qc_by_threshold %>%
  mutate(
    threshold_f = factor(threshold, levels = sort(unique(threshold))),
    sample = fct_inorder(sample),
    genotype = ifelse(is.na(genotype), "unknown", genotype),
    fraction = ifelse(is.na(fraction), "unknown", fraction),
    qc_group = ifelse(is.na(qc_group), "unknown", qc_group)
  )

p1 <- plot_df %>%
  ggplot(aes(x = threshold_f, y = domain_interpeak, group = sample, color = genotype)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ qc_group, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Domain / interpeak across BASIN thresholds",
    x = "Drop threshold",
    y = "Domain / interpeak"
  )

p2 <- plot_df %>%
  ggplot(aes(x = threshold_f, y = mean_anchor_interpeak, group = sample, color = genotype)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ qc_group, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Mean anchor / mean interpeak across BASIN thresholds",
    x = "Drop threshold",
    y = "Mean anchor / mean interpeak"
  )

ggsave(file.path(sec_dir, "threshold_domain_interpeak.png"), p1, width = 12, height = 6, dpi = 300)
ggsave(file.path(sec_dir, "threshold_anchor_interpeak.png"), p2, width = 12, height = 6, dpi = 300)

baseline_bar_df <- qc %>%
  filter(set_type == "baseline") %>%
  select(sample, qc_group, fraction, genotype, median_domain, median_anchor, median_interpeak) %>%
  pivot_longer(
    cols = c(median_domain, median_anchor, median_interpeak),
    names_to = "region_class",
    values_to = "median_signal"
  ) %>%
  mutate(
    region_class = recode(
      region_class,
      median_domain = "domain",
      median_anchor = "anchor",
      median_interpeak = "interpeak"
    ),
    sample = fct_inorder(sample)
  )

p3 <- baseline_bar_df %>%
  ggplot(aes(x = sample, y = median_signal, fill = region_class)) +
  geom_col(position = "dodge") +
  facet_grid(qc_group ~ genotype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Baseline median signal by sample and region class",
    x = "Sample",
    y = "Median signal"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file.path(sec_dir, "baseline_median_signal_grouped_barplot.png"), p3, width = 18, height = 8, dpi = 300)

baseline_ratio_bar_df <- qc %>%
  filter(set_type == "baseline") %>%
  select(sample, qc_group, fraction, genotype, domain_interpeak, anchor_interpeak, anchor_domain) %>%
  pivot_longer(
    cols = c(domain_interpeak, anchor_interpeak, anchor_domain),
    names_to = "ratio_name",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_name = recode(
      ratio_name,
      domain_interpeak = "domain/interpeak",
      anchor_interpeak = "anchor/interpeak",
      anchor_domain = "anchor/domain"
    ),
    sample = fct_inorder(sample)
  )

p_baseline_ratio <- baseline_ratio_bar_df %>%
  ggplot(aes(x = sample, y = ratio_value, fill = ratio_name)) +
  geom_col() +
  facet_grid(ratio_name ~ qc_group + genotype, scales = "free_y", space = "free_x") +
  theme_bw() +
  labs(
    title = "Baseline region-class ratios by sample",
    x = "Sample",
    y = "Ratio"
  ) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file.path(sec_dir, "baseline_ratio_grouped_barplot.png"),
       p_baseline_ratio, width = 18, height = 10, dpi = 300)

primary_bar_df <- qc %>%
  filter(set_type == "scaled", abs(threshold - primary_drop) < 1e-8) %>%
  select(sample, qc_group, fraction, genotype, median_domain, median_anchor, median_interpeak) %>%
  pivot_longer(
    cols = c(median_domain, median_anchor, median_interpeak),
    names_to = "region_class",
    values_to = "median_signal"
  ) %>%
  mutate(
    region_class = recode(
      region_class,
      median_domain = "domain",
      median_anchor = "anchor",
      median_interpeak = "interpeak"
    ),
    sample = fct_inorder(sample)
  )

p4 <- primary_bar_df %>%
  ggplot(aes(x = sample, y = median_signal, fill = region_class)) +
  geom_col(position = "dodge") +
  facet_grid(qc_group ~ genotype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = paste0("Primary-drop (", primary_drop, ") median signal by sample and region class"),
    x = "Sample",
    y = "Median signal"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file.path(sec_dir, "primary_drop_median_signal_grouped_barplot.png"), p4, width = 18, height = 8, dpi = 300)

ratio_bar_df <- qc %>%
  filter(set_type == "scaled", abs(threshold - primary_drop) < 1e-8) %>%
  select(sample, qc_group, fraction, genotype, domain_interpeak, anchor_interpeak, anchor_domain) %>%
  pivot_longer(
    cols = c(domain_interpeak, anchor_interpeak, anchor_domain),
    names_to = "ratio_name",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_name = recode(
      ratio_name,
      domain_interpeak = "domain/interpeak",
      anchor_interpeak = "anchor/interpeak",
      anchor_domain = "anchor/domain"
    ),
    sample = fct_inorder(sample)
  )

p5 <- ratio_bar_df %>%
  ggplot(aes(x = sample, y = ratio_value, fill = ratio_name)) +
  geom_col() +
  facet_grid(ratio_name ~ qc_group + genotype, scales = "free_y", space = "free_x") +
  theme_bw() +
  labs(
    title = paste0("Primary-drop (", primary_drop, ") region-class ratios by sample"),
    x = "Sample",
    y = "Ratio"
  ) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file.path(sec_dir, "primary_drop_ratio_grouped_barplot.png"),
       p5, width = 18, height = 10, dpi = 300)

retained_anchor_chrom_plot_df <- retained_anchor_chrom_fractions %>%
  filter(!is.na(threshold)) %>%
  left_join(
    manifest %>% select(sample, qc_group, fraction, genotype),
    by = "sample"
  ) %>%
  select(sample, qc_group, fraction, genotype, threshold,
         frac_retained_anchor_rows_autosomes,
         frac_retained_anchor_rows_X,
         frac_retained_anchor_rows_Y) %>%
  pivot_longer(
    cols = c(frac_retained_anchor_rows_autosomes,
             frac_retained_anchor_rows_X,
             frac_retained_anchor_rows_Y),
    names_to = "chrom_class",
    values_to = "fraction_of_retained_anchor_rows"
  ) %>%
  mutate(
    chrom_class = recode(
      chrom_class,
      frac_retained_anchor_rows_autosomes = "autosomes",
      frac_retained_anchor_rows_X = "X",
      frac_retained_anchor_rows_Y = "Y"
    ),
    threshold_f = factor(threshold, levels = sort(unique(threshold))),
    sample = fct_inorder(sample)
  )

p6 <- retained_anchor_chrom_plot_df %>%
  ggplot(aes(x = sample, y = fraction_of_retained_anchor_rows, fill = chrom_class)) +
  geom_col(position = "fill") +
  facet_grid(threshold_f ~ qc_group + genotype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Retained anchor chromosome composition by BASIN threshold",
    x = "Sample",
    y = "Fraction of retained anchor rows",
    fill = "Chromosome class"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file.path(sec_dir, "retained_anchor_chromosome_fraction_by_threshold.png"),
       p6, width = 18, height = 12, dpi = 300)
EOF

########################################
# STEP 5: RUN R SCRIPT
########################################

echo "Checking R syntax..."
Rscript -e "parse(file='${OUTDIR}/run_basin_qc.R')"

Rscript "${OUTDIR}/run_basin_qc.R" \
  "${OUTDIR}" \
  "${MANIFEST_TSV}" \
  "${SCALES_TSV}" \
  "${SIGNAL_SUMMARY_TSV}" \
  "${PRIMARY_DROP}" \
  "${PREFERRED_DROPS}" \
  "${MANUAL_QC_GROUPS}" \
  "${RETAINED_BARS_TSV}"

echo "QC complete: ${OUTDIR}"
