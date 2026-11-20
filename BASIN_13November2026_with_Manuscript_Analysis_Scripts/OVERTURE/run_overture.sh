#!/usr/bin/env bash
#SBATCH --job-name=overture_driver
#SBATCH --output=overture_driver.%j.out
#SBATCH --error=overture_driver.%j.err
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00

# ==============================================================================
# run_overture.sh — OVERTURE master wrapper
#
# Sets shared configuration in one place and submits all four OVERTURE steps
# as a Slurm job chain. Each step waits for the previous to succeed before
# starting (afterok dependencies).
#
# Usage:
#   sbatch run_overture.sh           # dry run — prints the job chain, submits nothing
#   sbatch run_overture.sh --submit  # submit all jobs to Slurm
#
# Before running:
#   1. Edit the CONFIG section below.
#   2. For universal panel:    verify igg_samples_used/igg_samples_universal.tsv
#                              has your local BAM paths filled in.
#      For experiment-specific: point MANIFEST at your own igg_samples.tsv.
#   3. If starting from FASTQs, set RUN_ALIGNMENT=true and fill in FASTQ_DIR
#      and BT2_INDEX. Otherwise set RUN_ALIGNMENT=false and point MANIFEST
#      directly at your pre-aligned BAMs.
#
# Note on the alignment step:
#   align_igg_bowtie2.sh performs stringent Cutadapt trimming and end-to-end Bowtie2 alignment.
#   If you have already aligned your IgG files — or prefer your own alignment
#   pipeline — set RUN_ALIGNMENT=false. The only requirement is that BAMs are
#   sorted, indexed, and aligned with consistent parameters across all samples.
#   For CUT&RUN/CUT&Tag, the retained workflow uses Bowtie2 settings (--no-mixed --no-discordant,
#   -I 10 -X 700) are strongly recommended.
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIG — edit everything in this section
# ==============================================================================

# --- Path to this OVERTURE directory (REQUIRED) ---
# Set this to the absolute path of your OVERTURE folder.
# Do not use relative paths — Slurm resolves paths from the spool directory,
# not from where you submitted the job.
SCRIPT_DIR="/path/to/OVERTURE"

# --- Genome ---
GENOME_ID="dm6"

# --- Run mode ---
# "universal"  : use the 28-dataset public compendium (igg_samples_universal.tsv)
# "experiment" : use your own IgG controls
MODE="universal"

# --- Output root ---
# All step outputs land under this directory
OUTROOT="$PWD/overture_run_${GENOME_ID}_${MODE}_$(date +%Y%m%d_%H%M%S)"

# --- Manifest ---
# For universal mode, use the distributed manifest (update bam paths first).
# For experiment-specific mode, point this at your own igg_samples.tsv.
MANIFEST="$PWD/igg_samples_used/igg_samples_universal.tsv"
# MANIFEST="$PWD/step2_bam_to_bedgraph/templates/igg_samples.tsv"  # experiment-specific

# --- Step 1: Alignment ---
# Set to true only if starting from FASTQs and want this script to align them.
# If false, Step 1 is skipped entirely and Step 2 reads BAMs from MANIFEST.
# Step 1 produces full sorted BAMs at OUTROOT/alignment/bam/.
# Step 2 reads BAM paths from MANIFEST; update the manifest before running.
RUN_ALIGNMENT=false
FASTQ_DIR="/path/to/your/fastqs"
BT2_INDEX="/path/to/bowtie2/index/dm6"

# --- Step 2: BAR-discovery coverage ---
# Geometry is fixed at binSize=1, uncentered, unextended. MAPQ is configurable.
MIN_MAPQ=0

# --- Step 3: SEACR grid ---
# Space-separated FDR thresholds. The default covers a broad range for discovery.
# Typical useful range: 0.0005 – 0.01
FDR_LIST="0.01 0.05 0.005 0.0005"

# --- Step 4: Panel evaluation ---
# Directory of CPM bigWigs used to benchmark candidate panels.
# These should be from a validation experiment with biological replicates and
# multiple genotypes (so replicate CV and genotype shift can be calculated).
# For the universal dm6 panel, this was a 6-mark L1 CUT&TAG dataset.
EVAL_BW_DIR="/path/to/CPM_bigwigs_for_evaluation"

# Reference genotype for genotype-shift calculation
REF_GENOTYPE="yw"

# Comma-separated histone marks present in BW_DIR to evaluate against
EVAL_MARKS="H3K36me2,H3K36me3,H3K27ac,panH4,RNAPIIS5P,RNAPIIS2P"

# --- Step 4: Panel evaluation thresholds ---
EVAL_MIN_REGIONS=200       # minimum regions per panel to be evaluated
EVAL_MAX_REGIONS=5000      # maximum regions; raise for mammalian genomes (e.g. 10000)
EVAL_MAX_CV=0.12           # maximum mean replicate CV; lower = more reproducible
EVAL_DROP_FRAC="median_drop20"  # drop fraction used for CV and shift calculations

# ==============================================================================
# END CONFIG
# ==============================================================================

SUBMIT=false
[[ "${1:-}" == "--submit" ]] && SUBMIT=true

echo "================================================================"
echo "  OVERTURE — BAR discovery"
echo "  Mode     : $MODE"
echo "  Genome   : $GENOME_ID"
echo "  Outroot  : $OUTROOT"
echo "  Manifest : $MANIFEST"
echo "  Submit   : $SUBMIT"
echo "================================================================"
echo

mkdir -p "$OUTROOT/logs"
mkdir -p "$OUTROOT/bedgraph_run"
mkdir -p "$OUTROOT/panel_eval_out"

# Helper: submit or dry-run
submit_job() {
  local label="$1"; shift
  if [[ "$SUBMIT" == true ]]; then
    local jid
    jid=$(sbatch --parsable "$@")
    echo "[SUBMITTED] $label → job $jid" >&2
    echo "$jid"
  else
    echo "[DRY RUN]  $label" >&2
    echo "  Would run: sbatch $*" >&2
    echo "dry_run"
  fi
}

JID_ALIGN="none"
JID_BG="none"
JID_GRID="none"
JID_EVAL="none"

# ------------------------------------------------------------------------------
# Step 1 — Alignment (optional)
# ------------------------------------------------------------------------------
if [[ "$RUN_ALIGNMENT" == true ]]; then
  echo "--- Step 1: Alignment ---"
  [[ -d "$FASTQ_DIR" ]]      || { echo "ERROR: FASTQ_DIR not found: $FASTQ_DIR"; exit 1; }
  [[ -f "${BT2_INDEX}.1.bt2" || -f "${BT2_INDEX}.1.bt2l" ]] \
    || echo "WARN: Bowtie2 index not found at $BT2_INDEX — verify before submitting."

  JID_ALIGN=$(submit_job "Step 1: align" \
    --job-name=overture_align \
    --output="$OUTROOT/logs/step1_align_%j.out" \
    --error="$OUTROOT/logs/step1_align_%j.err" \
    "$SCRIPT_DIR/step1_align/align_igg_bowtie2.sh" \
      --proj "$OUTROOT" \
      --fastq_dir "$FASTQ_DIR" \
      --bt2_index_prefix "$BT2_INDEX"
  )
  echo
else
  echo "--- Step 1: Alignment skipped (RUN_ALIGNMENT=false) ---"
  echo "  Reads BAMs directly from manifest: $MANIFEST"
  echo
fi

# ------------------------------------------------------------------------------
# Step 2 — BAM → CPM bedGraphs
# ------------------------------------------------------------------------------
echo "--- Step 2: BAM → CPM bedGraphs ---"

# Write a config env for this run
CONFIG_ENV="$OUTROOT/config_this_run.env"
cat > "$CONFIG_ENV" << EOF
GENOME_ID="${GENOME_ID}"
RUN_TAG="${MODE}"
OUTDIR="${OUTROOT}/bedgraph_run"
MIN_MAPQ="${MIN_MAPQ}"
LINK_LOCAL_INPUTS=true
RUN_SEACR=false
RUN_CONSENSUS=false
EOF

STEP2_ARGS=(
  --job-name=overture_bam_to_bg
  --output="$OUTROOT/logs/step2_bam_to_bg_%j.out"
  --error="$OUTROOT/logs/step2_bam_to_bg_%j.err"
)
[[ "$JID_ALIGN" != "none" && "$JID_ALIGN" != "dry_run" ]] && \
  STEP2_ARGS+=(--dependency="afterok:${JID_ALIGN}")

JID_BG=$(submit_job "Step 2: BAM → bedGraph" \
  "${STEP2_ARGS[@]}" \
  "$SCRIPT_DIR/step2_bam_to_bedgraph/run_overture_manifest.sh" \
    "$MANIFEST" "$CONFIG_ENV"
)
echo

# ------------------------------------------------------------------------------
# Step 3 — SEACR grid + consensus
# ------------------------------------------------------------------------------
echo "--- Step 3: SEACR grid + consensus ---"

STEP3_ARGS=(
  --job-name=overture_grid
  --output="$OUTROOT/logs/step3_grid_%j.out"
  --error="$OUTROOT/logs/step3_grid_%j.err"
)
[[ "$JID_BG" != "none" && "$JID_BG" != "dry_run" ]] && \
  STEP3_ARGS+=(--dependency="afterok:${JID_BG}")

JID_GRID=$(submit_job "Step 3: SEACR grid" \
  "${STEP3_ARGS[@]}" \
  --export=ALL,FDR_LIST="${FDR_LIST}",\
SEACR_SCRIPT="$SCRIPT_DIR/step3_seacr_grid/SEACR_IgG_peaks.sbatch",\
CONS_SCRIPT="$SCRIPT_DIR/step3_seacr_grid/SEACR_consensus.sbatch" \
  "$SCRIPT_DIR/step3_seacr_grid/run_seacr_grid.sh" \
    "${OUTROOT}/bedgraph_run"
)
echo

# ------------------------------------------------------------------------------
# Step 4 — Panel evaluation
#
# The grid driver (Step 3) spawns child jobs and exits immediately, so we
# cannot depend directly on JID_GRID for Step 4 — it would run before the
# SEACR and consensus work is done.
#
# Instead: the grid script writes the final summary job ID to
# bedgraph_run/summary_jid.txt. A lightweight launcher job (afterok:JID_GRID)
# reads that file and submits the real eval job with the correct dependency.
# ------------------------------------------------------------------------------
echo "--- Step 4: Panel evaluation (via launcher) ---"

[[ "$EVAL_BW_DIR" == "/path/to/CPM_bigwigs_for_evaluation" ]] && \
  echo "WARN: EVAL_BW_DIR is still the placeholder — update before submitting Step 4."

SUMMARY_JID_FILE="${OUTROOT}/bedgraph_run/summary_jid.txt"

LAUNCHER_SCRIPT="$OUTROOT/logs/launch_step4.sh"
cat > "$LAUNCHER_SCRIPT" << LAUNCHER
#!/usr/bin/env bash
#SBATCH --job-name=overture_step4_launcher
#SBATCH --output=${OUTROOT}/logs/step4_launcher_%j.out
#SBATCH --error=${OUTROOT}/logs/step4_launcher_%j.err
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

set -euo pipefail

# Wait for summary_jid.txt to be written by the grid script
for i in \$(seq 1 24); do
  [[ -f "${SUMMARY_JID_FILE}" ]] && break
  echo "Waiting for summary_jid.txt... (\$i/24)"
  sleep 5
done

[[ -f "${SUMMARY_JID_FILE}" ]] || {
  echo "ERROR: summary_jid.txt never appeared at ${SUMMARY_JID_FILE}"
  exit 1
}

SUMMARY_JID=\$(cat "${SUMMARY_JID_FILE}")
echo "Grid summary job ID: \$SUMMARY_JID"
echo "Submitting Step 4 with dependency afterok:\$SUMMARY_JID"

sbatch --dependency="afterok:\${SUMMARY_JID}" \
  --job-name=overture_panel_eval \
  --output="${OUTROOT}/logs/step4_eval_%j.out" \
  --error="${OUTROOT}/logs/step4_eval_%j.err" \
  --time=12:00:00 \
  --mem=12G \
  --cpus-per-task=1 \
  -p general \
  --export=ALL \
  --wrap="export BW_DIR='${EVAL_BW_DIR}'; \
export REGION_DIR='${OUTROOT}/bedgraph_run/consensus_grid_collected'; \
export OUTDIR='${OUTROOT}/panel_eval_out'; \
export REF_GENOTYPE='${REF_GENOTYPE}'; \
export MARKS='${EVAL_MARKS}'; \
export MIN_REGIONS='${EVAL_MIN_REGIONS}'; \
export EVAL_MIN_REGIONS='${EVAL_MIN_REGIONS}'; \
export EVAL_MAX_REGIONS='${EVAL_MAX_REGIONS}'; \
export EVAL_MAX_CV='${EVAL_MAX_CV}'; \
export EVAL_DROP_FRAC='${EVAL_DROP_FRAC}'; \
bash '${SCRIPT_DIR}/step4_panel_evaluation/eval_panels.sbatch'"

echo "Step 4 submitted."
LAUNCHER

chmod +x "$LAUNCHER_SCRIPT"

STEP4_LAUNCHER_ARGS=(
  --job-name=overture_step4_launcher
  --output="$OUTROOT/logs/step4_launcher_%j.out"
  --error="$OUTROOT/logs/step4_launcher_%j.err"
)
[[ "$JID_GRID" != "none" && "$JID_GRID" != "dry_run" ]] && \
  STEP4_LAUNCHER_ARGS+=(--dependency="afterok:${JID_GRID}")

JID_EVAL=$(submit_job "Step 4: panel eval launcher" \
  "${STEP4_LAUNCHER_ARGS[@]}" \
  "$LAUNCHER_SCRIPT"
)
echo

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo "================================================================"
if [[ "$SUBMIT" == true ]]; then
  echo "  Job chain submitted."
  echo ""
  [[ "$RUN_ALIGNMENT" == true ]] && echo "  Step 1 (align)        : job $JID_ALIGN"
  echo "  Step 2 (bedgraph)     : job $JID_BG"
  echo "  Step 3 (grid)         : job $JID_GRID"
  echo "  Step 4 launcher       : job $JID_EVAL"
  echo "  Step 4 (eval)         : submitted by launcher after grid completes"
  echo ""
  echo "  Logs : $OUTROOT/logs/"
  echo "  Final candidate panels : $OUTROOT/bedgraph_run/consensus_grid_collected/"
  echo "  Evaluation results     : $OUTROOT/panel_eval_out/"
  echo ""
  echo "  Note: Step 4 is submitted dynamically by the launcher job once the"
  echo "  grid summary job ID is known. Check squeue for 'overture_panel_eval'."
else
  echo "  DRY RUN complete — no jobs submitted."
  echo "  Review the config above, then run with --submit to launch."
fi
echo "================================================================"
