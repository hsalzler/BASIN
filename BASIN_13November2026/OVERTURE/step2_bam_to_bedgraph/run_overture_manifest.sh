#!/usr/bin/env bash
#SBATCH --job-name=overture_bam_to_bg
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/basin_univ_manifest_%j.out
#SBATCH --error=logs/basin_univ_manifest_%j.err

# ==============================================================================
# OVERTURE Step 2 — BAM to CPM bedGraph — BAM-ONLY
#
# Requires a TSV manifest with (at minimum) these columns (order doesn't matter):
#   sample_id   layout   bam   genome_id
# Optional:
#   source  (if present, must be BAM or blank)
#
# Example header:
#   sample_id  source  layout  bam  genome_id
#
# Run:
#   sbatch run_overture_manifest.sh templates/igg_samples.tsv config_overture.env
#
# Notes:
# - No SRA download
# - No FASTQ handling
# - No trimming/alignment
# - Produces CPM bedGraphs via bamCoverage
# - Submits SEACR job if SEACR sbatch script is present in current working dir
#   (or set SEACR_SCRIPT path in env file)
# - Optional: auto-submit consensus after SEACR (set RUN_CONSENSUS=true)
# ==============================================================================

set -euo pipefail
mkdir -p logs

# Recommended: ensure THREADS matches Slurm allocation
export THREADS="${SLURM_CPUS_PER_TASK:-8}"

MANIFEST="${1:-templates/igg_samples.tsv}"
ENVFILE="${2:-config_overture.env}"

die() {
  echo "ERROR: $*" >&2
  exit 2
}

norm_upper() { echo "${1:-}" | tr '[:lower:]' '[:upper:]'; }

# --- Load config ---
[[ -s "$ENVFILE" ]] || die "Config env file not found: $ENVFILE"
# shellcheck disable=SC1090
source "$ENVFILE"

# --- Required config vars (minimal) ---
: "${OUTDIR:?Config must define OUTDIR}"
: "${GENOME_ID:?Config must define GENOME_ID}"
: "${MIN_MAPQ:=0}"
: "${LINK_LOCAL_INPUTS:=true}"
: "${TMPBASE:=${SLURM_TMPDIR:-/tmp}}"

# Canonical BAR-discovery geometry is immutable. Legacy geometry variables are
# accepted only when harmless so old configs fail clearly rather than silently
# changing the measurement model.
if [[ "${BAMCOVERAGE_BINSIZE:-1}" != "1" ]]; then
  die "BAMCOVERAGE_BINSIZE is no longer configurable: OVERTURE BAR discovery requires binSize=1."
fi
if [[ -n "${BAMCOVERAGE_EXTRA:-}" ]]; then
  die "BAMCOVERAGE_EXTRA is no longer allowed. BAR discovery requires uncentered, unextended coverage; MIN_MAPQ is the supported read filter."
fi
if ! [[ "${MIN_MAPQ}" =~ ^[0-9]+$ ]]; then
  die "MIN_MAPQ must be a non-negative integer (got: ${MIN_MAPQ})."
fi

# Optional script paths
: "${SEACR_SCRIPT:=SEACR_IgG_peaks.sbatch}"
: "${CONSENSUS_SCRIPT:=SEACR_consensus.sbatch}"
: "${RUN_SEACR:=true}"
: "${RUN_CONSENSUS:=true}"

# --- Output subdirs ---
FASTQ_DIR="$OUTDIR/fastq"       # kept for compatibility; unused
BAM_DIR="$OUTDIR/bam"
BG_DIR="$OUTDIR/bedgraph_cpm"
SEACR_DIR="$OUTDIR/seacr"
CONS_DIR="$OUTDIR/consensus"
REPORT_DIR="$OUTDIR/report"

mkdir -p "$FASTQ_DIR" "$BAM_DIR" "$BG_DIR" "$SEACR_DIR/logs" "$CONS_DIR" "$REPORT_DIR"

cat > "$REPORT_DIR/coverage_contract.tsv" << EOF
key	value
contract_version	BASIN_CANONICAL_COVERAGE_V1
canonical_signal	aligned_read_coverage
canonical_bin_size	1
canonical_center_reads	false
canonical_extend_reads	false
normalization	CPM
min_mapq	${MIN_MAPQ}
EOF

# --- Modules (best effort) ---
if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  module load samtools 2>/dev/null || true
  module load deeptools 2>/dev/null || true
  module load bedtools 2>/dev/null || true
  module load SEACR/1.3 2>/dev/null || module load seacr/1.3 2>/dev/null || true
fi

echo "== OVERTURE Step 2 — BAM to CPM bedGraph =="
echo "Manifest: $MANIFEST"
echo "Config:   $ENVFILE"
echo "Genome:   $GENOME_ID"
echo "Outdir:   $OUTDIR"
echo "TMPBASE:  $TMPBASE"
echo "THREADS:  $THREADS"
echo "MIN_MAPQ: $MIN_MAPQ"
echo

[[ -s "$MANIFEST" ]] || die "Manifest not found or empty: $MANIFEST"

# ------------------------------------------------------------------------------
# Robust manifest parsing (HEADER-NAME based; order independent; fail-fast)
#
# We parse the manifest TSV by header *names* (not column positions) so reordering
# or adding/removing columns won’t break the pipeline.
#
# This block emits a clean, normalized 5-column stream:
#   sample_id \t source \t layout \t bam \t genome_id
#
# Key behaviors / hardening:
# - Skips the header row automatically.
# - Strips Windows CRLF line endings (removes trailing '\r') on every line.
# - Trims leading/trailing whitespace from all extracted fields.
# - Requires: sample_id, layout, bam header columns (exact names).
# - Optional: source, genome_id (may be absent or blank).
# - Skips rows with blank sample_id after trimming.
#
# Downstream, we freeze this stream to OUTDIR/manifest.used.tsv and run ALL work
# from that frozen file so the exact input universe is auditable and reproducible.
# ------------------------------------------------------------------------------


manifest_stream() {
  awk -F'\t' '
    BEGIN{OFS="\t"}
    NR==1{
      # strip CRLF on header line
      gsub(/\r$/, "", $0)

      for(i=1;i<=NF;i++){
        gsub(/^[ \t]+|[ \t]+$/,"",$i)
        h[$i]=i
      }
      if(!("sample_id" in h)) { print "MISSING_REQUIRED_COLUMN\tsample_id" > "/dev/stderr"; exit 3 }
      if(!("layout"    in h)) { print "MISSING_REQUIRED_COLUMN\tlayout"    > "/dev/stderr"; exit 3 }
      if(!("bam"       in h)) { print "MISSING_REQUIRED_COLUMN\tbam"       > "/dev/stderr"; exit 3 }
      next
    }
    NR>1{
      # strip CRLF on every line (VERY common culprit)
      gsub(/\r$/, "", $0)

      sid = $(h["sample_id"])
      src = ("source" in h) ? $(h["source"]) : ""
      lay = $(h["layout"])
      bam = $(h["bam"])
      gid = ("genome_id" in h) ? $(h["genome_id"]) : ""

      # trim outer whitespace
      gsub(/^[ \t]+|[ \t]+$/, "", sid)
      gsub(/^[ \t]+|[ \t]+$/, "", src)
      gsub(/^[ \t]+|[ \t]+$/, "", lay)
      gsub(/^[ \t]+|[ \t]+$/, "", bam)
      gsub(/^[ \t]+|[ \t]+$/, "", gid)

      if(sid=="") next
      print sid, src, lay, bam, gid
    }
  ' "$MANIFEST"
}

USED_MANIFEST="$OUTDIR/manifest.used.tsv"
manifest_stream | awk 'NF>0' > "$USED_MANIFEST"
echo "[INFO] Frozen manifest written: $USED_MANIFEST"

# Fail-fast: duplicate sample_id check
dups="$(awk -F'\t' '{c[$1]++} END{for(k in c) if(c[k]>1) print k}' "$USED_MANIFEST")"
if [[ -n "$dups" ]]; then
  echo "[ERROR] Duplicate sample_id(s) in manifest.used.tsv:" >&2
  echo "$dups" >&2
  exit 2
fi


# ------------------------------------------------------------------------------
# Main processing loop
# ------------------------------------------------------------------------------
while IFS=$'\t' read -r sample_id source layout bam genome_id; do
  # skip blank lines
  [[ -n "${sample_id:-}" ]] || continue

  # normalize
  source="$(norm_upper "${source:-BAM}")"
  layout="$(norm_upper "${layout:-}")"
  layout="$(echo "$layout" | tr -d '[:space:]')"  # trims tabs/spaces

  # allow blank genome_id -> assume GENOME_ID
  genome_id="${genome_id:-$GENOME_ID}"

  # validate genome
  if [[ "$genome_id" != "$GENOME_ID" ]]; then
    die "$sample_id has genome_id=$genome_id but config GENOME_ID=$GENOME_ID (one genome per run)."
  fi

  # validate source if provided
  if [[ -n "${source:-}" && "$source" != "BAM" ]]; then
    die "$sample_id has source=$source; BAM-only mode expects source=BAM (or blank)."
  fi

  # validate layout
  if [[ "$layout" != "PE" && "$layout" != "SE" ]]; then
    die "$sample_id has layout=$layout; must be PE or SE."
  fi

  # validate BAM
  [[ -n "${bam:-}" ]] || die "$sample_id bam path is empty."
  [[ -f "$bam" ]] || die "$sample_id bam not found: $bam"

  echo "==> Processing: $sample_id  (layout=$layout)"
  echo "    bam: $bam"

  # Link BAM into BAM_DIR (keeps downstream paths consistent)
  BAM_OUT="$BAM_DIR/${sample_id}.bam"
  if [[ -s "$BAM_OUT" ]]; then
    echo "    BAM already present: $BAM_OUT"
  else
    if [[ "${LINK_LOCAL_INPUTS}" == "true" ]]; then
      ln -sf "$(readlink -f "$bam")" "$BAM_OUT"
      echo "    Linked -> $BAM_OUT"
    else
      cp -av "$bam" "$BAM_OUT"
      echo "    Copied -> $BAM_OUT"
    fi
  fi

  # Ensure BAM index exists alongside BAM_OUT
  if [[ -f "${BAM_OUT}.bai" ]]; then
    : # ok
  else
    # If original BAM has an index, link it; else build one
    if [[ -f "${bam}.bai" ]]; then
      ln -sf "$(readlink -f "${bam}.bai")" "${BAM_OUT}.bai" || true
      echo "    Linked index -> ${BAM_OUT}.bai"
    elif [[ -f "${bam%.bam}.bai" ]]; then
      ln -sf "$(readlink -f "${bam%.bam}.bai")" "${BAM_OUT}.bai" || true
      echo "    Linked index -> ${BAM_OUT}.bai"
    else
      echo "    Index missing; creating index for $BAM_OUT"
      samtools index "$BAM_OUT"
    fi
  fi

  # bamCoverage CPM bedgraph
  BG_OUT="$BG_DIR/${sample_id}.CPM.bedgraph"
  if [[ -s "$BG_OUT" ]]; then
    echo "    Bedgraph exists, skipping: $BG_OUT"
  else
    echo "    bamCoverage -> CPM bedgraph"
tmp_bg="${BG_OUT}.tmp.$$"
bamCoverage -b "$BAM_OUT" -o "$tmp_bg" \
  --outFileFormat bedgraph \
  --binSize 1 \
  --minMappingQuality "$MIN_MAPQ" \
  --normalizeUsing CPM \
  --numberOfProcessors "$THREADS"

mv -f "$tmp_bg" "$BG_OUT"
  fi

  echo
done < "$USED_MANIFEST"

echo "[CHECK] Verifying bedgraph completeness vs manifest..."
missing=0
while IFS=$'\t' read -r sid _ _ _ _; do
  [[ -n "$sid" ]] || continue
  bg="$BG_DIR/${sid}.CPM.bedgraph"
  if [[ ! -s "$bg" ]]; then
    echo "[MISSING] $sid -> expected $bg" >&2
    missing=1
  fi
done < "$USED_MANIFEST"

if [[ "$missing" -ne 0 ]]; then
  die "One or more manifest samples did not produce CPM bedgraphs. Fix upstream before SEACR."
fi

echo "[OK] All manifest samples produced CPM bedgraphs."

# ------------------------------------------------------------------------------
# Submit SEACR
# ------------------------------------------------------------------------------
if [[ "${RUN_SEACR}" != "true" ]]; then
  echo "RUN_SEACR=false; skipping SEACR submission."
  echo "CPM bedgraphs are ready at: $BG_DIR"
  exit 0
fi

if [[ ! -f "$SEACR_SCRIPT" ]]; then
  echo "NOTE: SEACR sbatch script not found:"
  echo "  $SEACR_SCRIPT"
  echo
  echo "CPM bedgraphs are ready at:"
  echo "  $BG_DIR"
  echo
  echo "To run SEACR manually (example):"
  echo "  sbatch --output=\"$SEACR_DIR/logs/seacr_%j.out\" --error=\"$SEACR_DIR/logs/seacr_%j.err\" \\"
  echo "    --export=ALL,SEARCH_DIR=\"$BG_DIR\",MATCH_TOKEN=\"\" \\"
  echo "    \"$SEACR_SCRIPT\""
  exit 0
fi

echo "==> Submitting SEACR via: $SEACR_SCRIPT"

SEACR_JOBID="$(
  sbatch --output="$SEACR_DIR/logs/seacr_%j.out" --error="$SEACR_DIR/logs/seacr_%j.err" \
    --export=ALL,SEARCH_DIR="$BG_DIR",MATCH_TOKEN="",BASE_RELAXED_DIR="$SEACR_DIR/SEACR_peaks_relaxed_CPM",BASE_STRINGENT_DIR="$SEACR_DIR/SEACR_peaks_stringent_CPM",BASE_LOG_DIR="$SEACR_DIR/seacr_logs" \
    "$SEACR_SCRIPT" | awk '{print $4}'
)"

echo "Submitted SEACR job: $SEACR_JOBID"
echo

# ------------------------------------------------------------------------------
# Optional: Submit consensus after SEACR completes
# ------------------------------------------------------------------------------
if [[ "${RUN_CONSENSUS:-false}" == "true" ]]; then
  [[ -f "$CONSENSUS_SCRIPT" ]] || die "Consensus script not found: $CONSENSUS_SCRIPT"

  echo "==> Submitting consensus after SEACR completes: $CONSENSUS_SCRIPT"
  CONS_JOBID="$(
    sbatch --dependency="afterok:${SEACR_JOBID}" \
      --output="$CONS_DIR/consensus_%j.out" --error="$CONS_DIR/consensus_%j.err" \
      --export=ALL,SEACR_DIR="$SEACR_DIR",CONS_OUTDIR="$CONS_DIR",PATTERN="*.bed",MODE="${CONS_MODE:-merged}",MERGE_DISTANCE="${CONS_MERGE_DISTANCE:-0}" \
      "$CONSENSUS_SCRIPT" | awk '{print $4}'
  )"
  echo "Submitted consensus job: $CONS_JOBID"
fi
