#!/usr/bin/env bash
#SBATCH --job-name=overture_grid
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=seacr_grid_driver.%j.out
#SBATCH --error=seacr_grid_driver.%j.err
set -euo pipefail

OUTDIR="${1:?Usage: $0 /path/to/OUTDIR}"

# ==============================================================================
# run_seacr_consensus_grid.sh
#
# Purpose
# -------
# Run a small “grid search” over SEACR FDR fractions using the CPM bedGraphs you
# already generated in an upstream BASIN run folder (OUTDIR).
#
# For each FDR fraction in FDR_LIST:
#   1) Submit ONE SEACR job (your existing SEACR sbatch runs BOTH relaxed+stringent)
#      and FORCE its outputs into:
#         OUTDIR/seacr_grid/fdr0p01/
#         OUTDIR/seacr_grid/fdr0p05/
#         OUTDIR/seacr_grid/fdr0p005/
#         OUTDIR/seacr_grid/fdr0p0005/
#
#   2) Submit ONE consensus job (dependent on SEACR job completion) using your
#      updated consensus script, which runs consensus separately for:
#         relaxed peaks  -> OUTDIR/consensus_grid/fdr.../relaxed_IgG/
#         stringent peaks-> OUTDIR/consensus_grid/fdr.../stringent_IgG/
#
#   3) Submit ONE summary+collector job (dependent on all consensus jobs) that:
#      - Finds all common_{80,90,100}pct.bed files across the grid
#            (includes filtered variants like *_filter_away_lt50bp.bed)
#      - Computes useful region stats (counts, length distribution, total bp, N50, etc.)
#      - Copies/collects each consensus BED into a single folder with standardized names
#      - Writes a single TSV summary table
#
# Assumptions
# -----------
# Your upstream run created:
#   BG_DIR = OUTDIR/bedgraph_cpm
# and it contains only the bedgraphs you want SEACR to use.
#
# Run instructions
# ----------------
# 1) Make executable:
#    chmod +x slurm/run_seacr_consensus_grid.sh
#
# 2) Run (from anywhere):
#    sbatch run_seacr_consensus_grid.sh /path/to/your/OUTDIR
#
# Optional overrides (environment variables):
# -------------------------------------------
#   FDR_LIST="0.01 0.05 0.005 0.0005"
#   MATCH_TOKEN=""              # recommended (no filename filtering) for IgG-only grid
#   THREADS=8                   # how many threads SEACR wrapper should use
#   SEACR_SCRIPT=/path/to/SEACR_IgG_peaks.sbatch
#   CONS_SCRIPT=/path/to/SEACR_consensus.sbatch
#
##   # If your consensus script supports length-filtered copies:
#   MAKE_LEN_FILTERED=yes|no
#   MIN_LEN=50
#   LEN_FILTER_SUFFIX=_filter_away_lt50bp
#
# Notes on safety / avoiding overwrites
# -------------------------------------
# This script writes each FDR run into its own subfolder under OUTDIR, so runs do
# not overwrite one another. If you re-run this script, it will submit jobs again
# into those same subfolders; that could overwrite results for that FDR. If you
# want “never overwrite”, see the SKIP_IF_EXISTS option below.
# ==============================================================================

# Fixed by your upstream convention:
BG_DIR="${BG_DIR:-$OUTDIR/bedgraph_cpm}"

# Defaults you asked for:
FDR_LIST="${FDR_LIST:-0.01 0.05 0.005 0.0005}"

# Match token (blank means “don’t filter by filename”)
MATCH_TOKEN="${MATCH_TOKEN-}"

# Threads to pass into SEACR wrapper (if it uses THREADS)
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

# Scripts (default: expect them in current directory)
SEACR_SCRIPT="${SEACR_SCRIPT:-SEACR_IgG_peaks.sbatch}"
CONS_SCRIPT="${CONS_SCRIPT:-SEACR_consensus.sbatch}"

# If true, skip submitting SEACR+consensus for an FDR if consensus outputs already exist
SKIP_IF_EXISTS="${SKIP_IF_EXISTS:-false}"

die() { echo "[GRID] ERROR: $*" >&2; exit 2; }

[[ -d "$OUTDIR" ]] || die "OUTDIR not found: $OUTDIR"
[[ -d "$BG_DIR"  ]] || die "BG_DIR not found: $BG_DIR (expected OUTDIR/bedgraph_cpm)"
[[ -f "$SEACR_SCRIPT" ]] || die "SEACR script not found: $SEACR_SCRIPT"
[[ -f "$CONS_SCRIPT"  ]] || die "Consensus script not found: $CONS_SCRIPT"

mkdir -p \
  "$OUTDIR/seacr_grid" \
  "$OUTDIR/consensus_grid" \
  "$OUTDIR/consensus_grid_collected" \
  "$OUTDIR/consensus_grid_logs"

echo "== SEACR+CONSENSUS GRID =="
echo "OUTDIR      : $OUTDIR"
echo "BG_DIR      : $BG_DIR"
echo "FDR_LIST    : $FDR_LIST"
echo "MATCH_TOKEN : '${MATCH_TOKEN}'"
echo "THREADS     : $THREADS"
echo "SEACR_SCRIPT: $SEACR_SCRIPT"
echo "CONS_SCRIPT : $CONS_SCRIPT"
echo "SKIP_IF_EXISTS: $SKIP_IF_EXISTS"
echo

echo "Bedgraphs present:"
ls -1 "$BG_DIR"/*.bedgraph 2>/dev/null | wc -l || true
if [[ -n "$MATCH_TOKEN" ]]; then
  echo "Bedgraphs matching token '$MATCH_TOKEN':"
  ls -1 "$BG_DIR"/*.bedgraph 2>/dev/null | grep -i -- "$MATCH_TOKEN" | wc -l || true
fi
echo

# ----------------------------------------------------------------------
# Helper: make a safe tag like fdr0p005 from 0.005
# ----------------------------------------------------------------------
fdr_tag() {
  local f="$1"
  echo "fdr$(echo "$f" | sed 's/+//g; s/\./p/g')"
}

SEACR_JOBIDS=()
CONS_JOBIDS=()

for fdr in $FDR_LIST; do
  tag="$(fdr_tag "$fdr")"
  seacr_root="$OUTDIR/seacr_grid/$tag"
  cons_root="$OUTDIR/consensus_grid/$tag"

  mkdir -p "$seacr_root/logs" "$seacr_root/seacr_logs" "$cons_root/logs"

  # Optional skip logic (if consensus already exists)
  if [[ "$SKIP_IF_EXISTS" == "true" ]]; then
    if [[ -s "$cons_root/relaxed_IgG/common_80pct.bed" && -s "$cons_root/stringent_IgG/common_80pct.bed" ]]; then
      echo "[SKIP] $tag: consensus outputs already exist in $cons_root"
      continue
    fi
  fi

  echo "==> [$tag] Submitting SEACR (relaxed+stringent)"
  echo "    SEACR outputs -> $seacr_root"
  echo "    Consensus     -> $cons_root"

  # NOTE: This assumes your SEACR wrapper honors FDR_FRAC.
  seacr_jid="$(
    sbatch --parsable \
      --output="$seacr_root/logs/seacr_%j.out" --error="$seacr_root/logs/seacr_%j.err" \
      --export=ALL,SEARCH_DIR="$BG_DIR",MATCH_TOKEN="",THREADS="$THREADS",FDR_FRAC="$fdr",\
BASE_RELAXED_DIR="$seacr_root/SEACR_peaks_relaxed_CPM",\
BASE_STRINGENT_DIR="$seacr_root/SEACR_peaks_stringent_CPM",\
BASE_LOG_DIR="$seacr_root/seacr_logs" \
      "$SEACR_SCRIPT"
  )"
  echo "    SEACR jobid: $seacr_jid"
  SEACR_JOBIDS+=("$seacr_jid")

  echo "==> [$tag] Submitting CONSENSUS after SEACR (afterok:$seacr_jid)"
  cons_jid="$(
    sbatch --parsable \
      --dependency="afterok:${seacr_jid}" \
      --output="$cons_root/logs/consensus_%j.out" --error="$cons_root/logs/consensus_%j.err" \
      --export=ALL,SEACR_DIR="$seacr_root",CONS_OUTDIR="$cons_root",PATTERN="*.bed",MODE="merged",MERGE_DISTANCE="0" \
      "$CONS_SCRIPT"
  )"
  echo "    CONS jobid: $cons_jid"
  CONS_JOBIDS+=("$cons_jid")
  echo
done

if [[ "${#CONS_JOBIDS[@]}" -eq 0 ]]; then
  echo "[GRID] No consensus jobs submitted (everything skipped?)."
  echo "[GRID] Nothing to summarize."
  exit 0
fi

# ----------------------------------------------------------------------
# Create a one-off summary sbatch inside the run folder and submit it.
# Now filtered-aware: collects common_*pct*.bed (includes suffixes).
# ----------------------------------------------------------------------
SUMMARY_SCRIPT_PATH="$OUTDIR/consensus_grid_logs/consensus_grid_summary.sbatch"
cat > "$SUMMARY_SCRIPT_PATH" << 'SBATCH'
#!/usr/bin/env bash
#SBATCH --job-name=cons_grid_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=consensus_grid_summary.%j.out
#SBATCH --error=consensus_grid_summary.%j.err
set -euo pipefail

CONS_GRID_DIR="${CONS_GRID_DIR:?export CONS_GRID_DIR=/path/to/consensus_grid}"
COLLECT_DIR="${COLLECT_DIR:?export COLLECT_DIR=/path/to/consensus_grid_collected}"
SUMMARY_OUT="${SUMMARY_OUT:?export SUMMARY_OUT=/path/to/consensus_grid_summary.tsv}"

if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  module load python 2>/dev/null || true
fi
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found" >&2; exit 1; }

mkdir -p "$COLLECT_DIR"
mkdir -p "$(dirname "$SUMMARY_OUT")"

python3 - << 'PY'
import os, re, sys, math, statistics, shutil
from pathlib import Path

cons_grid_dir = Path(os.environ["CONS_GRID_DIR"]).resolve()
collect_dir   = Path(os.environ["COLLECT_DIR"]).resolve()
summary_out   = Path(os.environ["SUMMARY_OUT"]).resolve()

def n50(lengths):
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2
    s = 0
    for L in sorted(lengths, reverse=True):
        s += L
        if s >= half:
            return L
    return 0

def pctile(lengths, p):
    if not lengths:
        return 0
    xs = sorted(lengths)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return xs[int(k)]
    d0 = xs[f] * (c - k)
    d1 = xs[c] * (k - f)
    return int(round(d0 + d1))

# Discover consensus beds (unfiltered + filtered variants)
# e.g. consensus_grid/fdr*/(relaxed_IgG|stringent_IgG)/common_80pct.bed
#      consensus_grid/fdr*/.../common_80pct_filter_away_lt50bp.bed
bed_files = sorted(cons_grid_dir.glob("fdr*/**/common_*pct*.bed"))
if not bed_files:
    print(f"ERROR: no consensus BEDs found under {cons_grid_dir}", file=sys.stderr)
    sys.exit(1)

rows = []
for bed in bed_files:
    parts = bed.parts

    fdr = next((p for p in parts if p.startswith("fdr")), "fdrNA")

    if "relaxed_IgG" in parts:
        stringency = "relaxed"
    elif "stringent_IgG" in parts:
        stringency = "stringent"
    else:
        stringency = "unknown"

    # Parse: common_{pct}pct{variant}.bed
    # variant is "" or something like "_filter_away_lt50bp"
    m = re.search(r"^common_(\d+)pct(.*)\.bed$", bed.name)
    if not m:
        continue
    pct = int(m.group(1))
    variant = m.group(2) or ""  # includes leading underscore if present

    lengths = []
    n = 0
    total_bp = 0

    with bed.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                start = int(cols[1]); end = int(cols[2])
            except ValueError:
                continue
            L = max(0, end - start)
            lengths.append(L)
            n += 1
            total_bp += L

    if n == 0:
        meanL = medL = sdL = minL = maxL = p90 = p95 = n50L = 0
    else:
        meanL = total_bp / n
        medL  = statistics.median(lengths)
        sdL   = statistics.pstdev(lengths) if n > 1 else 0.0
        minL  = min(lengths)
        maxL  = max(lengths)
        p90   = pctile(lengths, 0.90)
        p95   = pctile(lengths, 0.95)
        n50L  = n50(lengths)

    # Deterministic collected name that preserves variant
    out_name = f"IgG_consensus_{fdr}_{stringency}_{pct}pct{variant}.bed"
    out_path = collect_dir / out_name
    shutil.copyfile(bed, out_path)

    # Also emit explicit chr/nochr versions to prevent downstream name mismatches.
    # These are line-for-line identical except for column 1.
    out_chr_path = collect_dir / out_name.replace(".bed", ".chr.bed")
    out_nochr_path = collect_dir / out_name.replace(".bed", ".nochr.bed")

    def _write_chr_variants(in_bed: Path, out_chr: Path, out_nochr: Path) -> None:
        with in_bed.open() as r, out_chr.open("w") as w_chr, out_nochr.open("w") as w_no:
            for line in r:
                if not line.strip() or line.startswith("#"):
                    # Preserve comment/blank lines as-is
                    w_chr.write(line)
                    w_no.write(line)
                    continue
                parts = line.rstrip("\n").split("\t")
                if not parts:
                    w_chr.write(line)
                    w_no.write(line)
                    continue
                chrom = parts[0]
                # chr version
                if chrom.startswith("chr"):
                    chrom_chr = chrom
                else:
                    chrom_chr = "chr" + chrom
                # no-chr version
                if chrom.startswith("chr"):
                    chrom_no = chrom[3:]
                else:
                    chrom_no = chrom
                parts_chr = [chrom_chr] + parts[1:]
                parts_no  = [chrom_no]  + parts[1:]
                w_chr.write("\t".join(parts_chr) + "\n")
                w_no.write("\t".join(parts_no) + "\n")

    _write_chr_variants(out_path, out_chr_path, out_nochr_path)

    rows.append({
        "fdr": fdr,
        "stringency": stringency,
        "pct": pct,
        "variant": variant.lstrip("_") if variant else "",
        "n_regions": n,
        "total_bp": total_bp,
        "mean_len": round(meanL, 3) if n else 0,
        "median_len": medL,
        "sd_len": round(sdL, 3) if n else 0,
        "min_len": minL,
        "p90_len": p90,
        "p95_len": p95,
        "n50_len": n50L,
        "max_len": maxL,
        "source_bed": str(bed),
        "collected_bed": str(out_path),
        "collected_chr_bed": str(out_chr_path),
        "collected_nochr_bed": str(out_nochr_path),
    })

# Also collect per-run summaries if present (optional but helpful)
# consensus_grid/fdr*/relaxed_IgG/consensus_len_summary*.tsv etc.
summary_files = sorted(cons_grid_dir.glob("fdr*/**/consensus_len_summary*.tsv"))
for s in summary_files:
    # Preserve uniqueness by prefixing path-derived tags
    parts = s.parts
    fdr = next((p for p in parts if p.startswith("fdr")), "fdrNA")
    if "relaxed_IgG" in parts:
        stringency = "relaxed"
    elif "stringent_IgG" in parts:
        stringency = "stringent"
    else:
        stringency = "unknown"
    out_name = f"IgG_{fdr}_{stringency}__{s.name}"
    out_path = collect_dir / out_name
    try:
        shutil.copyfile(s, out_path)
    except Exception:
        pass

rows.sort(key=lambda r: (r["fdr"], r["stringency"], r["pct"], r["variant"]))

cols = [
    "fdr","stringency","pct","variant",
    "n_regions","total_bp","mean_len","median_len","sd_len",
    "min_len","p90_len","p95_len","n50_len","max_len",
    "source_bed","collected_bed","collected_chr_bed","collected_nochr_bed"
]
with summary_out.open("w") as w:
    w.write("\t".join(cols) + "\n")
    for r in rows:
        w.write("\t".join(str(r[c]) for c in cols) + "\n")

print(f"Wrote summary: {summary_out}")
print(f"Collected beds: {collect_dir}")
print(f"Also collected {len(summary_files)} per-run summary TSVs (if any).")
PY
SBATCH
chmod +x "$SUMMARY_SCRIPT_PATH"

deps="$(IFS=,; echo "${CONS_JOBIDS[*]}")"

echo "==> Submitting SUMMARY after all consensus jobs succeed (afterok:$deps)"
summary_jid="$(
  sbatch --parsable \
    --dependency="afterok:${deps}" \
    --output="$OUTDIR/consensus_grid_logs/summary_%j.out" --error="$OUTDIR/consensus_grid_logs/summary_%j.err" \
    --export=ALL,CONS_GRID_DIR="$OUTDIR/consensus_grid",COLLECT_DIR="$OUTDIR/consensus_grid_collected",SUMMARY_OUT="$OUTDIR/consensus_grid_summary.tsv" \
    "$SUMMARY_SCRIPT_PATH"
)"
echo "SUMMARY jobid: $summary_jid"

# Write the summary job ID to a file so the OVERTURE wrapper can
# depend on it for Step 4 — the wrapper must wait for the summary
# job (which runs after all consensus jobs) not the grid driver itself.
echo "$summary_jid" > "$OUTDIR/summary_jid.txt"
echo "Summary job ID written to: $OUTDIR/summary_jid.txt"
echo
echo "[DONE] Submitted grid."
echo "  SEACR outputs:     $OUTDIR/seacr_grid/fdr*/"
echo "  Consensus outputs: $OUTDIR/consensus_grid/fdr*/(relaxed_IgG|stringent_IgG)/common_*pct*.bed"
echo "  Collected beds:    $OUTDIR/consensus_grid_collected/"
echo "  Summary TSV:       $OUTDIR/consensus_grid_summary.tsv"
echo
echo "Tip: start with a one-threshold smoke test:"
echo "  FDR_LIST=\"0.01\" MATCH_TOKEN=\"\" $0 \"$OUTDIR\""
