#!/usr/bin/env bash
#SBATCH --job-name=overture_panel_eval
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=overture_panel_eval.%j.out
#SBATCH --error=overture_panel_eval.%j.err

# ==============================================================================
# OVERTURE Step 4 — Panel evaluation
#
# Benchmarks candidate BAR (BASIN Anchor Region) panels produced by Step 3.
# For each candidate BED file, extracts signal from a set of CUT&RUN/CUT&TAG
# CPM bigWig files and computes trimmed-median background summaries across a
# range of drop fractions. Outputs panel-evaluation diagnostics and a summary
# of passing candidate panels for comparison and final panel selection.
#
# Usage (Slurm):
#   sbatch --export=ALL,BW_DIR=/path/to/bigwigs,REGION_DIR=/path/to/candidate_beds \
#     eval_panels.sbatch
#
# Or set variables in the environment before submitting:
#   export BW_DIR=/path/to/CPM_bigwigs
#   export REGION_DIR=/path/to/candidate_region_files
#   sbatch eval_panels.sbatch
#
# Required inputs:
#   BW_DIR      Directory containing CPM bigWig files (.bw) for evaluation.
#               Must also contain bw_metadata.tsv (see templates/bw_metadata_template.tsv).
#   REGION_DIR  Directory containing candidate BAR panel BED files from Step 3.
#               These are the consensus BEDs collected into
#               consensus_grid_collected/ by the grid summary job.
#
# Optional settings (override via environment or edit below):
#   OUTDIR        Where to write output TSVs (default: BW_DIR/panel_eval_out)
#   BW_META       Path to bigWig metadata TSV (default: BW_DIR/bw_metadata.tsv)
#   REF_GENOTYPE  Reference genotype for genotype-shift calculation (default: yw)
#   MARKS         Comma-separated list of histone marks to evaluate against
#   DROP_LIST     Comma-separated list of top-fraction removal values to test
#                 (e.g., 0.10 = remove top 10% of regions by signal)
#   MIN_REGIONS   Minimum regions required per panel to be evaluated (default: 200)
#
# Outputs (in OUTDIR):
#   panel_trimmed_medians.tsv    — per-sample trimmed-median background value
#                                  for each panel × drop-fraction combination
#   replicate_stability_trimmed.tsv — CV across replicates for each panel ×
#                                  mark × genotype × drop-fraction combination
#   genotype_shift_trimmed.tsv   — delta vs. reference genotype for each panel ×
#                                  mark × drop-fraction combination
#
# Panel selection guidance:
#   Good panels have:
#     - enough regions for a stable median estimate
#     - low replicate CV (lower = more reproducible)
#     - stable behavior at moderate trimming levels
#   Genotype-dependent medBAR shifts are retained as a diagnostic only and are
#   not minimized during automated panel selection. Such shifts may reflect
#   genuine background titration associated with global changes in specific signal.
#   Panels requiring extreme trimming (>50%) to stabilize may be contaminated
#   with true signal and should be inspected carefully.
#   See OVERTURE.md for a full description of the selection criteria.
#
# bw_metadata.tsv format (tab-separated, required columns):
#   bw_path   sample_id   genotype   mark   rep   center_reads   extend_reads   bin_size   min_mapq   normalization
# Validation tracks must be uncentered, unextended, binSize=1 CPM. All rows in
# one evaluation must use the same min_mapq.
#
# Dependencies: bigWigAverageOverBed (UCSC tools), python3, pandas
# ==============================================================================

set -euo pipefail

# =========================
# Required inputs
# =========================
BW_DIR="${BW_DIR:?ERROR: BW_DIR not set. Export BW_DIR=/path/to/CPM_bigwigs}"
REGION_DIR="${REGION_DIR:?ERROR: REGION_DIR not set. Export REGION_DIR=/path/to/candidate_beds}"

# =========================
# Optional settings
# =========================
BW_META="${BW_META:-${BW_DIR}/bw_metadata.tsv}"
OUTDIR="${OUTDIR:-${BW_DIR}/panel_eval_out}"

REF_GENOTYPE="${REF_GENOTYPE:-yw}"
MARKS="${MARKS:-H3K36me2,H3K36me3,H3K27ac,panH4,RNAPIIS5P,RNAPIIS2P}"
DROP_LIST="${DROP_LIST:-0.10,0.20,0.35,0.50}"
MIN_REGIONS="${MIN_REGIONS:-200}"

mkdir -p "$OUTDIR"

# =========================
# Module loading
# =========================
if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  module load ucsctools 2>/dev/null || true
  module load python 2>/dev/null || true
fi

# =========================
# Dependency checks
# =========================
command -v bigWigAverageOverBed >/dev/null 2>&1 || {
  echo "ERROR: bigWigAverageOverBed not found on PATH." >&2
  echo "Install UCSC tools or add their bin directory to PATH." >&2
  exit 2
}
command -v python3 >/dev/null 2>&1 || {
  echo "ERROR: python3 not found on PATH." >&2
  exit 2
}

[[ -d "$BW_DIR"     ]] || { echo "ERROR: BW_DIR not found: $BW_DIR"         >&2; exit 2; }
[[ -d "$REGION_DIR" ]] || { echo "ERROR: REGION_DIR not found: $REGION_DIR" >&2; exit 2; }
[[ -s "$BW_META"    ]] || { echo "ERROR: bw_metadata.tsv not found or empty: $BW_META" >&2
                             echo "       Expected at: $BW_META" >&2
                             echo "       See templates/bw_metadata_template.tsv for format." >&2
                             exit 2; }

echo "=== OVERTURE Panel Evaluation ==="
echo "BW_DIR        : $BW_DIR"
echo "REGION_DIR    : $REGION_DIR"
echo "BW_META       : $BW_META"
echo "OUTDIR        : $OUTDIR"
echo "REF_GENOTYPE  : $REF_GENOTYPE"
echo "MARKS         : $MARKS"
echo "DROP_LIST     : $DROP_LIST"
echo "MIN_REGIONS   : $MIN_REGIONS"
echo "bigWigAvgOverBed: $(command -v bigWigAverageOverBed)"
echo

# Export settings for Python
export BW_DIR REGION_DIR BW_META OUTDIR REF_GENOTYPE MARKS DROP_LIST MIN_REGIONS

# Wait for REGION_DIR to be populated — filesystem sync on the cluster can
# lag behind job completion, so we poll briefly before starting evaluation.
echo "Checking for BED files in REGION_DIR..."
for i in $(seq 1 18); do
  count=$(ls "$REGION_DIR"/*.bed 2>/dev/null | wc -l)
  if [[ "$count" -gt 0 ]]; then
    echo "Found $count BED files — proceeding."
    break
  fi
  echo "  Waiting for BED files... attempt $i/18 (sleeping 10s)"
  sleep 10
done
count=$(ls "$REGION_DIR"/*.bed 2>/dev/null | wc -l)
[[ "$count" -gt 0 ]] || {
  echo "ERROR: No BED files appeared in $REGION_DIR after 3 minutes. Exiting." >&2
  exit 2
}

python3 - << 'PY'
import os, subprocess, statistics, hashlib, atexit
from pathlib import Path
import pandas as pd

BW_DIR      = Path(os.environ["BW_DIR"]).resolve()
REGION_DIR  = Path(os.environ["REGION_DIR"]).resolve()
BW_META     = Path(os.environ["BW_META"]).resolve()
OUTDIR      = Path(os.environ["OUTDIR"]).resolve()

REF_GENOTYPE = os.environ.get("REF_GENOTYPE", "yw")
MARKS = {m.strip() for m in os.environ.get("MARKS", "").split(",") if m.strip()}
if not MARKS:
    raise RuntimeError("MARKS is empty after parsing.")
DROP_LIST = [float(x.strip()) for x in os.environ.get("DROP_LIST", "").split(",") if x.strip()]
if not DROP_LIST:
    raise RuntimeError("DROP_LIST is empty after parsing.")
MIN_REGIONS = int(os.environ.get("MIN_REGIONS", "200"))

OUTDIR.mkdir(parents=True, exist_ok=True)
tmp_tab = OUTDIR / "tmp_aob.tab"
atexit.register(lambda: tmp_tab.unlink(missing_ok=True))

meta = pd.read_csv(BW_META, sep="\t")
required_meta = {
    "bw_path", "sample_id", "genotype", "mark", "rep",
    "center_reads", "extend_reads", "bin_size", "min_mapq", "normalization",
}
missing_meta = sorted(required_meta - set(meta.columns))
if missing_meta:
    raise RuntimeError(
        "bw_metadata.tsv is missing required coverage-contract columns: "
        + ", ".join(missing_meta)
        + ". Use the distributed template and record validation-bigWig generation parameters."
    )

def norm_bool(x):
    return str(x).strip().lower()

if not meta[~meta["center_reads"].map(norm_bool).isin({"false", "0", "no"})].empty:
    raise RuntimeError("OVERTURE panel evaluation requires center_reads=false for every validation bigWig.")
if not meta[~meta["extend_reads"].map(norm_bool).isin({"false", "0", "no"})].empty:
    raise RuntimeError("OVERTURE panel evaluation requires extend_reads=false for every validation bigWig.")
if not meta[pd.to_numeric(meta["bin_size"], errors="coerce") != 1].empty:
    raise RuntimeError("OVERTURE panel evaluation requires bin_size=1 for every validation bigWig.")
if not meta[meta["normalization"].astype(str).str.upper() != "CPM"].empty:
    raise RuntimeError("OVERTURE panel evaluation requires CPM validation bigWigs.")
mapqs = sorted(set(pd.to_numeric(meta["min_mapq"], errors="raise").astype(int).tolist()))
if len(mapqs) != 1:
    raise RuntimeError(f"Validation bigWigs must use one consistent min_mapq; found {mapqs}.")
validation_mapq = mapqs[0]

meta = meta[meta["mark"].isin(MARKS)].copy()
if meta.empty:
    raise RuntimeError(f"No bigWigs selected. Check MARKS={MARKS} against bw_metadata.tsv.")

beds_all = sorted(REGION_DIR.glob("*.bed"))
if not beds_all:
    raise RuntimeError(f"No BED files found in {REGION_DIR}")

# Step 3 emits the native BED plus explicit .chr/.nochr coordinate variants.
# The variant matching the native convention can therefore be byte-identical to
# the native BED. De-duplicate exact-content copies before evaluation so the
# same biological panel cannot appear twice in panel_selection_summary.tsv.
def sha256_path(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

seen_hashes = {}
beds = []
for bed in beds_all:
    digest = sha256_path(bed)
    if digest in seen_hashes:
        print(f"  [DEDUP] {bed.name} is identical to {seen_hashes[digest].name}; evaluating once")
        continue
    seen_hashes[digest] = bed
    beds.append(bed)

print(f"Found {len(beds_all)} BED files ({len(beds)} unique by content) and {len(meta)} bigWig samples.")
print(f"Validation contract: CPM, binSize=1, center=false, extend=false, min_mapq={validation_mapq}")

with open(OUTDIR / "coverage_contract.tsv", "w") as fh:
    fh.write("key\tvalue\n")
    fh.write("contract_version\tBASIN_CANONICAL_COVERAGE_V1\n")
    fh.write("canonical_signal\taligned_read_coverage\n")
    fh.write("canonical_bin_size\t1\n")
    fh.write("canonical_center_reads\tfalse\n")
    fh.write("canonical_extend_reads\tfalse\n")
    fh.write("normalization\tCPM\n")
    fh.write(f"min_mapq\t{validation_mapq}\n")
    fh.write("allzero_rule\tremove_BAR_rows_zero_across_all_validation_samples_before_trimming\n")
    fh.write("trim_rule\tn_trim=round(N*drop_fraction);stable_descending_signal_order\n")

def trimmed_median(vals, drop_top_frac: float) -> float:
    """HARMONY-compatible trimmed median: drop round(N*d) highest values."""
    if not vals:
        return float("nan")
    xs = sorted(vals, reverse=True)
    n_drop = int(round(len(xs) * drop_top_frac))
    kept = xs[n_drop:]
    return statistics.median(kept) if kept else float("nan")

records = []
skipped_small = 0
col_idx = 4  # mean0 from bigWigAverageOverBed; uncovered bases count as zero

def bed_region_stats(bed_path):
    sizes = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            try:
                sizes.append(int(cols[2]) - int(cols[1]))
            except ValueError:
                continue
    return (len(sizes), statistics.median(sizes)) if sizes else (0, float("nan"))

for bed in beds:
    panel = bed.name
    panel_n_regions, panel_median_bp = bed_region_stats(bed)
    sample_vectors = []

    for _, row in meta.iterrows():
        subprocess.run(
            ["bigWigAverageOverBed", str(row["bw_path"]), str(bed), str(tmp_tab)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        vals = []
        with open(tmp_tab) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= col_idx:
                    continue
                try:
                    vals.append(float(cols[col_idx]))
                except ValueError:
                    vals.append(0.0)

        if len(vals) < MIN_REGIONS:
            skipped_small += 1
            continue
        if all(v == 0.0 for v in vals):
            print(f"  [SKIP] {panel} × {row['sample_id']}: all-zero signal (likely chr/nochr mismatch)")
            continue
        sample_vectors.append((row, vals))

    if not sample_vectors:
        continue
    lengths = {len(vals) for _, vals in sample_vectors}
    if len(lengths) != 1:
        raise RuntimeError(f"Inconsistent region counts returned for panel {panel}: {sorted(lengths)}")

    n_rows = next(iter(lengths))
    keep_mask = [any(vals[i] != 0.0 for _, vals in sample_vectors) for i in range(n_rows)]
    n_after_allzero = sum(keep_mask)
    if n_after_allzero < MIN_REGIONS:
        skipped_small += len(sample_vectors)
        print(f"  [SKIP] {panel}: only {n_after_allzero} regions remain after global all-zero filtering")
        continue

    for row, vals in sample_vectors:
        vals = [v for v, keep in zip(vals, keep_mask) if keep]
        rep_val = row.get("rep", None)
        rep = None if pd.isna(rep_val) else int(rep_val)
        rec = {
            "panel": panel,
            "panel_n_regions": panel_n_regions,
            "panel_n_regions_after_allzero_filter": n_after_allzero,
            "panel_median_bp": panel_median_bp,
            "sample_id": row["sample_id"],
            "mark": row["mark"],
            "genotype": row["genotype"],
            "rep": rep,
            "n_nonzero_regions": sum(v != 0.0 for v in vals),
            "min_mapq": validation_mapq,
        }
        for d in DROP_LIST:
            rec[f"median_drop{int(round(d*100)):02d}"] = trimmed_median(vals, d)
        records.append(rec)

if not records:
    raise RuntimeError("No panel × sample combinations survived evaluation.")

df = pd.DataFrame(records)
df.to_csv(OUTDIR / "panel_trimmed_medians.tsv", sep="\t", index=False)
print("Wrote panel_trimmed_medians.tsv")
print(f"Skipped panel/sample combinations for insufficient usable regions: {skipped_small}")

drop_cols = [c for c in df.columns if c.startswith("median_drop")]
rep_rows   = []
shift_rows = []

for drop_col in drop_cols:
    # Replicate CV: population SD / mean, within genotype+mark
    for (panel, mark, genotype), sub in df.groupby(["panel", "mark", "genotype"]):
        vals = sub[drop_col].dropna().tolist()
        if len(vals) < 2:
            continue
        mean = sum(vals) / len(vals)
        sd   = statistics.pstdev(vals)
        cv   = (sd / mean) if mean != 0 else float("nan")
        rep_rows.append({
            "drop":     drop_col,
            "panel":    panel,
            "mark":     mark,
            "genotype": genotype,
            "n":        len(vals),
            "mean":     mean,
            "cv":       cv
        })

    # Genotype shift diagnostic: genotype median – reference-genotype median,
    # per panel + mark. Reported for interpretation; not used for filtering or ranking.
    for (panel, mark), sub in df.groupby(["panel", "mark"]):
        ref = sub[sub["genotype"] == REF_GENOTYPE][drop_col].median()
        if pd.isna(ref):
            continue
        for genotype, sub2 in sub.groupby("genotype"):
            if genotype == REF_GENOTYPE:
                continue
            med = sub2[drop_col].median()
            shift_rows.append({
                "drop":                   drop_col,
                "panel":                  panel,
                "mark":                   mark,
                "genotype":               genotype,
                f"delta_vs_{REF_GENOTYPE}": med - ref
            })

pd.DataFrame(rep_rows).to_csv(
    OUTDIR / "replicate_stability_trimmed.tsv", sep="\t", index=False)
pd.DataFrame(shift_rows).to_csv(
    OUTDIR / "genotype_shift_trimmed.tsv", sep="\t", index=False)

print("Wrote replicate_stability_trimmed.tsv")
print("Wrote genotype_shift_trimmed.tsv")
print("Drop fractions evaluated:", sorted({int(c.replace("median_drop","")) for c in drop_cols}))

# ==============================================================================
# Panel selection summary
# ==============================================================================
print()
print("=== Panel selection summary ===")

# Thresholds
MIN_REGIONS_FILTER = int(os.environ.get("EVAL_MIN_REGIONS", "200"))
MAX_REGIONS_FILTER = int(os.environ.get("EVAL_MAX_REGIONS", "5000"))
MAX_CV_FILTER      = float(os.environ.get("EVAL_MAX_CV", "0.12"))
EVAL_DROP          = os.environ.get("EVAL_DROP_FRAC", "median_drop20")

# --- Panel size from medians TSV ---
panel_size = df.groupby("panel")[["panel_n_regions", "panel_median_bp"]].first().reset_index()

# --- Mean CV across marks and genotypes at the requested drop ---
cv_df = pd.DataFrame(rep_rows)
if cv_df.empty or EVAL_DROP not in set(cv_df["drop"]):
    raise RuntimeError(
        f"EVAL_DROP_FRAC={EVAL_DROP!r} was not produced. Available drops: "
        f"{sorted(set(cv_df['drop'])) if not cv_df.empty else []}"
    )
cv_at_drop = cv_df[cv_df["drop"] == EVAL_DROP].groupby("panel")["cv"].mean().reset_index()
cv_at_drop.columns = ["panel", "mean_cv"]

# --- Max absolute genotype shift across marks at the requested drop ---
# Diagnostic only: retained in the summary but not used for filtering or ranking.
shift_df = pd.DataFrame(shift_rows)
shift_col = f"delta_vs_{REF_GENOTYPE}"
if shift_df.empty or EVAL_DROP not in set(shift_df["drop"]):
    raise RuntimeError(
        f"EVAL_DROP_FRAC={EVAL_DROP!r} was not produced for genotype shifts. Available drops: "
        f"{sorted(set(shift_df['drop'])) if not shift_df.empty else []}"
    )
shift_at_drop = shift_df[shift_df["drop"] == EVAL_DROP].copy()

shift_at_drop["abs_shift"] = shift_at_drop[shift_col].abs()
max_shift = shift_at_drop.groupby("panel")["abs_shift"].max().reset_index()
max_shift.columns = ["panel", "max_abs_genotype_shift"]

# --- Merge all metrics ---
summary = panel_size.merge(cv_at_drop, on="panel", how="left")
summary = summary.merge(max_shift, on="panel", how="left")

# --- Apply filters ---
passing = summary[
    (summary["panel_n_regions"] >= MIN_REGIONS_FILTER) &
    (summary["panel_n_regions"] <= MAX_REGIONS_FILTER) &
    (summary["mean_cv"] <= MAX_CV_FILTER) 
].copy()

n_total   = len(summary)
n_passing = len(passing)
print(f"Panels evaluated : {n_total}")
print(f"Panels passing filters : {n_passing}")
print(f"  n_regions : {MIN_REGIONS_FILTER}–{MAX_REGIONS_FILTER}")
print(f"  mean CV   : ≤ {MAX_CV_FILTER} (at {EVAL_DROP})")

if passing.empty:
    print()
    print("WARNING: No panels passed all filters.")
    print("Consider relaxing thresholds or reviewing the raw TSVs.")
    summary.sort_values("mean_cv").to_csv(
        OUTDIR / "panel_selection_summary.tsv", sep="\t", index=False)
    print("Wrote panel_selection_summary.tsv (all panels, unfiltered)")
else:
    # Sort passing panels by replicate CV for convenient review.
    # Panel size is used as an eligibility criterion, not an optimization target.
    # Genotype shift is retained as a diagnostic and is not used for ranking.
    passing = passing.sort_values(
        ["mean_cv", "panel_n_regions"],
        ascending=[True, False]
    ).reset_index(drop=True)

    out_cols = [
        "panel",
        "panel_n_regions",
        "panel_median_bp",
        "mean_cv",
        "max_abs_genotype_shift"
    ]

    passing[out_cols].to_csv(
        OUTDIR / "panel_selection_summary.tsv", sep="\t", index=False
    )

    print()
    print("Passing candidate panels (sorted by mean replicate CV):")
    print(passing[out_cols].head(10).to_string(index=False))
    print()
    print("Wrote panel_selection_summary.tsv")
    print("No single panel is automatically designated as the recommended panel.")
    print("Genotype shift is reported for diagnostic interpretation only.")
    print("See OVERTURE README — Panel Selection for interpretation guidance.")
PY

echo
echo "[DONE] Panel evaluation complete."
echo "Outputs in: $OUTDIR"
echo
echo "Next step: review panel_selection_summary.tsv to select your final BAR panel."
echo "See OVERTURE.md — Panel Selection section for guidance on interpreting"
echo "replicate CV, genotype shift, and panel size."
ls -lh "$OUTDIR"
