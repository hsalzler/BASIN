#!/usr/bin/env bash
#
# SEACR_Consensus_IgG_Peak_Selection.sbatch
#
# Builds 80/90/100% consensus peaks (BED5) from SEACR peak BED files.
#
# Multi-mode (recommended):
#   sbatch --export=ALL,SEACR_DIR=/path/to/run/seacr,CONS_OUTDIR=/path/to/run/consensus SEACR_Consensus_IgG_Peak_Selection.sbatch
#   sbatch --export=ALL,RELAXED_DIR=/path/to/.../relaxed/IgG,STRINGENT_DIR=/path/to/.../stringent/IgG,CONS_OUTDIR=/path/to/run/consensus SEACR_Consensus_IgG_Peak_Selection.sbatch
#
# Single-mode (backwards compatible):
#   sbatch --export=ALL,WORKDIR=/path/to/peaks,OUTDIR=/path/to/out SEACR_Consensus_IgG_Peak_Selection.sbatch
#
# Outputs (BED5):
#   common_80pct.bed, common_90pct.bed, common_100pct.bed
#   columns: chrom start end peak_00001 support_count
#
#SBATCH --job-name=overture_consensus
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=consensus_peaks.%j.out
#SBATCH --error=consensus_peaks.%j.err

set -euo pipefail
export LC_ALL=C

# ------------------------------------------------------------------------------
# Optional: Length filtering of consensus outputs (post-process)
# ------------------------------------------------------------------------------
MAKE_LEN_FILTERED="${MAKE_LEN_FILTERED:-yes}"        # yes/no
MIN_LEN="${MIN_LEN:-50}"                             # integer bp
LEN_FILTER_SUFFIX="${LEN_FILTER_SUFFIX:-_filter_away_lt${MIN_LEN}bp}"

# Prefer gawk if available (more consistent); fall back to awk.
AWK="awk"
command -v gawk >/dev/null 2>&1 && AWK="gawk"

# ------------------------------------------------------------------------------
# Defaults (original single-run mode)
# ------------------------------------------------------------------------------
WORKDIR="${WORKDIR:-$PWD}"
OUTDIR="${OUTDIR:-$WORKDIR}"
PATTERN="${PATTERN:-*.bed}"
MODE="${MODE:-merged}"                  # minimal | merged
MERGE_DISTANCE="${MERGE_DISTANCE:-0}"   # used only when MODE=merged
ID_PREFIX="${ID_PREFIX:-peak_}"
ID_WIDTH="${ID_WIDTH:-5}"

# ------------------------------------------------------------------------------
# Multi-run inputs (optional)
# ------------------------------------------------------------------------------
SEACR_DIR="${SEACR_DIR:-}"
RELAXED_DIR="${RELAXED_DIR:-}"
STRINGENT_DIR="${STRINGENT_DIR:-}"
CONS_OUTDIR="${CONS_OUTDIR:-}"

if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  module load bedtools 2>/dev/null || true
  module load python 2>/dev/null || true
fi

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found on PATH." >&2; exit 1; }
command -v python3  >/dev/null 2>&1 || { echo "ERROR: python3 not found on PATH."  >&2; exit 1; }

# Derive relaxed/stringent dirs from SEACR_DIR if not explicitly set
if [[ -n "$SEACR_DIR" ]]; then
  [[ -n "$RELAXED_DIR"   ]] || RELAXED_DIR="${SEACR_DIR%/}/SEACR_peaks_relaxed_CPM/IgG"
  [[ -n "$STRINGENT_DIR" ]] || STRINGENT_DIR="${SEACR_DIR%/}/SEACR_peaks_stringent_CPM/IgG"
fi

# Determine mode:
MULTI_MODE=0
if [[ -n "$RELAXED_DIR" || -n "$STRINGENT_DIR" ]]; then
  MULTI_MODE=1
fi

# Choose output root
if [[ "$MULTI_MODE" -eq 1 ]]; then
  if [[ -z "$CONS_OUTDIR" ]]; then
    if [[ -n "$SEACR_DIR" ]]; then
      CONS_OUTDIR="$(cd "${SEACR_DIR%/}/.." 2>/dev/null && pwd)/consensus"
    else
      CONS_OUTDIR="$PWD/consensus"
    fi
  fi
  mkdir -p "$CONS_OUTDIR"
else
  mkdir -p "$OUTDIR"
fi

echo "Host=$(hostname) JobID=${SLURM_JOB_ID:-N/A}"
echo "PATTERN=$PATTERN"
echo "MODE=$MODE"
echo "MERGE_DISTANCE=$MERGE_DISTANCE"
echo "ID_PREFIX=$ID_PREFIX"
echo "ID_WIDTH=$ID_WIDTH"
echo "MAKE_LEN_FILTERED=$MAKE_LEN_FILTERED"
echo "MIN_LEN=$MIN_LEN"
echo "LEN_FILTER_SUFFIX=$LEN_FILTER_SUFFIX"
echo "MULTI_MODE=$MULTI_MODE"
if [[ "$MULTI_MODE" -eq 1 ]]; then
  echo "SEACR_DIR=${SEACR_DIR:-<unset>}"
  echo "RELAXED_DIR=${RELAXED_DIR:-<unset>}"
  echo "STRINGENT_DIR=${STRINGENT_DIR:-<unset>}"
  echo "CONS_OUTDIR=$CONS_OUTDIR"
else
  echo "WORKDIR=$WORKDIR"
  echo "OUTDIR=$OUTDIR"
fi
echo

TMPDIR="${SLURM_TMPDIR:-/tmp}"
PY="$TMPDIR/consensus_peaks.py"

cat > "$PY" << 'PYCODE'
#!/usr/bin/env python3
"""
Build 80/90/100% consensus peaks from BEDs in a directory.

Outputs BED5:
  chrom start end peak_00001 support_count

Support definition:
  A base is "supported" by a file if that file has any peak covering that base (>=1 bp overlap).

Implementation:
  - bedtools sort each BED
  - bedtools multiinter across all sorted BEDs to get segments annotated with nFiles
  - filter segments to min_support for each percentage
  - optionally bedtools merge (keeping max support)
  - assign peak IDs
  - write inputs.used.txt for reproducibility
"""

from __future__ import annotations

import argparse
import glob
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple

PCTS: Tuple[float, float, float] = (0.80, 0.90, 1.00)


def require_exe(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def run_capture(cmd: List[str]) -> str:
    try:
        proc = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return proc.stdout
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "Command failed:\n"
            f"  {' '.join(cmd)}\n\n"
            f"STDOUT:\n{e.stdout}\n\n"
            f"STDERR:\n{e.stderr}\n"
        ) from e


def list_beds(workdir: Path, pattern: str) -> List[Path]:
    hits = sorted(Path(p).resolve() for p in glob.glob(str(workdir / pattern)))
    hits = [p for p in hits if p.is_file() and p.suffix.lower() == ".bed"]
    return hits


def bedtools_sort(in_bed: Path, out_bed: Path) -> None:
    out_bed.write_text(run_capture(["bedtools", "sort", "-i", str(in_bed)]))


def bedtools_multiinter(sorted_beds: List[Path], out_path: Path) -> None:
    out_path.write_text(run_capture(["bedtools", "multiinter", "-i", *map(str, sorted_beds)]))


def filter_multiinter_to_bed4(multi_path: Path, min_support: int, out_path: Path) -> None:
    # multiinter columns: chrom start end nfiles list
    # emit BED4: chrom start end support
    with multi_path.open("r") as r, out_path.open("w") as w:
        for line in r:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, nfiles = parts[0], parts[1], parts[2], parts[3]
            try:
                support = int(nfiles)
            except ValueError:
                continue
            if support >= min_support:
                w.write(f"{chrom}\t{start}\t{end}\t{support}\n")


def bedtools_merge_bed4(in_bed4: Path, out_bed4: Path, merge_distance: int) -> None:
    # Merge segments within distance; keep max support across merged pieces.
    out_bed4.write_text(
        run_capture(
            [
                "bedtools",
                "merge",
                "-i",
                str(in_bed4),
                "-d",
                str(merge_distance),
                "-c",
                "4",
                "-o",
                "max",
            ]
        )
    )


def assign_ids_bed4_to_bed5(in_bed4: Path, out_bed5: Path, prefix: str, width: int) -> None:
    i = 0
    with in_bed4.open("r") as r, out_bed5.open("w") as w:
        for line in r:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            i += 1
            peak_id = f"{prefix}{i:0{width}d}"
            w.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{peak_id}\t{parts[3]}\n")


def threshold(pct: float, n: int) -> int:
    # ceil(pct*n), with small epsilon to avoid float edge cases
    return int(math.ceil((pct * n) - 1e-12))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build 80/90/100% consensus peaks (BED5 with peak IDs).")
    ap.add_argument("--dir", default=".", help="Directory containing BEDs (default: .)")
    ap.add_argument("--pattern", default="*.bed", help="Glob pattern (default: *.bed)")
    ap.add_argument("--outdir", default=None, help="Output directory (default: same as --dir)")
    ap.add_argument(
        "--mode",
        choices=("minimal", "merged"),
        default="merged",
        help="minimal keeps segments; merged merges (default: merged)",
    )
    ap.add_argument("--merge-distance", type=int, default=0, help="For merged mode: merge within distance (default: 0)")
    ap.add_argument("--id-prefix", default="peak_", help="Peak ID prefix (default: peak_)")
    ap.add_argument("--id-width", type=int, default=5, help="Zero-pad width (default: 5)")
    return ap.parse_args()


def main() -> int:
    require_exe("bedtools")

    args = parse_args()
    workdir = Path(args.dir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else workdir
    outdir.mkdir(parents=True, exist_ok=True)

    beds = list_beds(workdir, args.pattern)
    if not beds:
        raise RuntimeError(f"No BED files found in {workdir} matching '{args.pattern}'")

    # audit trail
    (outdir / "inputs.used.txt").write_text("\n".join(str(p) for p in beds) + "\n")

    n = len(beds)

    with tempfile.TemporaryDirectory(prefix="consensus_peaks_") as td:
        tmp = Path(td)

        sorted_beds: List[Path] = []
        for i, bed in enumerate(beds, start=1):
            s = tmp / f"sorted_{i}.bed"
            bedtools_sort(bed, s)
            sorted_beds.append(s)

        multi = tmp / "multiinter.bed"
        bedtools_multiinter(sorted_beds, multi)

        for pct in PCTS:
            thr = threshold(pct, n)
            pct_label = f"{int(round(pct * 100)):d}pct"

            filtered_bed4 = tmp / f"filtered_{pct_label}.bed4"
            filter_multiinter_to_bed4(multi, thr, filtered_bed4)

            bed4_for_ids = filtered_bed4
            if args.mode == "merged":
                merged_bed4 = tmp / f"merged_{pct_label}.bed4"
                bedtools_merge_bed4(filtered_bed4, merged_bed4, args.merge_distance)
                bed4_for_ids = merged_bed4

            out_path = outdir / f"common_{pct_label}.bed"
            assign_ids_bed4_to_bed5(bed4_for_ids, out_path, args.id_prefix, args.id_width)

            print(f"Wrote {out_path} (threshold {thr}/{n}, mode={args.mode})")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(1)
PYCODE

chmod +x "$PY"

filter_beds_in_dir() {
  local d="$1"
  local min_len="$2"
  local suffix="$3"

  shopt -s nullglob
  local beds=( "$d"/common_*pct.bed )
  shopt -u nullglob

  if [[ "${#beds[@]}" -eq 0 ]]; then
    echo "[WARN] No consensus beds found to filter in: $d" >&2
    return 0
  fi

  for bed in "${beds[@]}"; do
    [[ "$bed" == *"${suffix}.bed" ]] && continue
    local out="${bed%.bed}${suffix}.bed"
    if [[ -s "$out" ]]; then
      echo "[SKIP] filtered exists: $(basename "$out")"
      continue
    fi

    "$AWK" -v OFS="\t" -v m="$min_len" 'BEGIN{FS=OFS="\t"} ($3-$2) >= m {print $0}' \
      "$bed" > "$out"

    echo "[OK  ] wrote $(basename "$out")"
  done
}

summarize_beds_in_dir() {
  local d="$1"
  local pattern_glob="$2"
  local out_tsv="$3"

  shopt -s nullglob
  local beds=( "$d"/$pattern_glob )
  shopt -u nullglob

  {
    printf "file\tn\tbp1\tlt10\tlt50\tmean_len\tfrac_bp1\n"
    for bed in "${beds[@]}"; do
      [[ -e "$bed" ]] || continue
      "$AWK" -v f="$(basename "$bed")" 'BEGIN{FS=OFS="\t"}
        {len=$3-$2; n++; if(len==1)n1++; if(len<10)n10++; if(len<50)n50++; sum+=len}
        END{
          if(n==0) printf "%s\t0\t0\t0\t0\tNA\tNA\n", f
          else     printf "%s\t%d\t%d\t%d\t%d\t%.2f\t%.6f\n", f,n,n1,n10,n50,(sum/n),(n1/n)
        }' "$bed"
    done
  } > "$out_tsv"

  echo "[OK  ] wrote $out_tsv"
}

run_one() {
  local _workdir="$1"
  local _outdir="$2"
  local _label="$3"

  if [[ ! -d "$_workdir" ]]; then
    echo "[WARN] ${_label}: WORKDIR not found, skipping: $_workdir" >&2
    return 0
  fi

  mkdir -p "$_outdir"

  echo "---- CONSENSUS: ${_label} ----"
  echo "WORKDIR=$_workdir"
  echo "OUTDIR=$_outdir"
  echo

  # ---- chr-prefix normalization (temp only; originals untouched) ----
  local td normdir
  td="$(mktemp -d "${TMPDIR%/}/consensus_norm_${_label}.XXXXXX")"
  normdir="$td/norm"
  mkdir -p "$normdir"

  mapfile -t beds < <(find "$_workdir" -maxdepth 1 -type f -name "$PATTERN" | sort)
  if [[ "${#beds[@]}" -eq 0 ]]; then
    echo "[WARN] ${_label}: No BED files matching '$PATTERN' in $_workdir" >&2
    rm -rf "$td"
    return 0
  fi

  local chr_yes=0 chr_no=0 first
  for f in "${beds[@]}"; do
    first="$("$AWK" 'BEGIN{FS="\t"} $0!~/^#/ && NF>=3 {print $1; exit}' "$f" || true)"
    if [[ "$first" == chr* ]]; then
      chr_yes=$((chr_yes+1))
    else
      chr_no=$((chr_no+1))
    fi
  done

  local target="nochr"
  if [[ "$chr_yes" -gt "$chr_no" ]]; then
    target="chr"
  fi
  echo "[INFO] ${_label}: chr-prefix files=$chr_yes, nochr files=$chr_no -> normalizing to '$target'"

  for f in "${beds[@]}"; do
    local base out
    base="$(basename "$f")"
    out="$normdir/$base"
    if [[ "$target" == "nochr" ]]; then
      sed -E 's/^chr//' "$f" > "$out"
    else
      "$AWK" 'BEGIN{OFS="\t"} /^#/ {print; next} NF<3 {next} { if($1 !~ /^chr/) $1="chr"$1; print }' \
        "$f" > "$out"
    fi
  done

  python3 "$PY" \
    --dir "$normdir" \
    --pattern "$PATTERN" \
    --outdir "$_outdir" \
    --mode "$MODE" \
    --merge-distance "$MERGE_DISTANCE" \
    --id-prefix "$ID_PREFIX" \
    --id-width "$ID_WIDTH"

  if [[ "$MAKE_LEN_FILTERED" == "yes" ]]; then
    echo "[INFO] ${_label}: creating length-filtered consensus copies (len >= ${MIN_LEN} bp)"
    filter_beds_in_dir "$_outdir" "$MIN_LEN" "$LEN_FILTER_SUFFIX"
  fi

  summarize_beds_in_dir "$_outdir" "common_*pct.bed" "${_outdir}/consensus_len_summary.tsv"

  if [[ "$MAKE_LEN_FILTERED" == "yes" ]]; then
    summarize_beds_in_dir "$_outdir" "common_*pct${LEN_FILTER_SUFFIX}.bed" \
      "${_outdir}/consensus_len_summary${LEN_FILTER_SUFFIX}.tsv"
  fi

  rm -rf "$td"
}

if [[ "$MULTI_MODE" -eq 1 ]]; then
  if [[ -n "${RELAXED_DIR:-}" ]]; then
    run_one "$RELAXED_DIR" "$CONS_OUTDIR/relaxed_IgG" "RELAXED"
  fi
  if [[ -n "${STRINGENT_DIR:-}" ]]; then
    run_one "$STRINGENT_DIR" "$CONS_OUTDIR/stringent_IgG" "STRINGENT"
  fi

  echo
  echo "[DONE] Multi-mode consensus complete."
  echo "Outputs:"
  [[ -n "${RELAXED_DIR:-}"   ]] && echo "  $CONS_OUTDIR/relaxed_IgG/common_{80,90,100}pct.bed"
  [[ -n "${STRINGENT_DIR:-}" ]] && echo "  $CONS_OUTDIR/stringent_IgG/common_{80,90,100}pct.bed"
else
  run_one "$WORKDIR" "$OUTDIR" "SINGLE"

  echo
  echo "[DONE] Single-mode consensus complete."
  echo "Outputs: $OUTDIR/common_{80,90,100}pct.bed"
fi
