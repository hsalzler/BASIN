#!/usr/bin/env bash
#
#SBATCH -J overture_seacr
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH -o seacr_IgG_only.%j.out
#SBATCH -e seacr_IgG_only.%j.err

set -euo pipefail

# --- SEACR availability ---
if command -v module >/dev/null 2>&1; then
  module load SEACR/1.3 2>/dev/null || module load seacr/1.3 2>/dev/null || true
fi
SEACR="${SEACR:-SEACR_1.3.sh}"
if ! command -v "$SEACR" >/dev/null 2>&1; then
  echo "[ERROR] SEACR not found. Pass SEACR=/full/path/SEACR_1.3.sh" >&2
  exit 1
fi

# === CONFIG ===
SEARCH_DIR="${SEARCH_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
MATCH_TOKEN="${MATCH_TOKEN-}"
THREADS="${THREADS:-${SLURM_NTASKS:-4}}"
SORT_INPUT="${SORT_INPUT:-no}"             # yes/no
TMPDIR_BASE="${TMPDIR_BASE:-/tmp}"
CASE_INSENSITIVE="${CASE_INSENSITIVE:-1}"  # 1=case-insensitive
RECURSE="${RECURSE:-0}"                    # 0=top-level only

# --- Optional: chromosome filtering (recommended for Drosophila) ---
# Many public BAMs/bedGraphs contain random/unplaced contigs. If you are using
# IgG peaks as "universal" normalization anchors, those contigs can add noise
# and also create downstream mismatches.
#
# Default behavior here is to KEEP ONLY the standard Drosophila chromosomes,
# accepting either naming convention (with or without "chr" prefix).
#
# Set STANDARD_CHR_ONLY=no to disable filtering.
# Set CHR_FILTER_MODE to one of: both | chr | nochr
#   both  : keep either chrX or X, etc. (default; most robust)
#   chr   : keep only chrX/chr2L/... style
#   nochr : keep only X/2L/... style
STANDARD_CHR_ONLY="${STANDARD_CHR_ONLY:-yes}"   # yes/no
CHR_FILTER_MODE="${CHR_FILTER_MODE:-both}"       # both|chr|nochr
DROS_CHRS_CHR="${DROS_CHRS_CHR:-chrX,chrY,chr2L,chr2R,chr3L,chr3R,chr4,chrM}"
DROS_CHRS_NOCHR="${DROS_CHRS_NOCHR:-X,Y,2L,2R,3L,3R,4,M}"

# Output dirs: use your existing top-level dirs, we add IgG/ inside them
BASE_RELAXED_DIR="${BASE_RELAXED_DIR:-SEACR_peaks_relaxed_CPM}"
BASE_STRINGENT_DIR="${BASE_STRINGENT_DIR:-SEACR_peaks_stringent_CPM}"
BASE_LOG_DIR="${BASE_LOG_DIR:-seacr_logs}"
RELAXED_DIR="${BASE_RELAXED_DIR%/}/IgG"
STRINGENT_DIR="${BASE_STRINGENT_DIR%/}/IgG"
LOG_DIR="${BASE_LOG_DIR%/}/IgG"

# No-control cutoff (required when no IgG control track is used)
FDR_FRAC="${FDR_FRAC:-0.01}"  # e.g., 0.05 for more peaks

# Tag to append to output names, e.g., 0.01 -> frac0p01
FRAC_TAG="${FRAC_TAG:-$(awk -v x="$FDR_FRAC" 'BEGIN{gsub(/\./,"p",x); print "frac"x}')}"

# optional sanity check
if ! awk -v x="$FDR_FRAC" 'BEGIN{exit !(x>0 && x<=1)}'; then
  echo "[ERROR] FDR_FRAC must be in (0,1], got: $FDR_FRAC" >&2
  exit 1
fi

# === sanity ===
if [[ "$SORT_INPUT" == "yes" ]] && ! command -v bedtools >/dev/null 2>&1; then
  echo "[ERROR] SORT_INPUT=yes but bedtools not in PATH." >&2
  exit 1
fi

mkdir -p "$RELAXED_DIR" "$STRINGENT_DIR" "$LOG_DIR"

# Parallel
if command -v parallel >/dev/null 2>&1; then
  PARALLEL="parallel -0 -j ${THREADS} --halt now,fail=1 --xargs -n 1"
else
  PARALLEL="xargs -0 -I{} -P ${THREADS} bash -lc"
fi

echo "[INFO] Using SEARCH_DIR=${SEARCH_DIR}"
echo "[INFO] MATCH_TOKEN=${MATCH_TOKEN}  RECURSE=${RECURSE}  CASE_INSENSITIVE=${CASE_INSENSITIVE}"
echo "[INFO] THREADS=${THREADS}  SORT_INPUT=${SORT_INPUT}  FDR_FRAC=${FDR_FRAC} (${FRAC_TAG})"
echo "[INFO] STANDARD_CHR_ONLY=${STANDARD_CHR_ONLY}  CHR_FILTER_MODE=${CHR_FILTER_MODE}"
echo "[INFO] Outputs: relaxed=${RELAXED_DIR} stringent=${STRINGENT_DIR} logs=${LOG_DIR}"

run_seacr_single() {
  local exp_bg="$1"
  local base tmpd exp_in rel_prefix str_prefix log_r log_s

  [[ -s "$exp_bg" ]] || { echo "[WARN] Missing/empty: $exp_bg" >&2; return 0; }

  base="$(basename "$exp_bg")"
  base="${base%.bg.gz}"; base="${base%.bg}"
  base="${base%.bedgraph.gz}"; base="${base%.bedgraph}"

  rel_prefix="${RELAXED_DIR}/${base}.SEACR_relaxed_${FRAC_TAG}"
  str_prefix="${STRINGENT_DIR}/${base}.SEACR_stringent_${FRAC_TAG}"
  log_r="${LOG_DIR}/${base}.relaxed_${FRAC_TAG}.log"
  log_s="${LOG_DIR}/${base}.stringent_${FRAC_TAG}.log"

  if ls "${rel_prefix}"* >/dev/null 2>&1 && ls "${str_prefix}"* >/dev/null 2>&1; then
    echo "[SKIP] $base (outputs exist)"
    return 0
  fi

  tmpd="$(mktemp -d "${TMPDIR_BASE}/seacr_${base}.XXXXXX")"

  if [[ "$exp_bg" == *.gz ]]; then
    zcat "$exp_bg" > "${tmpd}/exp.bg"
  else
    cp -f "$exp_bg" "${tmpd}/exp.bg"
  fi

  # Optional: filter to standard Drosophila chromosomes only
  if [[ "$STANDARD_CHR_ONLY" == "yes" ]]; then
    # Build an awk keep-set in a way that's safe for commas/spaces
    keep_chr="$(echo "$DROS_CHRS_CHR" | tr ',' ' ' | xargs)"
    keep_nochr="$(echo "$DROS_CHRS_NOCHR" | tr ',' ' ' | xargs)"

    awk -v mode="$CHR_FILTER_MODE" -v keep_chr="$keep_chr" -v keep_nochr="$keep_nochr" '
      BEGIN{
        OFS="\t";
        n=split(keep_chr, a, /[[:space:]]+/);
        for(i=1;i<=n;i++) if(a[i]!="") K[a[i]]=1;
        n=split(keep_nochr, b, /[[:space:]]+/);
        for(i=1;i<=n;i++) if(b[i]!="") N[b[i]]=1;
      }
      {
        chr=$1;
        if(mode=="chr"){
          if(chr in K) print;
        } else if(mode=="nochr"){
          if(chr in N) print;
        } else {
          if((chr in K) || (chr in N)) print;
        }
      }
    ' "${tmpd}/exp.bg" > "${tmpd}/exp.stdchr.bg"

    if [[ ! -s "${tmpd}/exp.stdchr.bg" ]]; then
      echo "[WARN] After STANDARD_CHR_ONLY filtering, no lines remained for: $base" >&2
      echo "[HINT] Disable with STANDARD_CHR_ONLY=no or adjust DROS_CHRS_* lists." >&2
      rm -rf "$tmpd"
      return 0
    fi
    mv "${tmpd}/exp.stdchr.bg" "${tmpd}/exp.bg"
  fi

  if [[ "$SORT_INPUT" == "yes" ]]; then
    bedtools sort -i "${tmpd}/exp.bg" > "${tmpd}/exp.sorted.bg"
    mv "${tmpd}/exp.sorted.bg" "${tmpd}/exp.bg"
  fi
  exp_in="${tmpd}/exp.bg"

  # --- no-control: use numeric fraction + 'non' ---
  if ! ls "${rel_prefix}"* >/dev/null 2>&1; then
    echo "[RUN ] RELAXED (non, frac=$FDR_FRAC): $base"
    if ! "$SEACR" "$exp_in" "$FDR_FRAC" non relaxed "$rel_prefix" >"$log_r" 2>&1; then
      echo "[FAIL] RELAXED $base (see $log_r)" >&2
    fi
  else
    echo "[SKIP] RELAXED exists for $base"
  fi

  if ! ls "${str_prefix}"* >/dev/null 2>&1; then
    echo "[RUN ] STRINGENT (non, frac=$FDR_FRAC): $base"
    if ! "$SEACR" "$exp_in" "$FDR_FRAC" non stringent "$str_prefix" >"$log_s" 2>&1; then
      echo "[FAIL] STRINGENT $base (see $log_s)" >&2
    fi
  else
    echo "[SKIP] STRINGENT exists for $base"
  fi

  rm -rf "$tmpd"
}
export -f run_seacr_single
export RELAXED_DIR STRINGENT_DIR LOG_DIR SEACR SORT_INPUT TMPDIR_BASE FDR_FRAC FRAC_TAG

# Collect *IgG* bedGraphs from SEARCH_DIR
find_opts=(-type f)
[[ "$RECURSE" == "0" ]] && find_opts=(-maxdepth 1 "${find_opts[@]}")

name_flag="-name"
[[ "$CASE_INSENSITIVE" == "1" ]] && name_flag="-iname"

mapfile -d '' BG_LIST < <(
  find "${SEARCH_DIR}" "${find_opts[@]}" \
    \( \
      $name_flag "*${MATCH_TOKEN}*.bg" -o \
      $name_flag "*${MATCH_TOKEN}*.bg.gz" -o \
      $name_flag "*${MATCH_TOKEN}*.bedgraph" -o \
      $name_flag "*${MATCH_TOKEN}*.bedgraph.gz" \
    \) -print0 2>/dev/null || true
)

# Fallback globbing (top-level only)
if [[ ${#BG_LIST[@]} -eq 0 && "$RECURSE" == "0" ]]; then
  shopt -s nullglob
  [[ "$CASE_INSENSITIVE" == "1" ]] && shopt -s nocaseglob
  for f in \
    "${SEARCH_DIR}"/*"${MATCH_TOKEN}"*.bg \
    "${SEARCH_DIR}"/*"${MATCH_TOKEN}"*.bg.gz \
    "${SEARCH_DIR}"/*"${MATCH_TOKEN}"*.bedgraph \
    "${SEARCH_DIR}"/*"${MATCH_TOKEN}"*.bedgraph.gz; do
    [[ -f "$f" ]] && BG_LIST+=("$f")
  done
  [[ "$CASE_INSENSITIVE" == "1" ]] && shopt -u nocaseglob
  shopt -u nullglob
fi

if [[ ${#BG_LIST[@]} -eq 0 ]]; then
  echo "[WARN] No bedGraphs found in ${SEARCH_DIR} matching *${MATCH_TOKEN}*.{bg,bg.gz,bedgraph,bedgraph.gz}"
  echo "[HINT] ls -1 ${SEARCH_DIR}/*${MATCH_TOKEN}*.bg* ${SEARCH_DIR}/*${MATCH_TOKEN}*.bedgraph*"
  exit 0
fi

printf "[INFO] Found %d bedGraph(s):\n" "${#BG_LIST[@]}"
for p in "${BG_LIST[@]}"; do echo "  - $(basename "$p")"; done

printf '%s\0' "${BG_LIST[@]}" | \
  { if [[ "$PARALLEL" == parallel* ]]; then eval "$PARALLEL run_seacr_single"; else eval "$PARALLEL 'run_seacr_single {}'"; fi; }

echo "[DONE] SEACR (non) runs complete."
echo "[INFO] Relaxed dir:   $RELAXED_DIR"
echo "[INFO] Stringent dir: $STRINGENT_DIR"
echo "[INFO] Logs dir:      $LOG_DIR"
