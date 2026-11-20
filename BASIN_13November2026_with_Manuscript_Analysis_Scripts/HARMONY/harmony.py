#!/usr/bin/env python3
"""
HARMONY -- BASIN normalization pipeline (Pipeline 2).

Coverage contract
-----------------
BAR-based scaling is defined on one canonical coverage representation:

  * aligned read coverage (not fragment-extended)
  * no read centering
  * 1-bp bins
  * the run's configured minimum MAPQ

The canonical representation is intentionally NOT user-configurable for
centering, extension, or bin size. This keeps BAR quantification in HARMONY
compatible with the coverage geometry used to define BARs in OVERTURE and with
CODA's pre-normalization QC. MAPQ remains configurable and is recorded as run
metadata; it must be applied consistently to the canonical NONE and CPM tracks
within a run.

Scaling method (the core of BASIN)
----------------------------------
HARMONY scales every sample so that its signal is expressed relative to its own
empirical IgG background, measured at the BASIN Anchor Regions (BARs):

    scale_factor = 1 / median(raw canonical BAR coverage after trimming)

The scale factor is estimated from RAW (normalizeUsing=None) canonical coverage
and applied to a canonical scaled bigWig made from the same read representation.
This preserves the exact coordinate-wise interpretation:

    canonical_scaled(x) = canonical_raw(x) / median(canonical_raw_BAR)

HARMONY also writes a canonical CPM bigWig using the same geometry. That CPM
track is an ANALYTICAL QC input for CODA, not a freely configurable display
track.

Optional downstream bigWigs
---------------------------
After the BASIN scale factor has been estimated, HARMONY can additionally make
separate downstream CPM and/or BASIN-scaled bigWigs. These derivative tracks
can be centered without altering the canonical files used for normalization or
CODA. Read extension is permitted only on these explicitly labeled derivative
tracks and is warned about because extension makes integrated coverage depend
on inferred fragment length; such tracks should not be treated as directly
quantitative across libraries unless fragment-length distributions have been
shown to be comparable.

This runner does not use bigwigCompare. An earlier approach multiplied a CPM
bigWig by a CPM-derived scale factor via bigwigCompare; that double-corrected
for library size and silently rebinned the output. The raw-coverage route above
replaces it.
"""
import argparse
import csv
import datetime
import getpass
import hashlib
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
from pathlib import Path


COVERAGE_CONTRACT_VERSION = "BASIN_CANONICAL_COVERAGE_V1"
CANONICAL_BIN_SIZE = 1
CANONICAL_CENTER_READS = False
CANONICAL_EXTEND_READS = False
CANONICAL_QC_NORM = "CPM"


class TeeIO:
    """Write to both a stream (e.g., sys.stdout) and a file handle."""
    def __init__(self, stream, fh):
        self.stream = stream
        self.fh = fh

    def write(self, data):
        self.stream.write(data)
        self.fh.write(data)
        self.fh.flush()
        self.stream.flush()

    def flush(self):
        self.stream.flush()
        self.fh.flush()


def die(msg: str, code: int = 2):
    print(f"[HARMONY] ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def norm_colname(x: str) -> str:
    x = x.strip()
    if x.startswith("#"):
        x = x[1:].strip()
    if (x.startswith("'") and x.endswith("'")) or (x.startswith('"') and x.endswith('"')):
        x = x[1:-1]
    return x.strip()


def shlex_quote(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"


def sh(cmd, log_path=None):
    """Run a command and tee stdout/stderr into Python stdout."""
    if isinstance(cmd, list):
        printable = " ".join([shlex_quote(str(x)) for x in cmd])
    else:
        printable = cmd

    if log_path:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a") as f:
            f.write(printable + "\n")

    res = subprocess.run(
        cmd,
        shell=isinstance(cmd, str),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if res.stdout:
        print(res.stdout, end="")
    if res.returncode != 0:
        raise RuntimeError(f"Command failed ({res.returncode}): {printable}")
    return res.stdout


def read_yaml_minimal(path: Path) -> dict:
    """
    Minimal YAML reader for simple key: value configs used here.
    Supports booleans, ints, floats, strings, null, and simple lists [a,b].
    """
    cfg = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            continue

        key, val = line.split(":", 1)
        key = key.strip()
        val = val.strip()

        if "#" in val and not (val.startswith("'") or val.startswith('"')):
            val = val.split("#", 1)[0].strip()

        def parse_scalar(x):
            if x in ("null", "None", ""):
                return None
            if x.lower() in ("true", "false"):
                return x.lower() == "true"
            if re.fullmatch(r"-?\d+", x):
                return int(x)
            if re.fullmatch(r"-?\d+\.\d+", x):
                return float(x)
            if (x.startswith("'") and x.endswith("'")) or (x.startswith('"') and x.endswith('"')):
                return x[1:-1]
            return x

        if val.startswith("[") and val.endswith("]"):
            inner = val[1:-1].strip()
            if not inner:
                cfg[key] = []
            else:
                parts = [p.strip() for p in inner.split(",")]
                cfg[key] = [parse_scalar(p) for p in parts]
        else:
            cfg[key] = parse_scalar(val)

    return cfg


def write_text(path: Path, txt: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(txt)


def now_stamp():
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def get_slurm_jobid():
    return os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_ARRAY_JOB_ID")


def make_run_id(cfg: dict) -> str:
    rid = cfg.get("run_id", "auto")
    prefix = cfg.get("run_id_prefix", None)

    if rid and str(rid) != "auto":
        return str(rid)

    user = getpass.getuser()
    host = socket.gethostname().split(".")[0]
    job = get_slurm_jobid() or "NA"
    stamp = now_stamp()
    if not stamp:
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    if prefix and bool(cfg.get("append_run_stamp_to_prefix", True)):
        return f"{prefix}_{stamp}_{user}_{host}_{job}"

    return f"{stamp}_{user}_{host}_{job}"


def which_or_none(exe: str):
    return shutil.which(exe)


def tool_version(exe: str, args=("--version",)):
    p = which_or_none(exe)
    if not p:
        return f"{exe}\tNOT_FOUND"
    out = subprocess.run([exe, *args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout.strip()
    first = out.splitlines()[0] if out else ""
    return f"{exe}\t{first}"


def resolve_from_cfg(cfg_path: Path, value, kind="file") -> Path | None:
    """Resolve a path from YAML relative to the config file location."""
    if value in (None, "", "null"):
        return None
    p = Path(str(value)).expanduser()
    if not p.is_absolute():
        p = (cfg_path.parent / p).resolve()
    else:
        p = p.resolve()
    return p


def path_within_run(run_dir: Path, value, default_leaf: str) -> Path:
    """
    Paths like 'outputs/bigwigs_scaled' are intended to live under the run dir.
    Absolute paths are honored as-is.
    """
    if value in (None, "", "null"):
        return run_dir / default_leaf
    p = Path(str(value)).expanduser()
    if p.is_absolute():
        return p.resolve()
    return (run_dir / p).resolve()


def read_samples(samples_tsv: Path):
    rows = []
    with open(samples_tsv, newline="") as f:
        reader = csv.DictReader(filter(lambda r: not r.lstrip().startswith("#"), f), delimiter="\t")
        required = ["sample_id", "bam_path", "antibody", "genotype", "sex", "replicate"]
        for c in required:
            if c not in reader.fieldnames:
                die(f"samples.tsv missing required column: {c}")
        for r in reader:
            if not r.get("sample_id"):
                continue
            rows.append(r)
    if not rows:
        die("No samples found in samples.tsv")
    return rows


def count_dup_flagged(bam: Path) -> int:
    cmd = ["samtools", "view", "-c", "-f", "1024", str(bam)]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr.strip() or "samtools view failed")
    return int(out.stdout.strip() or "0")


def ensure_bai(bam: Path, commands_log: Path):
    bai1 = bam.with_suffix(bam.suffix + ".bai")  # sample.bam.bai
    bai2 = bam.with_suffix(".bai")               # sample.bai
    if bai1.exists() or bai2.exists():
        return
    sh(["samtools", "index", str(bam)], log_path=commands_log)


def markdup_samtools(in_bam: Path, out_bam: Path, commands_log: Path):
    """
    Conservative duplicate marking:
      sort by name -> fixmate -> sort by coord -> markdup -> index
    """
    tmp = out_bam.parent / f"{out_bam.stem}.tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    nsort = tmp / "name_sorted.bam"
    fixm = tmp / "fixmate.bam"
    psort = tmp / "pos_sorted.bam"

    sh(["samtools", "sort", "-n", "-o", str(nsort), str(in_bam)], log_path=commands_log)
    sh(["samtools", "fixmate", "-m", str(nsort), str(fixm)], log_path=commands_log)
    sh(["samtools", "sort", "-o", str(psort), str(fixm)], log_path=commands_log)
    sh(["samtools", "markdup", str(psort), str(out_bam)], log_path=commands_log)
    sh(["samtools", "index", str(out_bam)], log_path=commands_log)


def dedup_bam_flag1024(in_bam: Path, out_bam: Path, commands_log: Path):
    sh(["samtools", "view", "-b", "-F", "1024", "-o", str(out_bam), str(in_bam)], log_path=commands_log)
    sh(["samtools", "index", str(out_bam)], log_path=commands_log)


def validate_downstream_geometry(cfg: dict):
    """Validate the user-facing downstream coverage geometry contract.

    HARMONY defines these config flags by the *requested output representation*:
      center=false, extend=false -> aligned read footprints
      center=true,  extend=false -> read-length footprints moved to PE fragment centers
      center=false, extend=true  -> full fragment-spanning coverage
      center=true,  extend=true  -> unsupported/ambiguous by design

    deepTools requires --extendReads internally for the centered/read-length case so
    that paired-end fragment geometry is available to --centerReads. That internal
    implementation detail does not change the HARMONY output contract above.
    """
    center = bool(cfg.get("downstream_center_reads", True))
    extend = bool(cfg.get("downstream_extend_reads", False))
    if center and extend:
        die(
            "downstream_center_reads=true and downstream_extend_reads=true is intentionally unsupported. "
            "Use center=true/extend=false for read-length footprints centered on paired-end fragments, "
            "or center=false/extend=true for full fragment-spanning coverage."
        )
    return center, extend


def bamcoverage_cmd(
    bam: Path,
    bw: Path,
    cfg: dict,
    norm: str,
    scale_factor=None,
    frag_mode=None,
    representation="canonical",
):
    """Build a bamCoverage command under an explicit coverage representation.

    ``canonical`` is the immutable BASIN/CODA measurement model: 1-bp,
    uncentered, unextended aligned-read coverage. ``downstream`` is an optional
    derivative representation that may center and/or extend reads and may use a
    different output bin size. Downstream geometry never feeds BAR scaling or
    CODA.
    """
    if representation not in ("canonical", "downstream"):
        die(f"Unknown coverage representation: {representation}")

    if representation == "canonical":
        bin_size = CANONICAL_BIN_SIZE
        center_reads = CANONICAL_CENTER_READS
        extend_reads = CANONICAL_EXTEND_READS
    else:
        bin_size = int(cfg.get("downstream_bin_size", CANONICAL_BIN_SIZE))
        if bin_size < 1:
            die("downstream_bin_size must be >= 1")
        center_reads, extend_reads = validate_downstream_geometry(cfg)

    cmd = [
        "bamCoverage",
        "-b", str(bam),
        "-o", str(bw),
        "--binSize", str(bin_size),
        "--minMappingQuality", str(cfg.get("min_mapq", 0)),
    ]

    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    cmd += ["--numberOfProcessors", str(nproc)]

    # deepTools only has a fragment-length reference for --centerReads when
    # read extension is enabled. With --centerReads + --extendReads, deepTools
    # first infers the paired-end fragment and then returns a read-length
    # footprint centered on that fragment; it does NOT leave the full fragment
    # extended in the output. Therefore a BASIN "centered.unextended"
    # derivative intentionally passes both CLI flags.
    #
    # By contrast, downstream_extend_reads=true with centering disabled produces
    # genuinely fragment-spanning (fragment-length-weighted) coverage.
    if center_reads:
        cmd += ["--centerReads", "--extendReads"]
    elif extend_reads:
        cmd.append("--extendReads")

    if frag_mode == "lt":
        cmd += ["--maxFragmentLength", str(int(cfg.get("split_threshold_bp", 120)) - 1)]
    elif frag_mode == "gt":
        cmd += ["--minFragmentLength", str(int(cfg.get("split_threshold_bp", 120)))]

    norm = str(norm).upper()
    if norm == "CPM":
        cmd += ["--normalizeUsing", "CPM"]
    elif norm == "RPGC":
        eg = cfg.get("effective_genome_size_dm6")
        if not eg:
            die("RPGC output requested but effective_genome_size_dm6 is not set")
        cmd += ["--normalizeUsing", "RPGC", "--effectiveGenomeSize", str(eg)]
    elif norm == "NONE":
        pass
    else:
        die(f"Unknown normalization: {norm}")

    if scale_factor is not None:
        cmd += ["--scaleFactor", str(scale_factor)]

    return cmd


def multi_bigwig_summary(bws: list[Path], bed: Path, out_prefix: Path, commands_log: Path):
    out_npz = str(out_prefix) + ".npz"
    out_tab = str(out_prefix) + ".rawCounts.tsv"
    cmd = [
        "multiBigwigSummary", "BED-file",
        "--BED", str(bed),
        "--bwfiles", *[str(b) for b in bws],
        "--outFileName", out_npz,
        "--outRawCounts", out_tab,
    ]
    sh(cmd, log_path=commands_log)
    return Path(out_tab)


def median(vals):
    vals = sorted(vals)
    n = len(vals)
    if n == 0:
        return None
    mid = n // 2
    if n % 2 == 1:
        return float(vals[mid])
    return (float(vals[mid - 1]) + float(vals[mid])) / 2.0


def compute_scales_per_sample(
    raw_counts_tsv: Path,
    drop_fractions,
    out_scales_tsv: Path,
    retained_bars_tsv: Path | None = None,
):
    """Compute BASIN scales and optionally record the exact retained BAR rows.

    Contract shared with CODA:
      1. Remove rows that are zero across ALL samples.
      2. For each sample, rank remaining BARs by signal descending.
      3. n_drop = round(N * drop_fraction) using Python's round() semantics.
      4. Retain rows after the first n_drop and take their median.
    """
    with open(raw_counts_tsv, newline="") as f:
        reader = csv.reader(f, delimiter="	")
        header = next(reader)

        sample_idxs = [i for i, h in enumerate(header) if ".bw" in h]
        if not sample_idxs:
            die(f"Could not detect any bigWig columns in rawCounts header: {header}")

        sample_cols = [norm_colname(header[i]) for i in sample_idxs]
        per = {sample_col: [] for sample_col in sample_cols}
        parse_fail = {sample_col: 0 for sample_col in sample_cols}

        n_rows_total = 0
        n_rows_kept = 0
        n_rows_allzero = 0

        for row_order, row in enumerate(reader, start=1):
            n_rows_total += 1
            if len(row) <= max(sample_idxs):
                for sample_col in sample_cols:
                    parse_fail[sample_col] += 1
                n_rows_allzero += 1
                continue

            row_vals = []
            for sample_col, j in zip(sample_cols, sample_idxs):
                try:
                    v = float(row[j])
                except Exception:
                    parse_fail[sample_col] += 1
                    v = 0.0
                row_vals.append(v)

            if all(v == 0.0 for v in row_vals):
                n_rows_allzero += 1
                continue

            n_rows_kept += 1
            chrom = row[0] if len(row) > 0 else ""
            start = row[1] if len(row) > 1 else ""
            end = row[2] if len(row) > 2 else ""
            for sample_col, v in zip(sample_cols, row_vals):
                per[sample_col].append((row_order, chrom, start, end, v))

    out_scales_tsv.parent.mkdir(parents=True, exist_ok=True)
    retained_fh = None
    retained_writer = None
    if retained_bars_tsv is not None:
        retained_bars_tsv.parent.mkdir(parents=True, exist_ok=True)
        retained_fh = open(retained_bars_tsv, "w", newline="")
        retained_writer = csv.writer(retained_fh, delimiter="	")
        retained_writer.writerow([
            "sample_colname", "drop_fraction", "chr", "start", "end",
            "raw_signal", "rank_desc", "n_rows_after_allzero_filter", "n_trim",
        ])

    try:
        with open(out_scales_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="	")
            w.writerow([
                "sample_colname", "drop_fraction", "n_peaks_total",
                "n_peaks_used_after_allzero_filter", "n_kept_after_drop",
                "median_kept", "scale_factor", "n_parse_fail",
                "n_rows_allzero_filtered",
            ])

            for sample_col, entries in per.items():
                n_used = len(entries)
                entries_sorted = sorted(entries, key=lambda x: x[4], reverse=True)
                for d in drop_fractions:
                    d = float(d)
                    drop_n = int(round(n_used * d)) if n_used > 0 else 0
                    kept_entries = entries_sorted[drop_n:] if n_used > 0 else []
                    kept_vals = [x[4] for x in kept_entries]
                    med = median(kept_vals)
                    sf = "NA" if (med is None or med == 0) else (1.0 / med)
                    w.writerow([
                        sample_col, f"{d:.2f}", n_rows_total, n_rows_kept,
                        len(kept_entries), med if med is not None else "NA", sf,
                        parse_fail[sample_col], n_rows_allzero,
                    ])
                    if retained_writer is not None:
                        for rank_desc, entry in enumerate(kept_entries, start=drop_n + 1):
                            _, chrom, start, end, v = entry
                            retained_writer.writerow([
                                sample_col, f"{d:.2f}", chrom, start, end, v,
                                rank_desc, n_used, drop_n,
                            ])
    finally:
        if retained_fh is not None:
            retained_fh.close()


def write_versions(path: Path):
    lines = [
        tool_version("samtools", ("--version",)),
        tool_version("bamCoverage", ("--version",)),
        tool_version("multiBigwigSummary", ("--version",)),
        "python\t" + sys.version.splitlines()[0],
    ]
    write_text(path, "\n".join(lines) + "\n")


def write_environment(path: Path):
    keys = [
        "SLURM_JOB_ID",
        "SLURM_ARRAY_JOB_ID",
        "SLURM_ARRAY_TASK_ID",
        "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE",
        "TMPDIR",
    ]
    lines = [
        f"date\t{datetime.datetime.now().isoformat()}",
        f"user\t{getpass.getuser()}",
        f"host\t{socket.gethostname()}",
        f"pwd\t{os.getcwd()}",
        f"platform\t{platform.platform()}",
    ]
    for k in keys:
        if os.environ.get(k):
            lines.append(f"{k}\t{os.environ.get(k)}")
    write_text(path, "\n".join(lines) + "\n")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_key_value_tsv(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    with open(path, newline="") as f:
        for row in csv.reader(f, delimiter="	"):
            if len(row) < 2 or not row[0] or row[0].startswith("#"):
                continue
            if row[0] == "key" and row[1] == "value":
                continue
            out[row[0]] = row[1]
    return out


def write_coverage_contract(path: Path, cfg: dict, universal_bed: Path, dedup_mode: str):
    rows = [
        ("contract_version", COVERAGE_CONTRACT_VERSION),
        ("canonical_signal", "aligned_read_coverage"),
        ("canonical_bin_size", str(CANONICAL_BIN_SIZE)),
        ("canonical_center_reads", str(CANONICAL_CENTER_READS).lower()),
        ("canonical_extend_reads", str(CANONICAL_EXTEND_READS).lower()),
        ("canonical_qc_norm", CANONICAL_QC_NORM),
        ("min_mapq", str(int(cfg.get("min_mapq", 0)))),
        ("dedup_mode", dedup_mode),
        ("bar_file", str(universal_bed)),
        ("bar_sha256", sha256_file(universal_bed)),
        ("allzero_rule", "remove_BAR_rows_zero_across_all_samples_before_trimming"),
        ("trim_rule", "n_trim=round(N*drop_fraction);stable_descending_signal_order"),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="	")
        w.writerow(["key", "value"])
        w.writerows(rows)


def downstream_geometry_tag(cfg: dict) -> str:
    center = bool(cfg.get("downstream_center_reads", True))
    extend = bool(cfg.get("downstream_extend_reads", False))
    center_tag = "centered" if center else "uncentered"
    extend_tag = "extended_fragmentweighted" if extend else "unextended"
    return f"{center_tag}.{extend_tag}"


def build_sample_contexts(samples, cfg_path: Path):
    contexts = []
    for r in samples:
        bam = resolve_from_cfg(cfg_path, r["bam_path"])
        if bam is None:
            die(f"Missing bam_path for sample {r['sample_id']}")
        contexts.append({
            "row": r,
            "sample_id": r["sample_id"],
            "original_bam": bam,
            "prepared_bam": None,
            "prepared_mode": None,
        })
    return contexts


def prepare_bam_for_sample(ctx: dict, cfg: dict, prep_bam_dir: Path, commands_log: Path, reuse_prep_dirs: list | None = None) -> Path:
    """
    Lazily prepare the BAM actually required for coverage/scaling for this sample.

    Reuse order:
      1. already resolved in current context
      2. matching prepared BAM already present in THIS run's bams_prepared/
      3. matching prepared BAM present in any EXTERNAL reuse directory, in the
         order given by reuse_prep_dirs (e.g. the config 'bams_prepared' path
         first, then a prior run's bams_prepared/ if reuse_baseline_from_run is
         also set)
      4. build a new prepared BAM in THIS run

    Reusing prepared BAMs is the expensive step to avoid: markdup + dedup are
    skipped entirely whenever a matching prepared BAM is found. Point the config
    'bams_prepared' key at any directory of prepared BAMs (for example a prior
    HARMONY run's bams_prepared/ directory) to reuse them. Filenames must follow
    the convention below, derived from the ORIGINAL bam stem:
        keep   mode:  <stem>.markdup.bam
        remove mode:  <stem>.markdup.dedup.bam  or  <stem>.dedup.bam
    """
    if ctx["prepared_bam"] is not None:
        return ctx["prepared_bam"]

    bam = ctx["original_bam"]
    sid = ctx["sample_id"]

    if not bam.exists():
        die(f"BAM not found for sample {sid}: {bam}")

    dedup_mode = str(cfg.get("dedup_mode", "remove")).lower()
    if dedup_mode not in ("remove", "keep"):
        die("dedup_mode must be 'remove' or 'keep'")

    reuse_prep_dirs = reuse_prep_dirs or []

    # Expected prepared BAM names derived from the ORIGINAL bam stem
    current_markdup = prep_bam_dir / f"{bam.stem}.markdup.bam"
    current_dedup_from_markdup = prep_bam_dir / f"{bam.stem}.markdup.dedup.bam"
    current_dedup_direct = prep_bam_dir / f"{bam.stem}.dedup.bam"

    # ------------------------------------------------------------
    # KEEP mode
    # ------------------------------------------------------------
    if dedup_mode == "keep":
        # Reuse current-run prepared BAM if present
        if current_markdup.exists():
            print(f"[HARMONY] Reusing current-run markdup BAM for sample {sid}: {current_markdup}")
            ctx["prepared_bam"] = current_markdup
            ctx["prepared_mode"] = "keepdup"
            return ctx["prepared_bam"]

        # Reuse an external prepared BAM if present (bams_prepared / prior run)
        for ext_dir in reuse_prep_dirs:
            ext_markdup = ext_dir / f"{bam.stem}.markdup.bam"
            if ext_markdup.exists():
                print(f"[HARMONY] Reusing external markdup BAM for sample {sid}: {ext_markdup}")
                ctx["prepared_bam"] = ext_markdup
                ctx["prepared_mode"] = "keepdup_reused_external"
                return ctx["prepared_bam"]

        # Otherwise, if duplicate flags already exist in the original BAM, original BAM is valid
        if bool(cfg.get("require_bai", True)):
            ensure_bai(bam, commands_log)

        dup_count = count_dup_flagged(bam)
        if dup_count > 0:
            print(f"[HARMONY] Sample {sid} already has duplicate flags; using original BAM in keep mode: {bam}")
            ctx["prepared_bam"] = bam
            ctx["prepared_mode"] = "original_with_dupflags"
            return ctx["prepared_bam"]

        # Otherwise build markdup only if requested; else original BAM is still the keep-mode input
        if bool(cfg.get("mark_duplicates_if_missing", True)):
            tool = str(cfg.get("markdup_tool", "samtools")).lower()
            if tool != "samtools":
                die("picard markdup not implemented in this runner (set markdup_tool: samtools).")
            print(f"[HARMONY] Building markdup BAM for sample {sid}: {bam} -> {current_markdup}")
            markdup_samtools(bam, current_markdup, commands_log)
            ctx["prepared_bam"] = current_markdup
            ctx["prepared_mode"] = "keepdup_built"
            return ctx["prepared_bam"]

        print(f"[HARMONY] Using original BAM for sample {sid} in keep mode: {bam}")
        ctx["prepared_bam"] = bam
        ctx["prepared_mode"] = "original"
        return ctx["prepared_bam"]

    # ------------------------------------------------------------
    # REMOVE mode
    # ------------------------------------------------------------
    # Prefer current-run dedup BAMs
    for cand in (current_dedup_from_markdup, current_dedup_direct):
        if cand.exists():
            print(f"[HARMONY] Reusing current-run deduplicated BAM for sample {sid}: {cand}")
            ctx["prepared_bam"] = cand
            ctx["prepared_mode"] = "dedup"
            return ctx["prepared_bam"]

    # Then external dedup BAMs (bams_prepared / prior run), in priority order
    for ext_dir in reuse_prep_dirs:
        for cand in (ext_dir / f"{bam.stem}.markdup.dedup.bam", ext_dir / f"{bam.stem}.dedup.bam"):
            if cand.exists():
                print(f"[HARMONY] Reusing external deduplicated BAM for sample {sid}: {cand}")
                ctx["prepared_bam"] = cand
                ctx["prepared_mode"] = "dedup_reused_external"
                return ctx["prepared_bam"]

    # If we get here, we need to build a new prepared BAM
    if bool(cfg.get("require_bai", True)):
        ensure_bai(bam, commands_log)

    dup_count = count_dup_flagged(bam)
    bam_mark = bam

    if dup_count == 0 and bool(cfg.get("mark_duplicates_if_missing", True)):
        tool = str(cfg.get("markdup_tool", "samtools")).lower()
        if tool != "samtools":
            die("picard markdup not implemented in this runner (set markdup_tool: samtools).")
        print(f"[HARMONY] Building markdup BAM for sample {sid}: {bam} -> {current_markdup}")
        markdup_samtools(bam, current_markdup, commands_log)
        bam_mark = current_markdup
    elif dup_count > 0:
        print(f"[HARMONY] Sample {sid} already has duplicate flags; using original BAM as marked input: {bam}")

    # Decide final dedup output name based on what bam_mark actually is
    if bam_mark.name.endswith(".markdup.bam"):
        final_dedup = prep_bam_dir / f"{bam.stem}.markdup.dedup.bam"
    else:
        final_dedup = prep_bam_dir / f"{bam.stem}.dedup.bam"

    if not final_dedup.exists() or bool(cfg.get("allow_overwrite", False)):
        print(f"[HARMONY] Building deduplicated BAM for sample {sid}: {bam_mark} -> {final_dedup}")
        dedup_bam_flag1024(bam_mark, final_dedup, commands_log)
    else:
        print(f"[HARMONY] Reusing current-run deduplicated BAM for sample {sid}: {final_dedup}")

    ctx["prepared_bam"] = final_dedup
    ctx["prepared_mode"] = "dedup_built"
    return ctx["prepared_bam"]

def load_scales(scales_tsv: Path):
    """
    Load scale information, keeping NA samples visible.
    Returns:
      scales_by_sample_col: dict[str, dict[str, float|None]]
    """
    scales = {}
    with open(scales_tsv, newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            sample_col = norm_colname(row["sample_colname"])
            drop = f"{float(row['drop_fraction']):.2f}"
            raw = row.get("scale_factor")
            sf = None
            if raw not in (None, "", "NA"):
                try:
                    sf = float(raw)
                except ValueError:
                    sf = None
            scales.setdefault(sample_col, {})[drop] = sf
    return scales


def match_scale_key(scales_norm: dict, bw_path: Path):
    bw_abs = norm_colname(str(bw_path.resolve()))
    bw_str = norm_colname(str(bw_path))
    bw_base = norm_colname(bw_path.name)

    for cand in (bw_abs, bw_str, bw_base):
        if cand in scales_norm:
            return cand

    matches = [
        k for k in scales_norm.keys()
        if k.endswith("/" + bw_base) or k.endswith("\\" + bw_base) or k == bw_base
    ]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError(f"Multiple matching scale columns for basename {bw_base}: {matches[:5]} ...")
    return None


def main():
    ap = argparse.ArgumentParser(description="HARMONY normalization pipeline.")
    ap.add_argument("--config", required=True, help="Path to config.yaml")
    args = ap.parse_args()

    cfg_path = Path(args.config).resolve()
    if not cfg_path.exists():
        die(f"Config not found: {cfg_path}")
    cfg = read_yaml_minimal(cfg_path)

    # Legacy global geometry keys are accepted only when they equal the
    # canonical BASIN contract. They no longer control canonical coverage.
    if cfg.get("bin_size") not in (None, CANONICAL_BIN_SIZE):
        die("Legacy bin_size can no longer change canonical BASIN coverage; use downstream_bin_size.")
    if bool(cfg.get("center_reads", False)):
        die("Legacy center_reads=true is no longer allowed on the canonical path; use downstream_center_reads.")
    if bool(cfg.get("extend_reads", False)):
        die("Legacy extend_reads=true is unsafe for BAR scaling; use downstream_extend_reads only for derivatives.")
    legacy_baseline_norm = cfg.get("baseline_track_norm")
    if legacy_baseline_norm is not None and str(legacy_baseline_norm).upper() != CANONICAL_QC_NORM:
        die("baseline_track_norm is no longer configurable for CODA input; canonical baseline QC is always CPM.")
    if cfg.get("keep_cpm_bigwigs") is False:
        die("Canonical CPM QC bigWigs are required by CODA and can no longer be disabled.")

    samples_path = resolve_from_cfg(cfg_path, cfg.get("samples_tsv", "samples.tsv"))
    if samples_path is None or not samples_path.exists():
        die(f"samples.tsv not found: {samples_path}")
    samples = read_samples(samples_path)
    sample_contexts = build_sample_contexts(samples, cfg_path)

    universal_bed = resolve_from_cfg(cfg_path, cfg.get("universal_bed"))
    if universal_bed is None or not universal_bed.exists():
        die(f"universal_bed not found: {universal_bed}")

    run_root = resolve_from_cfg(cfg_path, cfg.get("run_root", "outputs/runs"), kind="dir")
    if run_root is None:
        die("run_root could not be resolved")
    run_id = make_run_id(cfg)
    run_dir = run_root / run_id
    if run_dir.exists() and not bool(cfg.get("allow_overwrite", False)):
        die(f"Run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)

    baseline_bw_dir = path_within_run(run_dir, cfg.get("baseline_bw_dir"), "bigwigs_baseline")
    scaled_bw_dir = path_within_run(run_dir, cfg.get("scaled_bw_dir"), "bigwigs_scaled")
    downstream_bw_dir = path_within_run(run_dir, cfg.get("downstream_bw_dir"), "bigwigs_downstream")
    scales_dir = path_within_run(run_dir, cfg.get("scales_dir"), "scales")
    logs_dir = path_within_run(run_dir, cfg.get("logs_dir"), "logs")
    signal_summary_dir = path_within_run(run_dir, cfg.get("signal_summary_dir"), "signal_summary")
    prep_bam_dir = path_within_run(run_dir, cfg.get("prepared_bam_dir"), "bams_prepared")

    for d in [baseline_bw_dir, scaled_bw_dir, downstream_bw_dir, scales_dir, logs_dir, signal_summary_dir, prep_bam_dir]:
        d.mkdir(parents=True, exist_ok=True)

    if bool(cfg.get("tee_run_log", True)):
        run_log_path = logs_dir / "harmony_run.log"
        run_log_fh = open(run_log_path, "a")
        sys.stdout = TeeIO(sys.stdout, run_log_fh)
        sys.stderr = TeeIO(sys.stderr, run_log_fh)
        print(f"[HARMONY] TEE logging to: {run_log_path}")

    commands_log = run_dir / "run_manifest" / "commands.sh"
    commands_log.parent.mkdir(parents=True, exist_ok=True)
    write_text(commands_log, "#!/usr/bin/env bash\nset -euo pipefail\n\n")

    # Optional prior run baseline reuse
    reuse_baseline_from_run = resolve_from_cfg(cfg_path, cfg.get("reuse_baseline_from_run"))
    old_baseline_dir = None
    if reuse_baseline_from_run is not None:
        old_baseline_dir = reuse_baseline_from_run / "bigwigs_baseline"
        if not old_baseline_dir.exists():
            die(f"reuse_baseline_from_run was set, but baseline dir not found: {old_baseline_dir}")
        print(f"[HARMONY] Reusing baseline bigWigs from prior run: {reuse_baseline_from_run}")

    # External prepared-BAM reuse directories, searched in priority order by
    # prepare_bam_for_sample. The config 'bams_prepared' path is checked first
    # (this is the cluster-friendly way to avoid re-running markdup/dedup: point
    # it at any directory of prepared BAMs, e.g. a prior run's bams_prepared/),
    # then a prior run's bams_prepared/ if reuse_baseline_from_run is also set.
    reuse_prep_dirs = []

    bams_prepared_cfg = resolve_from_cfg(cfg_path, cfg.get("bams_prepared"), kind="dir")
    if bams_prepared_cfg is not None:
        if not bams_prepared_cfg.exists():
            die(f"bams_prepared was set, but directory not found: {bams_prepared_cfg}")
        print(f"[HARMONY] Reusing prepared BAMs from bams_prepared when available: {bams_prepared_cfg}")
        reuse_prep_dirs.append(bams_prepared_cfg)

    old_prep_bam_dir = None
    if reuse_baseline_from_run is not None:
        old_prep_bam_dir = reuse_baseline_from_run / "bams_prepared"
        if old_prep_bam_dir.exists():
            print(f"[HARMONY] Reusing prepared BAMs from prior run when available: {old_prep_bam_dir}")
            reuse_prep_dirs.append(old_prep_bam_dir)
        else:
            print(f"[HARMONY] Prior run has no bams_prepared directory; prepared BAMs will be rebuilt if needed.")

    if bool(cfg.get("write_run_manifest", True)):
        resolved_cfg = dict(cfg)
        resolved_cfg["samples_tsv"] = str(samples_path)
        resolved_cfg["universal_bed"] = str(universal_bed)
        resolved_cfg["run_dir"] = str(run_dir)
        if reuse_baseline_from_run is not None:
            resolved_cfg["reuse_baseline_from_run"] = str(reuse_baseline_from_run)
        cfg_lines = [f"{k}: {v}" for k, v in resolved_cfg.items()]
        write_text(run_dir / "run_manifest" / "config.resolved.yaml", "\n".join(cfg_lines) + "\n")

        with open(run_dir / "run_manifest" / "samples.resolved.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            cols = ["sample_id", "bam_path", "antibody", "genotype", "sex", "replicate"]
            w.writerow(cols)
            for r in samples:
                w.writerow([r.get(c, "") for c in cols])

        write_versions(run_dir / "run_manifest" / "versions.txt")
        write_environment(run_dir / "run_manifest" / "environment.txt")

    # Preflight
    for exe in ["samtools", "bamCoverage", "multiBigwigSummary"]:
        if not which_or_none(exe):
            die(f"Required tool not found in PATH: {exe}")

    dedup_mode = str(cfg.get("dedup_mode", "remove")).lower()
    if dedup_mode not in ("remove", "keep"):
        die("dedup_mode must be 'remove' or 'keep'")
    dup_tag = "dedup" if dedup_mode == "remove" else "keepdup"

    make_downstream_cpm = bool(cfg.get("make_downstream_cpm_bigwigs", False))
    make_downstream_scaled = bool(cfg.get("make_downstream_scaled_bigwigs", False))
    downstream_center = bool(cfg.get("downstream_center_reads", True))
    downstream_extend = bool(cfg.get("downstream_extend_reads", False))
    if make_downstream_cpm or make_downstream_scaled:
        downstream_center, downstream_extend = validate_downstream_geometry(cfg)
    if downstream_extend and (make_downstream_cpm or make_downstream_scaled):
        print(
            "[HARMONY] WARNING: downstream_extend_reads=true. Derivative filenames will be labeled "
            "extended_fragmentweighted. Extension makes integrated coverage depend on inferred fragment length."
        )

    contract_path = run_dir / "run_manifest" / "coverage_contract.tsv"
    expected_contract = {
        "contract_version": COVERAGE_CONTRACT_VERSION,
        "canonical_bin_size": str(CANONICAL_BIN_SIZE),
        "canonical_center_reads": "false",
        "canonical_extend_reads": "false",
        "canonical_qc_norm": CANONICAL_QC_NORM,
        "min_mapq": str(int(cfg.get("min_mapq", 0))),
        "dedup_mode": dedup_mode,
    }

    # A resumed/overwritten run may already contain canonical bigWigs. Never
    # trust those merely because their filenames match: validate the contract
    # that existed before this invocation.
    preexisting_contract = read_key_value_tsv(contract_path) if contract_path.exists() else {}
    preexisting_canonical_files = any(baseline_bw_dir.glob("*.NONE.bw")) or any(baseline_bw_dir.glob("*.CPM.bw"))
    if preexisting_canonical_files:
        if not preexisting_contract:
            if not bool(cfg.get("allow_unverified_legacy_baseline_reuse", False)):
                die(
                    "Current run directory already contains canonical baseline bigWigs but no prior coverage_contract.tsv. "
                    "Use a new run_id, delete/rebuild those files, or set allow_unverified_legacy_baseline_reuse: true only "
                    "after manually verifying their geometry/MAPQ/dedup settings."
                )
            print("[HARMONY] WARNING: reusing current-run legacy baseline bigWigs without a machine-readable coverage contract.")
        else:
            mismatches = [
                f"{k}: existing={preexisting_contract.get(k)!r}, requested={v!r}"
                for k, v in expected_contract.items()
                if preexisting_contract.get(k) != v
            ]
            if mismatches:
                die("Existing current-run baseline coverage contract does not match this invocation:\n  " + "\n  ".join(mismatches))

    write_coverage_contract(contract_path, cfg, universal_bed, dedup_mode)

    if reuse_baseline_from_run is not None:
        old_contract_path = reuse_baseline_from_run / "run_manifest" / "coverage_contract.tsv"
        old_contract = read_key_value_tsv(old_contract_path)
        if not old_contract:
            if not bool(cfg.get("allow_unverified_legacy_baseline_reuse", False)):
                die(
                    "Prior run has no coverage_contract.tsv. Set allow_unverified_legacy_baseline_reuse: true only "
                    "after manually confirming binSize=1, no centering, no extension, the same min_mapq and dedup_mode."
                )
            print("[HARMONY] WARNING: reusing legacy baseline bigWigs without a machine-readable coverage contract.")
        else:
            mismatches = [
                f"{k}: prior={old_contract.get(k)!r}, current={v!r}"
                for k, v in expected_contract.items()
                if old_contract.get(k) != v
            ]
            if mismatches:
                die("Prior-run baseline coverage contract does not match this run:\n  " + "\n  ".join(mismatches))

    # ------------------------------------------------------------------
    # Helper: obtain a baseline bigWig for a given normalization, honoring
    # the current-run and prior-run reuse policy. Returns the bigWig path,
    # or None if it could not be produced.
    # ------------------------------------------------------------------
    def obtain_baseline_bw(ctx, norm_label):
        sid = ctx["sample_id"]
        bw = baseline_bw_dir / f"{sid}.{norm_label}.bw"

        if bw.exists() and bool(cfg.get("skip_baseline_if_present", True)):
            print(f"[HARMONY] Reusing existing {norm_label} bigWig in current run: {bw}")
            return bw, "reused_current_run"

        if old_baseline_dir is not None:
            old_bw = old_baseline_dir / f"{sid}.{norm_label}.bw"
            if old_bw.exists():
                if bool(cfg.get("copy_reused_baseline_bw", True)):
                    shutil.copy2(old_bw, bw)
                    print(f"[HARMONY] Copied {norm_label} bigWig from prior run: {old_bw} -> {bw}")
                    return bw, f"copied_from:{old_bw}"
                else:
                    if bw.exists() or bw.is_symlink():
                        bw.unlink()
                    os.symlink(old_bw, bw)
                    print(f"[HARMONY] Symlinked {norm_label} bigWig from prior run: {old_bw} -> {bw}")
                    return bw, f"symlinked_from:{old_bw}"

        bam_use = prepare_bam_for_sample(ctx, cfg, prep_bam_dir, commands_log, reuse_prep_dirs)
        print(f"[HARMONY] Generating {norm_label} bigWig for sample {sid} from BAM: {bam_use}")
        sh(bamcoverage_cmd(bam_use, bw, cfg, norm=norm_label), log_path=commands_log)
        return bw, f"generated_from_bam:{bam_use}"

    # ------------------------------------------------------------------
    # Canonical baseline bigWigs
    #   * NONE (raw) -> BASIN scale basis
    #   * CPM        -> analytical baseline/QC input for CODA
    # Both are ALWAYS 1-bp, uncentered, and unextended.
    # ------------------------------------------------------------------
    scale_basis_paths = []
    cpm_qc_paths = []
    provenance_rows = []

    for ctx in sample_contexts:
        sid = ctx["sample_id"]
        none_bw, none_note = obtain_baseline_bw(ctx, "NONE")
        if none_bw.exists():
            scale_basis_paths.append(none_bw)
        cpm_bw, cpm_note = obtain_baseline_bw(ctx, "CPM")
        if cpm_bw.exists():
            cpm_qc_paths.append(cpm_bw)

        downstream_cpm_bw = ""
        if make_downstream_cpm:
            bam_use = prepare_bam_for_sample(ctx, cfg, prep_bam_dir, commands_log, reuse_prep_dirs)
            geom_tag = downstream_geometry_tag(cfg)
            ddir = downstream_bw_dir / "CPM"
            ddir.mkdir(parents=True, exist_ok=True)
            dbw = ddir / f"{sid}.CPM.{geom_tag}.bw"
            print(f"[HARMONY] Generating downstream CPM bigWig ({geom_tag}) for {sid}: {dbw}")
            sh(bamcoverage_cmd(bam_use, dbw, cfg, norm="CPM", representation="downstream"), log_path=commands_log)
            downstream_cpm_bw = str(dbw)

        provenance_rows.append([
            sid, str(ctx["original_bam"]), str(ctx["prepared_bam"]) if ctx["prepared_bam"] is not None else "",
            ctx["prepared_mode"] or "", str(none_bw) if none_bw.exists() else "", none_note or "",
            str(cpm_bw) if cpm_bw.exists() else "", cpm_note or "", downstream_cpm_bw,
        ])

    with open(run_dir / "run_manifest" / "sample_provenance.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="	")
        w.writerow([
            "sample_id", "original_bam", "prepared_bam", "prepared_mode",
            "scale_basis_bigwig", "scale_basis_note", "canonical_cpm_qc_bigwig",
            "canonical_cpm_qc_note", "downstream_cpm_bigwig",
        ])
        w.writerows(provenance_rows)

    if not scale_basis_paths:
        die("No canonical scale-basis (NONE) bigWigs are available for summarization.")
    if not cpm_qc_paths:
        die("No canonical CPM QC bigWigs are available for CODA.")

    none_summary_prefix = signal_summary_dir / "universal_peaks.NONE"
    none_counts_tsv = multi_bigwig_summary(scale_basis_paths, universal_bed, none_summary_prefix, commands_log)
    cpm_summary_prefix = signal_summary_dir / "universal_peaks"
    multi_bigwig_summary(cpm_qc_paths, universal_bed, cpm_summary_prefix, commands_log)

    # ------------------------------------------------------------------
    # Compute scales from the RAW (NONE) BAR coverage
    # ------------------------------------------------------------------
    drop_fracs = cfg.get("drop_fractions", [0.10, 0.20, 0.35, 0.50])
    scales_tsv = scales_dir / "scales_long.tsv"
    retained_bars_tsv = scales_dir / "retained_bars.tsv"
    compute_scales_per_sample(none_counts_tsv, drop_fracs, scales_tsv, retained_bars_tsv=retained_bars_tsv)

    scales_norm = load_scales(scales_tsv)

    # ------------------------------------------------------------------
    # Scaled bigWigs: produced from the BAM via
    #   bamCoverage --normalizeUsing None --scaleFactor <sf>
    # so output = raw_coverage * sf = raw / median(raw_BAR). Full resolution,
    # depth-invariant, no bigwigCompare.
    # ------------------------------------------------------------------
    for ctx in sample_contexts:
        sid = ctx["sample_id"]
        # Scale factors are keyed by the NONE bigWig column from multiBigwigSummary.
        none_bw = baseline_bw_dir / f"{sid}.NONE.bw"

        if not none_bw.exists():
            print(f"[HARMONY] WARNING: skipping sample {sid} because scale-basis bigWig is missing: {none_bw}")
            continue

        try:
            key = match_scale_key(scales_norm, none_bw)
        except RuntimeError as e:
            die(str(e))

        if key is None:
            if bool(cfg.get("skip_samples_with_na_scale", True)):
                print(
                    f"[HARMONY] WARNING: skipping sample {sid} because no matching scale entry was found.\n"
                    f"  Scale-basis bigWig: {none_bw}\n"
                    f"  Check rawCounts header in {none_counts_tsv}"
                )
                continue
            die(f"Could not find scale column for sample {sid} (scale-basis bigWig: {none_bw})")

        for d in drop_fracs:
            d_str = f"{float(d):.2f}"
            sf = scales_norm.get(key, {}).get(d_str)

            if sf is None:
                if bool(cfg.get("skip_samples_with_na_scale", True)):
                    print(f"[HARMONY] WARNING: skipping scaled bigWig for sample {sid}, drop {d_str} because scale_factor is NA")
                    continue
                die(f"Missing scale factor for sample {sid}, drop {d_str}")

            outdir_d = scaled_bw_dir / f"drop{d_str}"
            outdir_d.mkdir(parents=True, exist_ok=True)

            bam_use = prepare_bam_for_sample(ctx, cfg, prep_bam_dir, commands_log, reuse_prep_dirs)

            out_bw = outdir_d / f"{sid}.scaled.drop{d_str}.{dup_tag}.bw"
            print(
                f"[HARMONY] Generating scaled bigWig for sample {sid}, drop {d_str}, "
                f"scale {sf} from BAM (norm=NONE): {bam_use} -> {out_bw}"
            )
            sh(bamcoverage_cmd(bam_use, out_bw, cfg, norm="NONE", scale_factor=sf, representation="canonical"), log_path=commands_log)

            if make_downstream_scaled:
                geom_tag = downstream_geometry_tag(cfg)
                ddir = downstream_bw_dir / "BASIN" / f"drop{d_str}"
                ddir.mkdir(parents=True, exist_ok=True)
                dbw = ddir / f"{sid}.scaled.drop{d_str}.{dup_tag}.{geom_tag}.bw"
                print(f"[HARMONY] Generating downstream BASIN bigWig ({geom_tag}) for {sid}, drop {d_str}: {dbw}")
                sh(
                    bamcoverage_cmd(bam_use, dbw, cfg, norm="NONE", scale_factor=sf, representation="downstream"),
                    log_path=commands_log,
                )

            if bool(cfg.get("make_fragment_splits", False)):
                thr = int(cfg.get("split_threshold_bp", 120))
                bw_lt = outdir_d / f"{sid}.scaled.drop{d_str}.{dup_tag}.lt{thr}.bw"
                bw_gt = outdir_d / f"{sid}.scaled.drop{d_str}.{dup_tag}.gt{thr}.bw"
                print(f"[HARMONY] Generating fragment-split scaled bigWigs for sample {sid}, drop {d_str}")
                sh(bamcoverage_cmd(bam_use, bw_lt, cfg, norm="NONE", scale_factor=sf, frag_mode="lt"), log_path=commands_log)
                sh(bamcoverage_cmd(bam_use, bw_gt, cfg, norm="NONE", scale_factor=sf, frag_mode="gt"), log_path=commands_log)

    if bool(cfg.get("write_checksums", False)):
        chk = run_dir / "run_manifest" / "checksums.txt"
        lines = [
            f"{samples_path}\t{sha256_file(samples_path)}",
            f"{universal_bed}\t{sha256_file(universal_bed)}",
        ]
        write_text(chk, "\n".join(lines) + "\n")

    print(f"[HARMONY] DONE. Run directory: {run_dir}")


if __name__ == "__main__":
    main()
