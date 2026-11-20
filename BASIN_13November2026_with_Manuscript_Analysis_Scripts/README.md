# BASIN — Background-Anchored Signal Normalization with Integrated Quality Control for CUT&RUN and CUT&Tag

**Package release:** 2026-11-13. This package includes the 2026-09-23
coverage-contract patch; see `BASIN_VALIDATION_2026-09-23.txt` for its
validation record. End-to-end validation of this release on real data is
pending.

BASIN is a normalization framework for CUT&RUN and CUT&TAG chromatin profiling
data. It addresses the failure of CPM normalization in experiments where global
chromatin state differs across conditions — for example, histone mutants, drug
treatments, or developmental comparisons — by anchoring normalization to genomic
regions of consistent nonspecific IgG enrichment called BASIN Anchor Regions
(BARs). It can be used either as a replacement for spike-ins, or as an 
orthogonal method that is compatible with spike-ins. 

The key biological insight is that stronger antibody enrichment reduces
nonspecific IgG binding — so IgG signal is inversely related to true signal
strength, and can be used as a normalization reference without adding any
exogenous material. BASIN requires no spike-in, no additional library
preparation, and no paired IgG control beyond what is standard for CUT&RUN and
CUT&TAG experiments.

---

## The three pipelines

**OVERTURE** — BAR discovery (`OVERTURE/`)

Discovers BAR files from IgG CUT&RUN or CUT&TAG datasets. Supports universal
panel discovery (from diverse public IgG datasets) and experiment-specific
discovery (from your own IgG controls). Pre-built universal BAR files for
*Drosophila melanogaster* (dm6) are distributed with BASIN — most users will
not need to run OVERTURE at all. Mammalian anchor panels are currently in 
development. 

**HARMONY** — Normalization (`HARMONY/`)

The core pipeline. Takes aligned BAM files and a BAR file as input and produces
canonical BASIN-scaled bigWigs alongside a canonical CPM track used by CODA. The
BAR scale basis is fixed to 1-bp, uncentered, unextended aligned-read coverage;
this geometry is not user-configurable. The scale factor is derived from raw
IgG-background coverage at the BARs and applied to raw coverage, so normalization
is depth-invariant (identical samples sequenced at different depths produce
identical output). HARMONY can additionally generate centered and/or extended
**downstream derivative** tracks in a separate directory. Those derivatives are
never used for BAR scaling or CODA QC. Handles duplicate marking, scale factor
computation across multiple drop fractions, exact retained-BAR recording, and
full run provenance logging. Prepared (markdup/dedup) BAMs from a prior run can
be reused via the `bams_prepared` config key to skip the expensive duplicate-
handling step. This is the pipeline most users will run.

**CODA** — Post-normalization QC (`CODA/`)

Evaluates whether normalization succeeded, flags individual samples that may be
problematic, and issues experiment-level recommendations (pass / warning /
manual inspection required / hard fail). CODA is run from inside a completed
HARMONY run directory — the same directory containing `bigwigs_scaled/`,
`bigwigs_baseline/`, and `scales/`. It creates a `qc_[label]/` subdirectory
there with all outputs. CODA operates only on HARMONY's **canonical** CPM and
canonical BASIN-scaled bigWigs; optional centered/extended derivatives are ignored.
Current HARMONY runs also provide CODA the exact BAR identities retained at each
trim threshold. The two primary summaries (`experiment_qc_summary.tsv`,
`sample_qc_summary.tsv`) sit
at the top of `qc_[label]/`; all supporting tables and figures are written into a
`secondary_outputs/` subdirectory.

---

## Canonical coverage contract

BASIN separates the coverage representation used to **define and measure the
normalization reference** from optional representations used for downstream
visualization or analysis. The canonical BAR-normalization/QC path is:

- **OVERTURE BAR discovery:** 1-bp, uncentered, unextended CPM coverage.
- **HARMONY scale basis:** 1-bp, uncentered, unextended raw (`NONE`) coverage.
- **HARMONY CODA baseline:** 1-bp, uncentered, unextended CPM coverage.
- **HARMONY canonical BASIN output:** 1-bp, uncentered, unextended scaled raw coverage.

Optional downstream coverage tracks. HARMONY can also generate additional CPM 
and BASIN-scaled BigWigs with user-configurable centering, read extension, and 
bin size. Enable these outputs using make_downstream_cpm_bigwigs and/or 
make_downstream_scaled_bigwigs, and specify the desired geometry with 
downstream_center_reads, downstream_extend_reads, and downstream_bin_size. These 
derivative tracks are intended for downstream visualization and analysis; they 
do not affect BAR quantification, normalization factors, canonical BASIN outputs, 
or CODA quality control.

`min_mapq` remains configurable (default 0) and is recorded in run provenance;
it must be applied consistently to the canonical products within a HARMONY run.
The distributed dm6 panels were discovered at MAPQ 0.

Centering and extension are deliberately separated at the BASIN *output-contract*
level. For the requested `centered=true, extended=false` derivative, HARMONY must
internally pass both deepTools `--centerReads` and `--extendReads`: in deepTools,
`--extendReads` supplies paired-end fragment geometry, while `--centerReads` then
returns a read-length footprint centered on that fragment rather than retaining a
fragment-spanning footprint. By contrast, `centered=false, extended=true` passes
`--extendReads` alone and yields fragment-spanning coverage. BASIN rejects
`centered=true, extended=true` as ambiguous rather than silently mislabeling it.
These derivative choices cannot change the BAR scale basis or CODA's analytical
inputs.

Each current HARMONY run writes `run_manifest/coverage_contract.tsv`, and CODA
validates that contract before QC. Legacy HARMONY runs without a contract remain
readable for backward compatibility, but CODA prints a warning requiring manual
confirmation of the canonical geometry.

---

## Typical workflow

```
Your IgG BAMs + BAR file
        |
        v
  HARMONY (Pipeline 2)        <-- most users start here
  harmony.py
        |
        v
  HARMONY run directory
  ├── bigwigs_baseline/       canonical NONE + CPM
  ├── bigwigs_scaled/         canonical BASIN tracks (CODA input)
  ├── bigwigs_downstream/     optional centered/extended derivatives
  ├── scales/
  │   └── retained_bars.tsv   exact BARs retained per sample/drop
  ├── run_manifest/
  │   └── coverage_contract.tsv
  └── ...
        |
        v
  CODA (Pipeline 3)           <-- run from inside the HARMONY run directory
  coda.sh
        |
        v
  qc_[label]/                 <-- created inside the HARMONY run directory
  ├── experiment_qc_summary.tsv
  └── sample_qc_summary.tsv
```

If you need to build BAR files for a new organism or experimental system, run
OVERTURE first (see `OVERTURE/OVERTURE.md`). For dm6, use the distributed BAR
files in `HARMONY/templates/anchor_files/dm6/`.

---

## Distributed BAR files (dm6)

Three BAR file variants for *Drosophila melanogaster* (dm6) are in
`HARMONY/templates/anchor_files/dm6/`. All three were derived from 28 IgG
datasets spanning 8 independent laboratories, 6 tissues/cell types, and both
CUT&RUN and CUT&TAG. See `OVERTURE/OVERTURE.md` for guidance on which to use.

For most experiments, `...chr.bed` (the full panel) is recommended. It provides
the largest and most stable set of anchor regions.

---

## Quick start (HARMONY)

From the repository root, use the following steps. The launch script
expects to run from the `HARMONY/` directory; `config.yaml` resolves its input
paths relative to its own location.

```bash
cd HARMONY
mkdir -p logs
cp templates/config.yaml config.yaml
cp templates/samples.tsv samples.tsv

# Edit samples.tsv: replace example rows with your sample IDs, absolute BAM
# paths, antibody, genotype, sex, and replicate metadata.
# Edit config.yaml: set run_id_prefix and change universal_bed to:
#   templates/anchor_files/dm6/IgG_consensus_fdr0p005_relaxed_90pct_filter_away_lt50bp.chr.bed
# The copied config already has samples_tsv: samples.tsv.
# All other template settings can remain unchanged for a basic run.

sbatch slurm/run_harmony.sh config.yaml
```

The `universal_bed` path above is relative to the copied `config.yaml` in
`HARMONY/`. When using `templates/config.yaml` in place instead, the shorter
`anchor_files/dm6/...chr.bed` path in the template is correct. The example
manifest intentionally contains placeholder BAM paths, not bundled BAM data.

HARMONY will handle duplicate marking, canonical NONE/CPM bigWig generation,
BAR signal summarization, scale factor computation (at drop fractions 0.10,
0.20, 0.35, and 0.50 by default), exact retained-BAR recording, canonical scaled
bigWig generation, and any requested downstream derivatives automatically. All
outputs land in a timestamped run directory with full provenance logging.

---

## Repository structure

```
BASIN/
├── README.md                        this file
|
├── OVERTURE/                        Pipeline 1 — BAR discovery
|   ├── OVERTURE.md
|   ├── run_overture.sh              master wrapper
|   ├── step1_align/
|   ├── step2_bam_to_bedgraph/
|   ├── step3_seacr_grid/
|   ├── step4_panel_evaluation/
|   └── igg_samples_used/
|
├── HARMONY/                         Pipeline 2 — normalization
|   ├── HARMONY.md
|   ├── harmony.py
|   ├── slurm/
|   |   └── run_harmony.sh
|   └── templates/
|       ├── config.yaml
|       ├── samples.tsv
|       └── anchor_files/
|           └── dm6/                 distributed dm6 BAR files (start here)
|
└── CODA/                            Pipeline 3 — post-normalization QC
|   ├── CODA.md
|   ├── coda.sh

```

---

## Dependencies

| Tool | Role |
|------|------|
| Bowtie2 ≥2.4, cutadapt ≥4, samtools ≥1.15 | OVERTURE Step 1 (alignment) |
| deepTools ≥3 (bamCoverage) | OVERTURE Step 2, HARMONY |
| SEACR 1.3, bedtools ≥2.29, python3 ≥3.8 | OVERTURE Step 3 |
| bigWigAverageOverBed (UCSC tools), pandas ≥1.3 | OVERTURE Step 4 |
| samtools ≥1.17, deepTools ≥3, Python ≥3.10 (tested with 3.12.4) | HARMONY |
| R, tidyverse, ggplot2, deepTools, bedtools | CODA |

The provided Slurm launch scripts contain site-specific `module load` commands.
Adjust module names and versions for your cluster, or adapt the launchers to
use tools already on `PATH`; automatic fallback is not guaranteed.

---

## Citation

If you use BASIN, please cite:

Salzler et al. (2026) — preprint details to be added upon posting.

---

## Contact

Harmony Salzler, MateraLab, UNC Chapel Hill.
Questions and issues can be submitted via the repository issue tracker.
