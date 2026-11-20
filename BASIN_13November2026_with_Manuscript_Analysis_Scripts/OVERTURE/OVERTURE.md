# OVERTURE — BAR Discovery Pipeline

**OVERTURE** is Pipeline 1 of the BASIN suite. It discovers **BASIN Anchor Regions (BARs)** —
genomic regions of consistent nonspecific IgG enrichment that serve as normalization anchors
for HARMONY (Pipeline 2).

BARs are derived from IgG CUT&RUN or CUT&TAG data by calling peaks across a parameter grid,
filtering for consensus across samples, and benchmarking candidate panels for replicate
stability and genotype-dependent background behavior.

---

## Do you need to run OVERTURE?

**Most Drosophila users do not.** Pre-built universal BAR files for *Drosophila melanogaster* (dm6) are
distributed with BASIN and are ready to use directly in HARMONY. See
[`HARMONY/templates/anchor_files/dm6/`](../HARMONY/templates/anchor_files/dm6/).

Run OVERTURE only if you need to:

- Build BARs for a different organism
- Build an experiment-specific panel from your own IgG controls (e.g., a highly specialized
  system that diverges substantially from the universal dm6 compendium)
- Reproduce the universal dm6 BAR discovery from the original 28-dataset compendium

---

## Two discovery modes

OVERTURE supports two modes that use identical scripts. The only difference is which IgG
BAM files go into your manifest.

**Universal panel** uses a large, diverse set of public IgG datasets spanning many labs,
tissues, and assay types. This produces a robust, broadly applicable reference panel. The
distributed dm6 universal BARs were built from 28 datasets across 8 independent labs,
6 tissues, and both CUT&RUN and CUT&TAG (see `igg_samples_used/igg_samples_universal.tsv`
and Figure 3B of the BASIN paper).

**Experiment-specific panel** uses only the IgG controls from your own experiment, run across
multiple genotypes or conditions. This is appropriate when your system is highly unusual, or
when working with a new organism where no universal panel exists. For the experiment-specific
dm6 panels in Figure 3, IgG controls were collected from WT, Set2Δ, WT12xControl, and K36R
genotypes (3 replicates each where available).

Both modes go through the same four steps: align → bedGraphs → SEACR grid → panel evaluation.

---

## Quick start

The master wrapper script handles the entire job chain from a single config block.

```bash
# 1. Edit the CONFIG section at the top of run_overture.sh
#    The first thing to set is SCRIPT_DIR — the absolute path to your OVERTURE folder.
#    Slurm resolves paths from its spool directory, so relative paths will not work.
nano run_overture.sh

# 2. Dry run first — prints the job chain without submitting anything
sbatch run_overture.sh

# 3. When satisfied, submit
sbatch run_overture.sh --submit
```

The wrapper runs as a minimal Slurm job (1 CPU, 1G RAM, 10 minutes) that does
nothing except submit the four-step chain and exit. Each step starts only after
the previous one succeeds. The optional Step 4b medBAR genotype-shift typicality
diagnostic is run separately after Step 4 and does not require rerunning the
OVERTURE workflow.

---

## Workflow overview

```
IgG FASTQs (optional)        IgG BAMs (pre-aligned)
        |                            |
        v                            |
  Step 1: Align                      |
  align_igg_bowtie2.sh  _____________|
        |
        v
  Step 2: BAM to CPM bedGraphs
  run_overture_manifest.sh
        |
        v
  Step 3: SEACR grid + consensus
  run_seacr_grid.sh
  |-- SEACR_IgG_peaks.sbatch    (one job per FDR threshold)
  |-- SEACR_consensus.sbatch    (one job per FDR threshold, after SEACR)
        |
        v
  Candidate BED files in consensus_grid_collected/
        |
        v
  Step 4: Panel evaluation
  eval_panels.sbatch
        |
        v
  Optional Step 4b: medBAR genotype-shift typicality QC
  summarize_genotype_shift_typicality_BASE_R.R
        |
        v
  Select final BAR panel, use in HARMONY
```

---

## Step 1 — IgG alignment

**Script:** `step1_align/align_igg_bowtie2.sh`

This is the retained stringent Cutadapt + Bowtie2 end-to-end alignment workflow used to
prepare realigned IgG datasets for universal BAR discovery. It processes paired-end
FASTQs, applies adapter and quality trimming to both CUT&RUN and CUT&Tag samples,
and creates sorted and indexed full-length BAMs. It **does not** split alignments by
fragment length or use assay-specific alignment settings.

- **Cutadapt:** TruSeq and Nextera paired-end adapters, plus mosaic-end sequences;
  `-q 20 -m 25 -O 7 -e 0.10 --no-indels` by default.
- **Bowtie2:** `--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700`.
- **Reference:** dm6 Bowtie2 index (the study used FlyBase release 6.38).
- **Outputs:** `${PROJ}/alignment/bam/<sample>.sorted.bam` and BAM indexes; alignment
  summaries and Cutadapt logs are also produced.
- **CLI arguments:** `--proj`, `--fastq_dir`, `--bt2_index_prefix`; optional
  `--cores`, `--minlen`, `--qcut`, `--adapters_file`.

Run `bash step1_align/align_igg_bowtie2.sh --help` for its supported options.

The three Henikoff-laboratory L1 IgG datasets were analyzed from supplied aligned BAMs,
not realigned using this script.

In `run_overture.sh`, `RUN_ALIGNMENT=false` skips this optional step. When
`RUN_ALIGNMENT=true`, the script writes BAMs to `OUTROOT/alignment/bam/`.
**The manifest is not generated or rewritten automatically:** before submitting
Steps 2–4, make sure the BAM paths in `MANIFEST` point to the correct aligned
BAMs, including the supplied L1 BAMs where applicable. The downstream BAR-discovery
MAPQ filter is configured separately via `MIN_MAPQ` in the wrapper.

---

## Coverage geometry contract

OVERTURE defines BAR coordinates from an **uncentered, unextended, 1-bp aligned
read coverage representation**. HARMONY uses the same geometry when measuring
BAR signal and CODA uses it for baseline QC. This spatial geometry is fixed.
MAPQ is treated differently: it filters which reads contribute without changing
BAR coordinates, so it remains a recorded analysis parameter rather than a
hard-coded invariant.

---

## Step 2 — BAM to CPM bedGraphs

**Script:** `step2_bam_to_bedgraph/run_overture_manifest.sh`  
**Config:** `step2_bam_to_bedgraph/config_overture.env`

Step 2 generates the CPM bedGraphs that define the signal geometry used for BAR
discovery. To prevent a BAR panel from being discovered under one coverage
representation and quantified under another, OVERTURE now enforces a canonical
geometry:

- `--binSize 1`
- no `--centerReads`
- no `--extendReads`
- CPM normalization

These settings are not configurable. The old `BAMCOVERAGE_BINSIZE` and
`BAMCOVERAGE_EXTRA` escape hatches are rejected if they would change the
canonical geometry. `MIN_MAPQ` remains configurable (default `0`) and is written
to `report/coverage_contract.tsv` for provenance.

The manifest (`igg_samples.tsv`) contains `sample_id`, `layout` (PE or SE),
`bam` (absolute path to sorted/indexed BAM), and optionally `genome_id` and
`source`. See `step2_bam_to_bedgraph/templates/igg_samples.tsv`.

To run Step 2 independently:

```bash
sbatch step2_bam_to_bedgraph/run_overture_manifest.sh \
  my_igg_samples.tsv step2_bam_to_bedgraph/config_overture.env
```

Key config options are `GENOME_ID`, `RUN_TAG`, and `MIN_MAPQ`. MAPQ 0 matches the
distributed dm6 BAR discovery run. Higher MAPQ values may be used for new panel
discovery, but should be recorded and applied consistently when evaluating that
panel.

---

## Step 3 — SEACR grid and consensus

**Scripts:** `step3_seacr_grid/run_seacr_grid.sh`, `SEACR_IgG_peaks.sbatch`,
`SEACR_consensus.sbatch`

Runs SEACR peak calling across an FDR parameter grid and builds 80%, 90%, and 100% consensus
BED files for each FDR × stringency combination. This produces the library of candidate panels
that Step 4 evaluates.

SEACR is run in **no-control mode** (`non`) because the IgG controls are the samples
themselves — there is no paired IP track to use as a reference.

For each FDR threshold in `FDR_LIST` (default: `0.05 0.01 0.005 0.0005`), the grid submits a
SEACR job (relaxed and stringent), then a dependent consensus job that produces
`common_80pct.bed`, `common_90pct.bed`, and `common_100pct.bed` for each stringency, with and
without a 50 bp minimum-length filter. All collected candidate BEDs land in
`consensus_grid_collected/` with standardized names, for example:

```
IgG_consensus_fdr0p01_relaxed_90pct.bed
IgG_consensus_fdr0p005_stringent_80pct_filter_away_lt50bp.bed
```

The grid is intentionally broad — Figure 3C of the BASIN paper shows region counts across
this parameter space for both universal and experiment-specific runs. Step 4 identifies
candidate panels meeting evaluation criteria. These candidate panels are run 
outputs for evaluation and reproducibility; only the selected BAR panels intended 
for HARMONY use are distributed with BASIN.

To run Step 3 independently:

```bash
sbatch step3_seacr_grid/run_seacr_grid.sh /path/to/your/overture_run_dir
```

---

## Step 4 — Panel evaluation

**Script:** `step4_panel_evaluation/eval_panels.sbatch`

Benchmarks each candidate panel by extracting CPM bigWig signal over its regions using
`bigWigAverageOverBed` and computing **trimmed-median background summaries** across the
same default drop fractions used by HARMONY (`0.10, 0.20, 0.35, 0.50`). A custom
`DROP_LIST` can still be supplied for method-development sweeps. At the end of the
evaluation, a panel selection summary of passing panels is
automatically generated.

### Choosing validation/test datasets for panel evaluation

Step 4 requires the user to provide an independent set of CPM-normalized bigWigs
used to test candidate BAR panels. These evaluation datasets are not used to
discover BAR coordinates; they are used only to ask whether candidate panels
behave stably and representatively when applied to real experimental data.

A useful evaluation set should ideally include:

- biological replicates, so replicate stability can be measured;
- at least two genotypes or conditions, including a clearly defined reference
  genotype, so medBAR genotype shift can be calculated;
- multiple marks or assay targets when available, so panel behavior is not judged
  from a single chromatin feature;
- datasets generated and processed consistently across samples per test dataset;
- biologically relevant perturbations representative of the kinds of experiments
  in which the BAR panel is expected to be used.

The evaluation set does not need to be biologically "neutral." In fact, including
targets or genotypes with genuine global changes can be informative, because a
well-behaved panel should track the shared background response rather than behaving
as an isolated outlier. However, this also means that medBAR genotype shift and
genotype-shift typicality are properties of the evaluation context: changing the
targets, genotypes, or perturbations used for testing can change the observed shift
distributions and which panel appears most typical.

For a universal panel, the evaluation set should be independent of the IgG datasets
used for BAR discovery and should span more than one target if possible. For an
experiment-specific panel, evaluation can use the experiment's own signal datasets,
provided they contain biological replication and a reference condition.

Required inputs (set in `run_overture.sh` or exported before submitting manually):

- `BW_DIR` — a directory that must contain your `bw_metadata.tsv` file (see below).
  The CPM bigWig files themselves can live anywhere on disk — `bw_metadata.tsv` points
  to them by absolute path, so they do not need to physically be inside `BW_DIR`.
  These bigWigs should come from a validation experiment with biological 
  replicates across multiple genotypes, so that replicate stability and 
  genotype-dependent behavior of the candidate BAR panels can be assessed. 
  Genotype-dependent changes in medBAR are not assumed to represent panel 
  failure: global changes in specific signal can alter the fraction of reads 
  assigned to nonspecific background, producing the background titration that 
  BASIN is designed to measure.
  
For the universal dm6 panels, evaluation used a 6-target L1 CUT&TAG dataset 
(H3K36me2, H3K36me3, H3K27ac, pan-H4ac, RNAPII-S2P, and RNAPII-S5P). For 
experiment-specific panels, use bigWigs from your own dataset. The “most 
typical” universal panel referenced in the BASIN paper was then validated on 
additional CUT&RUN and CUT&TAG datasets, including wing-disc CUT&RUN (pan-H4ac, 
H3K9me3, and H3K27me3) and CUT&TAG datasets targeting H3K27me3, H3K4me1, and 
H3K36me3.

- `BW_DIR/bw_metadata.tsv` — a tab-separated file linking each bigWig to sample
  metadata **and its coverage-generation contract**. Required columns are
  `bw_path`, `sample_id`, `genotype`, `mark`, `rep`, `center_reads`,
  `extend_reads`, `bin_size`, `min_mapq`, and `normalization`. OVERTURE requires
  `center_reads=false`, `extend_reads=false`, `bin_size=1`, and
  `normalization=CPM` for every validation bigWig. `min_mapq` may differ from the
  distributed panel's MAPQ 0, but it must be the same across the validation set.
  This makes panel evaluation use the same spatial signal definition as BAR
  discovery and HARMONY's canonical scale basis. See the distributed template.

- `REGION_DIR` — the `consensus_grid_collected/` directory from Step 3.

The candidate panels include `.chr.bed` and `.nochr.bed` variants to accommodate
bigWigs with and without `chr`-prefixed chromosome names. Step 4 hashes candidate
BED files and evaluates byte-identical coordinate variants only once; variants that
still mismatch the bigWig chromosome convention produce all-zero signal and are
excluded with a logged warning.

**Panel selection thresholds** are set in the `run_overture.sh` CONFIG block and
passed through automatically. All thresholds can be adjusted for your organism and
experimental system:

- `EVAL_MIN_REGIONS` — minimum regions per panel (default: 200)
- `EVAL_MAX_REGIONS` — maximum regions per panel (default: 5000; raise for mammalian
  genomes, e.g. 10000)
- `EVAL_MAX_CV` — maximum mean replicate CV (default: 0.12)
- `EVAL_DROP_FRAC` — drop fraction used for panel evaluation (default:
  `median_drop20`)

### medBAR genotype shift

OVERTURE reports **medBAR genotype shift** as a diagnostic but does not use it as a
selection threshold or automated ranking criterion. For each mark and candidate
panel, genotype shift is the difference between the median trimmed medBAR for a
test genotype and the median trimmed medBAR for the configured reference genotype:

```text
medBAR genotype shift =
    median medBAR(test genotype) - median medBAR(reference genotype)
```

A negative value therefore means that BAR/background signal is lower in the test
genotype than in the reference genotype, whereas a positive value means that it is
higher. This is a shift in the **background medBAR distribution**, not a direct
measurement of the change in biological target signal.

The distinction matters because specific signal and nonspecific background are
compositionally related. If there is no true global difference in specific signal
between genotypes, a well-behaved BAR panel should generally show comparable
medBAR distributions between genotypes, apart from ordinary sampling and
experimental variation. By contrast, if a mark changes globally, the entire
background distribution can move even when the BARs themselves are not true target
sites. When a larger fraction of reads is devoted to specific signal, the relative
representation of nonspecific background can decrease; when global specific signal
decreases, nonspecific background can become relatively more represented. This
**background titration** is the behavior BASIN/HARMONY is intended to measure.

Accordingly, a coordinated medBAR shift observed across many independent passing
panels can be consistent with genuine dataset-wide biology rather than panel
failure. More concerning is a candidate panel whose shift is unusual or discordant
relative to the other passing panels. medBAR genotype shift should therefore be
interpreted across the panel ensemble and across marks rather than simply minimized
toward zero.

**What the trimmed-median measures:** OVERTURE now uses the same trimming
contract as HARMONY. For each candidate panel, rows that are zero across all
validation samples are removed first. For each sample and drop fraction,
`round(N × drop_fraction)` highest-signal BARs are removed and the median of the
remainder is calculated. The 0–70% evaluation range remains broader than the
usual HARMONY range so panel sensitivity to trimming is visible.

Output files (written to `OUTDIR`, default: `BW_DIR/panel_eval_out/`):

- `panel_trimmed_medians.tsv` — per-sample trimmed-median background for each
  panel/drop, including original panel size, panel size after the global all-zero
  filter, median region length, per-sample nonzero-region count, and validation
  MAPQ.
- `coverage_contract.tsv` — the geometry/MAPQ/trim contract verified for the
  validation bigWigs.
- `replicate_stability_trimmed.tsv` — coefficient of variation (CV) across replicates, per
  panel, mark, genotype, and drop fraction
- `genotype_shift_trimmed.tsv` — change in trimmed medBAR relative to the configured
  reference genotype, per panel, mark, and drop fraction. This is a diagnostic output
  and does not contribute to automated panel filtering or ranking.
- `panel_selection_summary.tsv` — summary of all panels passing the size and 
replicate-stability thresholds set in run_overture.sh. For each passing panel, 
the table reports panel size, median region length, mean replicate CV, and 
maximum absolute genotype shift at the selected evaluation drop fraction. Genotype 
shift is retained as a diagnostic for comparison among candidate panels but does 
not affect filtering. If no panels pass all filters, all panels are written so 
the thresholds can be reviewed manually. Only plain .bed panels are 
included — .chr.bed and .nochr.bed duplicates are excluded from the summary.

Example output TSVs from the original universal dm6 run are in
`step4_panel_evaluation/example_output/`.

To run Step 4 independently:

```bash
sbatch \
  --export=ALL,BW_DIR=/path/to/bigwigs,REGION_DIR=/path/to/consensus_grid_collected \
  step4_panel_evaluation/eval_panels.sbatch
```

When run independently, panel-selection thresholds default to the same values as
the master wrapper (`200`, `5000`, `0.12`, `median_drop20`) and can be overridden
with `EVAL_MIN_REGIONS`, `EVAL_MAX_REGIONS`, `EVAL_MAX_CV`, and `EVAL_DROP_FRAC`.
The requested `EVAL_DROP_FRAC` must actually exist in the evaluated drop set; the
script now fails rather than silently substituting another drop.

---

## Optional Step 4b — medBAR genotype-shift typicality QC

**Script:** `summarize_genotype_shift_typicality_BASE_R.R`

This optional post-evaluation diagnostic asks whether each **passing** panel follows
the central medBAR genotype-shift behavior of the passing-panel ensemble. It does
not change Step 4 pass/fail status and does not automatically select a panel.

The script is intentionally standalone so that it can be run after Step 4 without
rerunning OVERTURE. It uses base R only and requires no additional R packages.

From the `OVERTURE/` directory, the script can usually locate the most
recently modified valid panel-evaluation directory automatically:

```bash
Rscript summarize_genotype_shift_typicality_BASE_R.R
```

For reproducible analysis of a specific run, an evaluation directory can instead be
supplied explicitly:

```bash
Rscript summarize_genotype_shift_typicality_BASE_R.R \
  /path/to/overture_run/panel_eval_out
```

The script reads:

- `genotype_shift_trimmed.tsv`
- `panel_selection_summary.tsv`

and writes a `genotype_shift_typicality_qc/` directory containing:

- `panel_genotype_shift_typicality.tsv` — panel-level mean and median typicality
  distances and their ranks
- `genotype_shift_typicality_detail.tsv` — mark × genotype values used for the
  calculation
- `genotype_shift_passed_panels_signed.png` — signed medBAR genotype-shift
  distributions across passing panels, using a common y-axis across facets

### What "typicality" means

Typicality is calculated **within each mark × genotype comparison** so that marks
with intrinsically larger medBAR shifts do not dominate the summary. For each
comparison, the median signed shift across passing panels is treated as the ensemble
center and each panel's robust standardized distance is calculated as:

```text
typicality distance =
    | panel shift - ensemble median shift | / MAD
```

where MAD is the median absolute deviation across passing panels.

Two panel-level summaries are reported:

- **Mean typicality distance** — the mean standardized distance across all evaluated
  mark × genotype comparisons. Lower values indicate behavior that is more
  representative overall. Because this is a mean, occasional large deviations are
  penalized.
- **Median typicality distance** — the median standardized distance across
  comparisons. Lower values indicate behavior that is more representative in a
  typical comparison and are less affected by isolated outliers.

A low typicality distance does **not** mean that genotype shift itself is small.
For example, if a real global biological change causes the medBAR distribution to
shift for a particular mark, many passing panels may move together. A panel can
therefore show a large medBAR genotype shift while still having a low typicality
distance because it follows the ensemble.

Conversely, a panel can have only a modest absolute shift but a high typicality
distance if it behaves differently from the other passing panels.

Typicality is therefore a **secondary robustness diagnostic**, not a universal
property of a BAR panel and not an independent optimization target. Its value
depends on the validation datasets used for Step 4, including the marks, genotypes,
reference genotype, and biological perturbations represented. A different
evaluation set may change which panel appears most typical.

---

## Panel selection

After Step 4, review the panel-evaluation outputs together rather than treating any
single diagnostic as sufficient on its own.

**Panel size: ~300–5,000 regions.** Panels below 200 regions are automatically excluded.
Very large panels approaching the upper limit in smaller genomes may include marginal
peaks and warrant closer inspection.

**Replicate CV: ~0.07–0.12, lower is better.** Check `replicate_stability_trimmed.tsv`.
Filter to a moderate drop fraction (e.g., `median_drop20`) and look at the `cv` column across
panels. Panels with CV consistently above 0.15 across multiple marks should be deprioritized.
A panel whose CV stabilizes at low trim levels (10–20%) is preferable to one that requires
50%+ trimming to converge — the latter suggests contamination with real signal.

**Genotype-dependent medBAR behavior: diagnostic, not a minimization criterion.** Check
`genotype_shift_trimmed.tsv`. The medBAR genotype-shift value measures the change in
trimmed medBAR between each genotype and the configured reference genotype. It should not
be read as a direct measurement of biological target-signal change. When specific signal
changes globally, the relative representation of nonspecific background can change as
well, displacing the medBAR distribution through background titration. Genotype shifts
should therefore be interpreted in biological context and compared across candidate
panels and marks. Reproducible, broadly shared shifts may reflect dataset-level biology,
whereas unusually panel-specific or discordant behavior may indicate contamination by
genotype-responsive specific signal. The optional Step 4b typicality diagnostic provides
a robust way to distinguish representative ensemble behavior from panel-specific
deviations.

The goal is a balanced panel: large enough to provide a stable median estimate and
reproducible across biological replicates. Genotype-dependent medBAR behavior is retained
as a diagnostic rather than minimized, because genuine global changes in specific signal
can produce corresponding changes in background representation. As shown in Figure 3D,
relaxed SEACR mode produces larger panels and stringent mode produces smaller ones.

---

## The distributed dm6 universal BAR files

Three variants are distributed with BASIN, in `HARMONY/templates/anchor_files/dm6/`. All three
come from the same 28-dataset compendium and have been validated on independent CUT&RUN
and CUT&TAG datasets.

`...chr.bed` is the **recommended default** for most users. It includes all chromosomes
and histone gene regions, giving the largest and most stable panel. More anchor regions
mean more robust background estimation.

`...NO_HISTONE_GENE_PEAKS.bed` includes all chromosomes but has regions overlapping histone
gene clusters removed. Histone genes are highly transcribed and can have copy number
variation and mappability issues in some contexts. Use this if those factors are a concern
in your system.

`...AUTOSOMES_ONLY.bed` excludes sex chromosomes entirely. Use this when sex chromosome
composition is genuinely unknown or uncontrolled across samples — for example, a mixed-sex
population where chrX and chrY IgG signal might vary in ways that introduce normalization
bias. In most well-controlled experiments (including cell lines where sex is known), the
full chr.bed is preferable because the additional anchors on sex chromosomes improve
panel stability.

---

## The universal IgG dataset compendium

The 28 IgG datasets used to build the universal dm6 panels are documented in
`igg_samples_used/igg_samples_universal.tsv`. All are publicly available and GSM/SRR
accession numbers are preserved for reproducibility. Local cluster paths in the `bam` column
have been replaced with placeholders — update these with your local paths before running.

Dataset composition: 28 samples, 8 independent labs, 6 tissues/cell types (S2 cells, L1
larvae, ovary, larval brain, wing disc, embryo NC14, adult head), both CUT&RUN and CUT&TAG.

---

## File structure

```
OVERTURE/
├── OVERTURE.md                                 this file
├── run_overture.sh                             master wrapper (start here)
├── summarize_genotype_shift_typicality_BASE_R.R    optional Step 4b typicality QC
|
├── step1_align/
|   └── align_igg_bowtie2.sh                   trim + align FASTQs to BAMs
|
├── step2_bam_to_bedgraph/
|   ├── run_overture_manifest.sh               BAMs to CPM bedGraphs
|   ├── config_overture.env                    config template
|   └── templates/
|       └── igg_samples.tsv                    manifest template
|
├── step3_seacr_grid/
|   ├── run_seacr_grid.sh                      FDR grid driver
|   ├── SEACR_IgG_peaks.sbatch                 SEACR peak calling (relaxed + stringent)
|   └── SEACR_consensus.sbatch                 80/90/100% consensus BED builder
|
├── step4_panel_evaluation/
|   ├── eval_panels.sbatch                     trimmed-median panel benchmarking
|   ├── templates/
|   |   └── bw_metadata_template.tsv          bigWig metadata template
|   └── example_output/
|       ├── panel_trimmed_medians.tsv          example results from dm6 universal run
|       ├── replicate_stability_trimmed.tsv
|       └── genotype_shift_trimmed.tsv
└── igg_samples_used/
    └── igg_samples_universal.tsv              28-dataset compendium (accessions preserved)
```

---

## Dependencies

The included Slurm scripts are configured for an HPC module environment.
Adjust site-specific `module load` calls or adapt the launchers to use tools
already available on `PATH`; automatic fallback is not guaranteed.

- Step 1: Bowtie2 ≥2.4, cutadapt ≥4, samtools ≥1.15
- Step 2: samtools ≥1.15, deepTools ≥3 (bamCoverage)
- Step 3: SEACR 1.3, bedtools ≥2.29, python3 ≥3.8
- Step 4: bigWigAverageOverBed (UCSC tools), python3 ≥3.8, pandas ≥1.3
- Optional Step 4b: base R (`Rscript`); no additional R packages required

---

## Relation to the BASIN paper

OVERTURE corresponds to Figure 3 of Salzler et al. The BAR discovery workflow, parameter
grid, panel selection criteria, and validation against independent datasets are described in
the Methods under "BAR Discovery." The universal dm6 panels distributed with BASIN were
selected from this candidate grid and validated across independent CUT&RUN and CUT&TAG
datasets.
