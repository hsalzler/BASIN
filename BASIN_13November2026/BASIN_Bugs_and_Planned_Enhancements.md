# BASIN — Bugs and Planned Enhancements

## Resolved

- **Coverage geometry is now an enforced cross-pipeline contract (2026-09-23).**
  OVERTURE BAR discovery, HARMONY's raw NONE scale basis, HARMONY's canonical CPM
  track for CODA, and HARMONY's canonical BASIN-scaled output are fixed to 1-bp,
  uncentered, unextended coverage. The old free-form OVERTURE `BAMCOVERAGE_EXTRA`
  hook and HARMONY global `center_reads`/`extend_reads` switches can no longer
  alter normalization/QC geometry. HARMONY writes `coverage_contract.tsv`; CODA
  validates it. Centered and/or extended tracks are optional downstream
  derivatives in `bigwigs_downstream/` and are never read by CODA. MAPQ remains
  configurable, is applied consistently to canonical HARMONY products, and is
  recorded in the contract.

- **Exact retained BAR identities are now recorded and consumed.** HARMONY writes
  `scales/retained_bars.tsv` for every sample/drop after global all-zero filtering
  and trimming. CODA uses those identities for retained-BAR chromosome composition
  and the post-BASIN `median ≈ 1` integrity check. Legacy runs without the file use
  a corrected reconstruction path.

- **BAR trimming definitions were inconsistent across pipelines.** HARMONY used
  `round(N × drop)` after removing BAR rows that were zero across all samples,
  while CODA and OVERTURE evaluation used different rounding/zero-row rules. CODA
  and OVERTURE Step 4 now reproduce HARMONY's global all-zero filter plus
  `round(N × drop_fraction)` definition.

- **OVERTURE Step 4 could evaluate the wrong requested trim threshold.** The
  selection code tested whether `EVAL_DROP` was a dataframe column name instead
  of a value in the `drop` column, allowing a silent fallback to another drop. It
  now requires the requested drop to be present and fails clearly otherwise.

- **OVERTURE Step 4 accepted validation bigWigs without coverage provenance.**
  Validation metadata now records and verifies `center_reads`, `extend_reads`,
  `bin_size`, `min_mapq`, and normalization. Panel evaluation requires uncentered,
  unextended, 1-bp CPM validation tracks and one consistent MAPQ across the set.


- **Post-BASIN anchor convergence was measured on the untrimmed anchor set.** The
  check computed the scaled anchor median over the full ANCHOR_BED, including the
  top high-signal anchors that carry residual real signal — exactly the anchors
  the drop-fraction trim exists to exclude — so the metric was contaminated (and
  sat well above 1.0, e.g. ~1.32). It now summarizes over the trimmed kept set at
  the primary drop (top `primary_drop` fraction removed, reproducing HARMONY's
  kept set), so `primary_anchor_trimmed` and the group median/CV/z describe clean
  background and center on ~1.0 by construction. The untrimmed full-set median is
  still reported as `primary_anchor` for reference.

- **CODA baseline units were incorrectly treated as display-configurable.** An
  earlier patch allowed CPM or RPGC baseline bigWigs interchangeably. That is no
  longer permitted because CODA contains absolute baseline thresholds whose units
  must be stable. HARMONY now always produces a canonical 1-bp, uncentered,
  unextended **CPM** QC track, and CODA consumes CPM only. Optional RPGC or other
  user-facing representations belong in downstream derivatives, not CODA input.

- **CODA flagged "missing primary drop" for every sample (scale-table join-key
  regression).** When HARMONY's scale basis moved from the CPM to the raw (NONE)
  bigWigs, the columns in `scales_long.tsv` became `{sample_id}.NONE.bw`, but
  CODA's `normalize_sample_name` only stripped `.CPM.bw` / `.scaled…dedup.bw`
  suffixes, not `.NONE.bw`. The scale table therefore failed to join to the
  scaled bigWigs (`primary_scale_factor` NA for all samples). Added `.NONE.bw`
  stripping to `normalize_sample_name` so all three name forms reduce to the same
  sample id.

- **CODA died silently on a chromosome-naming mismatch.** STEP 1 filtered the
  domain/anchor/genome-sizes files with a bare `grep '^chr...'`. Under
  `set -euo pipefail`, a zero match (e.g. a non-`chr`-prefixed `chrom.sizes` or
  peak file) returned non-zero and aborted the whole job with no message in
  `.out` and an empty `.err`. Replaced with a `normalize_dm6_chr` helper that
  accepts both `chr2L` and bare `2L` naming (emitting `chr`-prefixed output to
  match the bigWigs) and, if nothing matches, exits with a clear error naming the
  offending file and showing its first lines. Also added an `ERR` trap so any
  future `set -e` failure reports its line number instead of vanishing.

- **HARMONY normalization was not depth-invariant (critical).** Scale factors
  were computed in CPM space and then applied to coverage via
  `bigwigCompare --operation first --scaleFactors sf:1`, which double-corrected
  for library size (a residual depth term survived, so more deeply sequenced
  samples read artifactually higher) and silently rebinned the output to
  `bigwigCompare`'s 50 bp default. HARMONY now derives the scale factor from raw
  (NONE) BAR coverage (`sf = 1/median(raw_BAR after trimming)`) and applies it on
  the BAM via `bamCoverage --normalizeUsing None --scaleFactor sf`, so the depth
  term cancels exactly and the output is at the configured bin size.
  `bigwigCompare` is no longer used. A canonical CPM bigWig is produced as a
  fixed-geometry analytical QC input for CODA; optional display derivatives are
  separate.

- **Prepared-BAM reuse across runs.** New `bams_prepared` config key reuses
  already markdup/dedup'd BAMs from any directory (e.g. a prior run's
  `bams_prepared/`), skipping the expensive duplicate-handling step.

- **CODA output clutter.** All secondary tables and figures now go to a
  `qc_[label]/secondary_outputs/` subdirectory; only `sample_qc_summary.tsv` and
  `experiment_qc_summary.tsv` remain at the top of `qc_[label]/`.

- **Post-BASIN anchor pipeline-integrity check added to CODA.** Measured on the
  canonical scaled bigWigs: per-sample retained-BAR medians and per-group CV,
  with sample- and experiment-level flags. Because the retained median is expected
  to be ~1 by construction, this is a mechanical reconciliation check rather than
  an independent biological validation metric
  (`post_basin_anchor_check.tsv`). Thresholds `POST_BASIN_ANCHOR_Z` (2.0) and
  `POST_BASIN_ANCHOR_CV_WARN` (0.15).

- **Mismatched default anchor file between HARMONY and CODA.** HARMONY's config
  defaulted to `...NO_HISTONE_GENE_PEAKS.bed` while CODA defaulted to
  `..._AUTOSOMES_ONLY.bed`, so a user accepting both defaults would compute scale
  factors and anchor QC over different BAR sets. CODA's default now matches
  HARMONY (`...NO_HISTONE_GENE_PEAKS.bed`), and the CODA config comment warns that
  the file must match HARMONY's `universal_bed`. This describes the 2026-09-23
  patch; the 2026-11-13 release instead recommends the full `...chr.bed` panel
  as the default in both pipelines. CODA reads raw medBAR from
  `scales_long.tsv` to rank drops by cross-replicate reproducibility of the
  background estimate (non-circular) plus a scale-factor elbow signal.
  Outputs `drop_stability_by_experiment.tsv`, `drop_stability_by_genotype.tsv`,
  and `drop_recommendation.tsv`. Partially addresses enhancement "anchor coverage
  QC / drop selection" below.


- **OVERTURE panel-selection duplicate coordinate variants.** Step 3 intentionally
  emits native, `.chr`, and `.nochr` BED variants. Step 4 now hashes candidate
  BED contents and evaluates byte-identical copies only once, preventing the
  native file and its matching explicit coordinate variant from appearing as
  duplicate biological panels in `panel_selection_summary.tsv`.

- **OVERTURE's default drop sweep was broader than the HARMONY production set.**
  Step 4 now defaults to `[0.10, 0.20, 0.35, 0.50]`, matching HARMONY. A custom
  `DROP_LIST` can still be supplied explicitly for method-development studies.

- **OVERTURE Step 4 left `tmp_aob.tab` behind.** The temporary
  `bigWigAverageOverBed` table is now registered for cleanup on normal Python
  exit.

- **CODA could duplicate the `replicate` column during summary joins.** The
  manifest's replicate value is now authoritative; right-hand join tables drop
  their duplicate `replicate` field before merging, preventing `.x`/`.y`
  variants in `sample_qc_summary.tsv`.

## Remaining bugs

No additional correctness bugs are currently tracked in this file as of the 
2026-09-23 patch. The patch was statically/regression tested as described in 
`BASIN_VALIDATION_2026-09-23.txt`.

## Enhancements

1. Parallelize `bigWigAverageOverBed` calls in `eval_panels.sbatch` (OVERTURE) —
   current serial implementation takes ~2–6 h; parallelization would reduce to
   ~10–20 min, important for mammalian genomes.
2. Generate only the relevant `.chr` or `.nochr` BED variant based on input BAM
   chromosome convention rather than all three variants — cuts panel count by
   two-thirds and speeds up evaluation.
3. Add mammalian BAR panels (hg38, mm10) using the Nordin et al. 2023 dataset
   compendium.
4. Add an SRA download and alignment wrapper for ingesting public IgG datasets.
5. Validate and distribute human and mouse anchor panels once mammalian data are
   processed.
6. Add an `experiment` column to the HARMONY manifest so it propagates into CODA,
   rather than setting experiment groups in the CODA config block.
7. Extend post-BASIN anchor-convergence QC from an experiment-wide warning to
   genotype-group-specific convergence metrics and flags, with groups containing
   fewer than three samples reported as not formally assessed.
8. Add a CODA anchor-convergence plot showing per-sample trimmed anchor signal at
   the primary BASIN threshold, with genotype-group convergence warnings
   annotated at the group level.
9. Incorporate genotype-shift concordance into OVERTURE panel evaluation. For each
   mark and biologically matched genotype comparison, calculate the panel-specific
   shift as log2(medBAR_genotype / medBAR_reference), estimate the consensus shift
   across candidate panels using the median, and quantify each panel's residual from
   that consensus. Panels with small residuals across marks and comparisons should be
   favored as providing scaling factors representative of the broader background
   response, without penalizing genuine genotype-dependent background titration.
   Initially report this as a diagnostic and empirically calibrate any weighting or
   threshold before incorporating it into automated panel ranking.
10. Add a script or some form of automated peak file generation for CODA inputs
   to streamline the transition from HARMONY to CODA.
11. Extend CODA region-set and chromosome handling beyond Drosophila dm6,
    with organism-specific reference files and validation before declaring
    support for additional genomes.
