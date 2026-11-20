# CODA: Post-Normalization Quality Control
## BASIN Pipeline 3 — Genotype-Aware Methods Description and Decision Framework

---

## CODA Revision

CODA retains all sample-level metrics, thresholds, and QC classifications from
CODA, but refines one experiment-level failed-QC rule. In the original decision
framework, an experiment was assigned Failed QC whenever at least one-third of
samples were flagged, regardless of how those flags were distributed. CODA2
distinguishes broad experiment-level instability from a coherent pattern
confined to a single genotype group:

- If at least one-third of samples are flagged **across two or more genotype
  groups**, the experiment receives **Failed QC**.
- If at least one-third of samples are flagged but the flags are **confined to
  one genotype group**, the experiment receives **Manual review**
  rather than Failed QC.

This revision does not suppress or remove any sample-level flag. It changes only
the experiment-level interpretation of flag distribution and adds explicit
summary columns documenting how many genotype groups contain flagged samples
and which groups they are.

## Overview

CODA is a post-normalization quality control framework that evaluates whether
BASIN normalization succeeded, identifies individual samples that may be
problematic, and issues experiment-level recommendations. It is run after
Pipeline 2 has produced scaled bigWig files, and it operates on two distinct
sets of outputs from that run: the baseline CPM-normalized bigWigs, and the
BASIN-scaled bigWigs at each drop threshold. These two sets of files serve
different and complementary roles in the QC framework. Most diagnostic metrics —
including interpeak diversion, genotype consistency, anchor support, and
baseline weirdness — are computed exclusively on the baseline CPM bigWigs,
because these reflect the pre-normalization signal landscape and are therefore
most informative for determining whether BASIN's assumptions are met. A second
set of metrics is computed on the scaled bigWigs to evaluate whether
normalization produced the expected outcome: specifically, whether samples are
internally consistent after scaling and whether any sample is an outlier within
its experiment group. The final QC recommendations integrate signals from both.

The core principle underlying CODA is that successful BASIN normalization
produces a predictable relationship between three classes of genomic regions:
domains (regions of high, biologically meaningful signal for the antibody being
used), anchor regions (BASIN anchor regions, or BARs — the IgG-derived
nonspecific enrichment sites that BASIN uses for normalization), and interpeak
regions (genomic intervals that fall outside domains, anchors, and repetitive
elements, and thus represent background). In a well-normalized experiment,
domain signal is high relative to background, anchor signal is modestly elevated
relative to background (reflecting the nonspecific enrichment that makes BARs
useful), and interpeak signal is low and consistent across samples. Deviations
from these expected relationships — whether at the level of individual samples
or entire experiment groups — are what CODA is designed to detect.

CODA is implemented as a SLURM-submittable bash script (`coda.sh`) that
generates and executes an R script using deepTools, bedtools, the tidyverse,
and ggplot2. All QC thresholds are configurable via environment variables,
with sensible defaults documented below. **This release of CODA supports only
Drosophila dm6 primary chromosomes** (chr2L, chr2R, chr3L, chr3R, chr4,
chrX, and chrY); OVERTURE and HARMONY can use organism-specific BAR panels,
but CODA has not been generalized to other genomes.

---

## Input Requirements

CODA requires the following outputs from a completed Pipeline 2 run:

- `bigwigs_baseline/`: canonical CPM bigWig files, one per sample
- `bigwigs_scaled/drop*/`: BASIN-scaled bigWig files at each drop threshold
- `scales/scales_long.tsv`: The table of per-sample, per-threshold scale factors
- `run_manifest/samples.resolved.tsv`: The resolved sample manifest from the
  Pipeline 2 run
- `signal_summary/universal_peaks.rawCounts.tsv`: The multiBigwigSummary output
  over BARs from the baseline bigWigs

In addition, CODA requires four reference paths for dm6:

- `DOMAIN_BED`: A BED file of regions of expected signal enrichment for the
  antibody being assessed — broad domains for diffuse chromatin marks, or a peak
  file for sharper targets such as transcription factors. This file is used to
  evaluate whether target signal is behaving as expected relative to background.
  It does not need to be perfect — it serves as a reference signal class, not a
  ground truth. No example domain BED is distributed. Supply a BED appropriate
  to the antibody and experiment (including H3K27me3 experiments).
  
  IMPORTANT: Peak/domain geometry requirement: If the supplied peak or domain 
  BED was generated empirically from coverage tracks, those tracks should use 
  the same canonical coverage geometry as CODA’s analytical 
  inputs: --centerReads disabled and --extendReads disabled. This avoids 
  mismatches between region boundaries and the signal representation used for CODA 
  domain/interpeak calculations.
  
- `ANCHOR_BED`: The BAR file used in Pipeline 2 (the universal Drosophila file
  is provided; see Pipeline 1 for generating organism-specific files).
- `GENOME_SIZES` and `REPEAT_GTF`: Chromosome size file and a
  repeat/transposable element annotation, used to define interpeak regions by
  exclusion.

---

## Quick start (Slurm; dm6 only)

Run CODA after HARMONY has finished, from the BASIN repository root. Provide
absolute paths to all four required reference files. In particular,
`ANCHOR_BED` **must be identical to the BAR file used by HARMONY** (see its
`run_manifest/config.resolved.yaml`). `DOMAIN_BED` should describe enrichment
regions for the antibody being evaluated. These files are not all bundled.

```bash
export DOMAIN_BED=/absolute/path/to/target_domains_or_peaks.bed
export ANCHOR_BED=/absolute/path/to/HARMONY/templates/anchor_files/dm6/IgG_consensus_fdr0p005_relaxed_90pct_filter_away_lt50bp.chr.bed
export REPEAT_GTF=/absolute/path/to/dm6_BDGP_rmsk_TE.gtf.gz
export GENOME_SIZES=/absolute/path/to/dm6.chrom.sizes

sbatch --export=ALL CODA/coda.sh /absolute/path/to/completed/HARMONY/run_directory
```

Replace every `/absolute/path/to/...` placeholder before submitting. The
example `ANCHOR_BED` is the recommended full dm6 BAR panel; if HARMONY used a
different panel, point CODA to that **same** panel. The provided Slurm script
loads site-specific environment modules; modify those lines if module names
or versions differ on your system. Results are written to
`qc_postrun_qc/` under the HARMONY run directory by default. You may set
`RUN_LABEL` to choose a different suffix.

---

## Region Classes

CODA computes signal statistics over three non-overlapping classes of genomic
regions for each sample:

**Target regions (referred to throughout as "domains").** Regions of expected
signal enrichment for the antibody being assessed. These are provided by the
user and are experiment-specific. For diffuse chromatin marks, broad domain
files are appropriate; for sharper targets such as transcription factors or
histone modifications with punctate enrichment, a peak file is suitable. The
pipeline is agnostic to the breadth of these regions — what matters is that they
capture where signal is expected to be high. No example domain BED is
distributed with CODA; users supply their own appropriate file.

**Anchors (BARs).** The same BASIN anchor regions used for normalization in
Pipeline 2. In a successful normalization, anchor signal should be consistent
across samples after scaling.

**Interpeak regions.** Genomic intervals that remain after excluding domains,
anchors, and repetitive elements, subject to a minimum size filter (default: 100
bp). These represent background — regions where no biologically meaningful
signal is expected. Interpeak signal serves as the key reference class for
evaluating normalization quality, because it should be low and consistent across
samples if normalization is working correctly.

Signal statistics (median, mean, fraction nonzero) are computed with
`multiBigwigSummary` over each region class for both canonical CPM bigWigs and
canonical BASIN-scaled bigWigs. Because both are fixed-geometry analytical
products, interpeak-diversion and ratio metrics cannot change merely because a
user requested a centered visualization track.

---

## Experiment Grouping

All QC metrics in CODA — including interpeak diversion z-scores, genotype
consistency, baseline weirdness, and experiment-level recommendations — are
computed within experiment groups rather than across the entire sample sheet.
This is a critical design feature, because samples from different experimental
contexts should never be evaluated against each other. Absolute signal levels,
background characteristics, and enrichment ratios differ systematically between
experiments run at different times, between different antibodies, and between
different fractions of the same assay (e.g., supernatant versus pellet in
CUT&RUN). Computing z-scores or consistency metrics across such groups would
produce meaningless results — samples from one experiment would appear as
outliers simply because of differences in experimental context rather than
sample quality.

CODA supports explicit user-defined grouping via the `MANUAL_QC_GROUPS`
environment variable. Groups are specified as a semicolon-separated list of
named groups, where each group contains a colon-separated group name followed by
a comma-separated list of sample IDs matching those in the sample manifest
exactly. For example:

```
MANUAL_QC_GROUPS="experiment1:sampleA,sampleB,sampleC;experiment2:sampleD,sampleE,sampleF"
```

Every sample in the manifest must be assigned to exactly one group. If samples
are missing from manual QC groups when specified, they will be ignored. 

If `MANUAL_QC_GROUPS` is left empty, all samples are placed in a single group
called `all_samples`. This is appropriate only when the entire sample sheet
represents a single coherent experiment — that is, samples that were processed
together and whose signal characteristics are expected to be directly
comparable.

Common situations that require explicit grouping include:

- Multiple experiments run months apart and processed in the same Pipeline 2
  run, which may have systematically different signal levels or background
  characteristics
- Supernatant and pellet fractions from the same CUT&RUN experiment, which have
  fundamentally different enrichment profiles and should be evaluated
  independently
- Different antibodies run together on the same sample sheet, which will have
  different domain/interpeak ratios by nature

When in doubt, it is better to define explicit groups than to rely on the
default single-group behavior. A group that is too broadly defined will produce
z-scores and consistency metrics that are difficult to interpret, and may mask
genuine sample-level problems or generate spurious flags.

---

## Replicate Count Requirements

CODA's statistical framework — including z-score-based outlier metrics,
MAD-based interpeak-diversion metrics, and genotype-consistency metrics — is
designed for experiments with at least three replicates per genotype group. When
this assumption is not met, CODA suppresses the affected metrics rather than
reporting values that appear meaningful but are not, and issues explicit
sample-level warnings and action-oriented recommendations in the experiment
summary.

**n = 1 (single sample per genotype group).** Within-group outlier statistics
are undefined or uninformative with a single sample, and dispersion-based
genotype-consistency metrics cannot be interpreted. The affected z-score and
consistency metrics are therefore suppressed and set to NA for that sample. The
sample retains its warning classification and the experiment receives a Pass
with caveat recommendation with a plain-English note specifying which genotype
group was affected. Review the sample QC summary to assess the affected sample
directly; CODA's statistical framework requires at least 3 replicates per group
for these comparative metrics to be meaningful.

**n = 2 (two samples per genotype group).** With only two samples, standardized
outlier scores are constrained by the two-point geometry and do not provide a
meaningful assessment of whether either sample is unusually different.
Dispersion estimates are similarly unstable. The affected z-score and
consistency metrics are therefore suppressed for genotype groups with n=2, and
the experiment receives a Pass with caveat recommendation noting that results
should be interpreted with caution.

**n ≥ 3 (three or more samples per genotype group).** Normal operation. All
metrics are computed and interpreted as described in the sections below.

The replicate count is evaluated per genotype group within each QC group, not
across the experiment as a whole. An experiment with three genotypes of three
replicates each will operate normally; the same experiment with one genotype
having only two replicates will suppress z-score metrics for that genotype group
and issue Pass with caveat at the experiment level, while the other genotype groups are evaluated normally.

Note that anchor support and dropout are not affected by replicate count — these
are computed per sample without reference to group-level statistics and will
always be reported regardless of n.

---

## QC Metrics

The following metrics are computed from the baseline CPM bigWigs. They capture
the pre-normalization signal landscape and are used to assess whether BASIN's
assumptions are met and whether individual samples are behaving as expected
before scaling is applied.

**Signal ratios.** Three ratios capture the relationship between region classes
in the baseline CPM bigWigs:

- *Domain/interpeak*: the ratio of median target region signal to median
  interpeak signal. This should be high in a successful experiment — enriched
  regions are elevated relative to background.
- *Anchor/interpeak*: the ratio of mean anchor signal to mean interpeak signal.
  Anchors should be modestly elevated above background.
- *Anchor/domain*: the ratio of anchor signal to target region signal. In a
  successful normalization, this should be relatively low — normalization
  anchors are not the primary signal.

**Anchor support.** The fraction of BAR regions with nonzero signal in the
baseline CPM bigWig for a given sample. Low anchor support (below 15% nonzero,
or fewer than 50 nonzero regions) indicates that BARs are not being sampled
adequately in the pre-normalization data, which compromises the reliability of
the scaling step.

**Baseline weirdness.** A sample is flagged as "baseline weird" if its median
interpeak signal or its mean anchor/interpeak ratio in the baseline CPM bigWig
has an ordinary within-experiment z-score greater than 2.5 in absolute value.
This identifies samples whose pre-normalization signal profile is anomalous
relative to the rest of the experiment. This metric is suppressed and set to NA
for samples in genotype groups with fewer than 3 replicates — see Replicate
Count Requirements.

**Baseline domain/interpeak outlier.** A sample is flagged as a baseline
domain/interpeak outlier if its baseline domain/interpeak ratio is unusually
high or low relative to the experiment group. This uses the same MAD-based
robust domain/interpeak z-score computed for interpeak diversion, but evaluates
both directions rather than only the low direction. This flag captures samples
whose pre-normalization signal-to-background relationship is unusual even when
anchor support is adequate and no dropout occurs.

The following metrics are computed from the BASIN-scaled bigWigs and assess
whether normalization produced the expected outcome.

**Post-normalization group outlier.** At the primary drop threshold, the median
domain signal and median interpeak signal of each sample's scaled bigWig are
converted to ordinary z-scores within the experiment group. A sample is flagged
as a post-normalization outlier if either z-score exceeds 2.0 in absolute value.
This detects samples that remain anomalous after scaling, which may indicate that
the scale factor was unreliable for that sample.

**Sex chromosome anchor skew.** The chromosome composition of the retained BAR
regions — those not trimmed during scale factor computation — is evaluated.
For current HARMONY runs, CODA reads the exact retained BAR identities from
`scales/retained_bars.tsv` rather than reconstructing them from scaled signal.
Legacy runs use the corrected reconstruction rule. The composition is evaluated
across drop thresholds. If a sample shows an unusual
enrichment of sex-chromosomal anchors relative to the experiment group (z-score
> 2.0), this is flagged as a warning. Sex chromosome skew can indicate that the
anchor trimming step is not behaving symmetrically across the genome, which may
affect the reliability of the scale factor.

**Dropout.** A sample is considered to have dropped out at a given threshold if
no valid scale factor could be computed — which occurs when too few BAR regions
remain after trimming and all-zero filtering to produce a reliable median.
Dropout is determined from the scales table rather than from the bigWigs
directly, and it is evaluated across all drop thresholds run in Pipeline 2.

---

## Interpeak Diversion

Interpeak diversion is a QC metric computed on the baseline CPM bigWigs that
identifies samples in which background signal is abnormally elevated relative to
domain signal. It is a relative, not absolute, measure — it captures samples
whose signal-to-background relationship is out of step with the rest of the
experiment group before any scaling has been applied. Because it operates on
baseline data, interpeak diversion reflects properties of the original sample
rather than artifacts of the normalization itself. All interpeak diversion
metrics are suppressed and set to NA for samples in genotype groups with fewer
than 3 replicates — see Replicate Count Requirements.

Three component flags are evaluated for each sample within its experiment group,
using robust z-scores (MAD-based):

1. *Interpeak high*: the sample's interpeak signal is elevated relative to the
   experiment-group median.
2. *Domain/interpeak low*: the sample's target region/interpeak ratio is reduced
   relative to the experiment-group median — enriched regions are weak compared
   to background.
3. *Anchor/domain high*: the sample's anchor/target region ratio is elevated
   relative to the experiment-group median — normalization anchors are
   dominating over the signal of interest.

CODA also records a related but separate baseline domain/interpeak outlier flag.
Unlike interpeak diversion, this flag is bidirectional: unusually low
domain/interpeak can indicate elevated background or weak domain enrichment,
whereas unusually high domain/interpeak can indicate a sample whose target signal
or background relationship is out of step with the rest of the cohort. This flag
does not contribute to the interpeak diversion signal count; it is treated as a
general baseline ratio outlier.

A sample is flagged for interpeak diversion (raw) when at least two of these
three conditions are satisfied simultaneously. Requiring convergence of multiple
independent indicators reduces false positives from noise in any single ratio.

However, some experiment contexts can produce spurious interpeak diversion flags
— specifically, experiments where background is genuinely very low (as is
typical for CUT&TAG), or where replicates are exceptionally tight (producing
inflated z-scores from small absolute differences). To address this, a
protection clause suppresses the interpeak diversion flag when the experiment-
group median interpeak signal is below 0.25 and the range of domain signal
across the group is below 0.30. The final `interpeak_diversion_flag` thus
reflects samples that show an abnormal background-to-signal pattern in a context
where that abnormality is meaningful.

The three output columns capture this logic explicitly:

- `interpeak_diversion_flag_raw`: mathematical detection without context
- `interpeak_diversion_protected`: flagged but likely benign given absolute
  signal context
- `interpeak_diversion_flag`: the final QC flag used for sample and experiment
  recommendations

Interpeak diversion complements the other QC metrics — it operates orthogonally
to dropout (which reflects anchor coverage failures) and genotype consistency
(which reflects within-genotype replicate agreement). A sample can pass all
other checks and still trigger interpeak diversion, and vice versa.

A note on peaky targets: for antibodies with sharp, punctate enrichment such as
transcription factors, the target region file will cover a smaller fraction of
the genome, and the domain/interpeak ratio may be inherently lower than for
broad chromatin marks. This does not affect the validity of the interpeak
diversion logic, which operates on relative differences between samples rather
than absolute signal levels. However, users should be aware that QC thresholds
developed primarily on broad-mark CUT&RUN data may benefit from adjustment when
applied to transcription factor experiments with very high signal-to-noise
ratios or very low genomic occupancy.

---

## Genotype Consistency

Within each experiment group, CODA evaluates whether replicates of the same
genotype are internally consistent with respect to their domain/interpeak and
anchor/interpeak ratios as measured in the baseline CPM bigWigs. This assessment
is deliberately performed on the pre-normalization data, so that genotype
inconsistency reflects sample-level differences already present before scaling.
The primary decision metric is the **max/min ratio spread**: a genotype group is
flagged for high spread when either the domain/interpeak or anchor/interpeak ratio
varies by at least the configured fold threshold across replicates (default
1.5-fold). CODA also reports the coefficient of variation of both ratios as a
complementary dispersion metric (`high_cv_flag`), but CV does **not** independently
determine the final `genotype_inconsistent_flag`.

As with interpeak diversion, a protection clause prevents false positives in
very low-background contexts. A high-spread genotype group is not called
inconsistent when median interpeak signal is very low and absolute domain signal
is tightly clustered, because the large ratio spread can then arise from small
denominator changes rather than meaningful replicate discordance.

If multiple genotype groups within an experiment show inconsistency, this
constitutes a more serious concern, because it suggests the problem is not
isolated to a single batch of samples.

All genotype consistency metrics are suppressed and set to NA for genotype
groups with fewer than 3 replicates — see Replicate Count Requirements.

---

## Class-Wide BASIN Assumption Violation

BASIN's normalization is predicated on the assumption that IgG nonspecific
enrichment at BARs is not systematically different between experimental
conditions. If a particular genotype or treatment class shows a coordinated
increase in interpeak signal alongside reduced domain/interpeak ratios — and
this pattern is seen in the majority of samples within that class — this may
indicate that the biological manipulation itself has altered the nonspecific
enrichment landscape. In such cases, the normalization may be working against
itself, scaling down samples with elevated background when the elevated
background is itself the signal.

CODA detects this by identifying genotype groups where the majority of samples
carry the interpeak diversion flag — itself computed on baseline CPM bigWigs —
the group's median baseline interpeak signal is above the experiment- group
median, and the group's baseline domain/interpeak ratio is below the experiment-
group median. When this pattern is detected in any genotype group, an
experiment-level flag is issued for manual inspection, because it implies that
one of BASIN's core assumptions may be violated for that class of samples.

---

## Post-BASIN Pipeline Integrity Check

HARMONY's canonical scale factor is `1 / median(raw coverage at retained BARs)`.
Because the canonical scaled bigWig uses the **same uncentered, unextended
coverage representation**, the retained-BAR median measured from that canonical
scaled output should be approximately 1.0.

For current HARMONY runs, CODA joins the canonical scaled BAR summary to the
exact retained coordinates recorded in `scales/retained_bars.tsv`. This avoids
any ambiguity from ties or rounding at the trim boundary. For legacy runs, CODA
reconstructs the set with the same global all-zero filter and
`round(N × drop_fraction)` rule used by HARMONY.

A large within-experiment CV of these retained-BAR medians (threshold
`POST_BASIN_ANCHOR_CV_WARN`, default 0.15) indicates a likely mechanical
problem such as the wrong BAR file, a scale factor that was not applied, or file
mismatch. It is not a biological quality metric.

The check intentionally uses only canonical scaled tracks under
`bigwigs_scaled/`. Centered or extended derivative tracks under
`bigwigs_downstream/` are never used; centering could move fixed read footprints
across BAR boundaries and would break the exact "median = 1 by construction"
identity even though the sample-level scalar itself remained valid.

Small departures from exactly 1 can arise from the independent
`multiBigwigSummary` measurement of the scaled files and chromosome filtering,
but CODA no longer attributes them to `bigWigAverageOverBed`; the current code
uses `multiBigwigSummary` for these regional summaries.

## Drop-Fraction Stability (advisory)

Choosing a drop fraction is a trade-off: trim too little and antibody-enriched
BARs inflate the background estimate; trim too much and you start removing
genuine background. CODA provides an advisory analysis to inform this choice,
built from the canonical depth-normalized CPM BAR summary (so it is independent
of the scaled output). CODA now applies the same global all-zero filter and
`round(N × drop_fraction)` trimming rule as HARMONY.

Two signals are reported per experiment group and drop fraction:

1. **Cross-replicate medBAR reproducibility.** For replicates of the same
   genotype, the depth-normalized background estimate (CPM medBAR) should agree. The drop that
   minimizes the within-genotype coefficient of variation of medBAR gives the
   most reproducible background. Computed when a genotype group has at least
   three replicates.
2. **Scale-factor elbow.** The per-sample relative change in the scale factor
   between consecutive drops. Where the curve flattens, additional trimming is
   removing genuine background rather than antibody-enriched anchors.

`secondary_outputs/drop_recommendation.tsv` reports, per experiment, the drop
that minimizes cross-replicate medBAR CV. This is **advisory only**: it
identifies the most reproducible background, but the authoritative choice should
be confirmed against a positive control (e.g. a known titration) showing that
genotype differences survive while replicate noise collapses. The supporting
tables are `drop_stability_by_experiment.tsv` and `drop_stability_by_genotype.tsv`.

---

## Sample-Level QC Recommendations

Each sample receives one of three classifications, based on a combination of
baseline and post-normalization metrics:

**ok.** No major issue detected under current QC rules.

**warning.** One or more of the following: the sample belongs to a genotype
group with fewer than 3 replicates (z-score metrics suppressed — see Replicate
Count Requirements); the sample dropped out at one or more BASIN thresholds
(i.e., no valid scale factor could be computed); baseline CPM signal behavior is
unusual relative to the experiment group (baseline weirdness); baseline 
domain/interpeak ratio is unusually high or low relative to the experiment 
group; anchor support is low in the baseline data; interpeak diversion is 
detected in the baseline data; the sample is an outlier in domain or interpeak 
signal among its experimentgroup in the scaled bigWigs at the primary drop 
threshold; or the retained anchor regions show unusual sex chromosome 
composition in the scaled data.

**suspect.** A sample is classified as suspect when two or more warning-level
conditions coincide — specifically, when dropout or baseline weirdness or
interpeak diversion is combined with low anchor support in the baseline data.
This combination suggests that the normalization is unlikely to have worked
correctly for this sample, because the BARs cannot provide a reliable reference
even before scaling is applied.

A plain-English note is generated for each sample explaining the specific basis
for its classification.

---

## Experiment-Level QC Recommendations

CODA issues one of four experiment-level recommendations describing the
appropriate next action. These recommendations are not grades assigned to the
experiment; they communicate what should happen next.

**Pass → proceed.** No major experiment-level QC trigger was detected.

**Pass with caveat → proceed while noting the limitation.** The experiment is
otherwise coherent, but one or more bounded limitations affect interpretation.
Examples include one or two warning-level samples or a genotype group with fewer
than three replicates, for which z-score-based metrics are suppressed. This
category does not indicate that normalization failed.

**Manual review → examine before proceeding.** Automated QC detected an issue
that requires scientific or technical judgment before the experiment is used.
Triggers include one or more suspect samples; three or more warning-level
samples, excluding warnings caused only by low replicate count; exactly one
genotype group with within-genotype inconsistency; a class-wide BASIN assumption
violation; two or more samples with interpeak diversion; post-BASIN anchor
nonconvergence; dropout requiring examination; or at least one-third of samples
flagged when those flags are confined to one genotype group. Manual review does
not itself mean that the experiment failed QC.

**Failed QC → do not proceed without remediation.** Strong evidence indicates
broad instability or a serious normalization problem. Triggers are unchanged:
dropout failed-QC thresholds reached; two or more genotype groups with
within-genotype inconsistency; domain/interpeak ratio spread exceeding 2.25
combined with genotype inconsistency; or at least one-third of samples flagged
across two or more genotype groups.

For each experiment, CODA records the specific decision trigger and a
plain-English note explaining the recommendation. It also reports
`n_flagged_genotype_groups` and `flagged_genotype_groups`, making the
distribution of sample-level flags explicit and auditable.


---

## Configurable Thresholds

All QC thresholds are exposed as environment variables and can be modified
before running CODA. The defaults reflect extensive empirical testing on
Drosophila CUT&RUN and CUT&TAG data and are recommended as starting points. 
The following table lists the key parameters, their defaults, and their roles.

**Baseline and anchor support**

| Parameter | Default | Role |
|---|---|---|
| `BASELINE_WEIRD_Z` | 2.5 | Ordinary within-experiment z-score threshold for anomalous baseline signal |
| `LOW_ANCHOR_FRAC_MIN` | 0.15 | Min fraction of BARs with nonzero signal |
| `LOW_ANCHOR_COUNT_MIN` | 50 | Min number of nonzero BAR regions |
| `BASELINE_RATIO_OUTLIER_Z` | 2.0 | Robust z-score for bidirectional baseline domain/interpeak outlier |

**Interpeak diversion**

| Parameter | Default | Role |
|---|---|---|
| `INTERPEAK_DIVERSION_Z` | 2.0 | Robust z-score for diversion components |
| `INTERPEAK_DIVERSION_SIGNALS_MIN` | 2 | Signals required to trigger flag |
| `INTERPEAK_DIVERSION_LOW_INTERPEAK_PROTECT_MAX` | 0.25 | Interpeak ceiling for protection |
| `INTERPEAK_DIVERSION_DOMAIN_RANGE_PROTECT_MAX` | 0.30 | Domain range ceiling for protection |

**Genotype consistency and post-normalization outliers**

| Parameter | Default | Role |
|---|---|---|
| `GENOTYPE_SPREAD_RATIO_WARN` | 1.5 | Max/min ratio spread threshold |
| `GENOTYPE_SPREAD_CV_WARN` | 0.25 | Advisory CV threshold reported as `high_cv_flag`; does not independently determine `genotype_inconsistent_flag` |
| `PRIMARY_OUTLIER_Z` | 2.0 | Ordinary within-experiment post-normalization outlier z-score |
| `SEX_CHROM_SKEW_Z` | 2.0 | Sex chromosome anchor skew z-score |
| `POST_BASIN_ANCHOR_CV_WARN` | 0.15 | Pipeline integrity: within-group CV of scaled trimmed anchor medians above which a group is flagged (indicates a mechanical issue, not biology) |

**Experiment-level Failed QC rules**

| Parameter | Default | Role |
|---|---|---|
| `DROP10_FAIL_FRAC` | 0.25 | Dropout fraction at drop0.10 for Failed QC |
| `DROP50_FAIL_FRAC` | 0.50 | Dropout fraction at drop0.50 for Failed QC |
| `FAIL_INCONSISTENT_GENOTYPE_GROUPS` | 2 | Inconsistent groups to trigger Failed QC |
| `FAIL_MAX_DOMAIN_RATIO_SPREAD` | 2.25 | Domain/interpeak spread for Failed QC |
| `FAIL_FRACTION_FLAGGED` | 0.333 | Flagged fraction evaluated together with flag distribution: Failed QC only when flags span at least 2 genotype groups; Manual review when confined to 1 genotype group |
| `WARNING_SAMPLE_MAX` | 2 | Max warning samples before escalating |

---

## Outputs

CODA produces the following files in a `qc_[label]/` directory within the run
folder. The two **primary** outputs sit at the top of `qc_[label]/`; all
**secondary** tables and figures are written into a `qc_[label]/secondary_outputs/`
subdirectory so the directory root stays uncluttered.

### Primary Outputs

These two files should be examined for every run. They consolidate all QC
signals into actionable per-sample and per-experiment summaries and are designed
to be the first and often only stop for routine quality assessment.

- `experiment_qc_summary.tsv`: One row per experiment group. Contains the
  experiment-level recommendation (Pass / Pass with caveat / Manual review /
  Failed QC), the specific decision trigger, and a plain-English note explaining
  the basis for the recommendation. CODA additionally reports
  `n_flagged_genotype_groups` and `flagged_genotype_groups`, which distinguish
  broadly distributed sample flags from a coherent pattern confined to one
  genotype group. The column currently named `recommended_drop` records the
  configured **primary drop** used for post-BASIN QC when that threshold is
  present; it is not the advisory drop recommendation from the drop-stability
  analysis. The advisory recommendation is reported separately in
  `secondary_outputs/drop_recommendation.tsv` as
  `suggested_drop_by_medbar_cv`. This is the top-level summary of whether the
  normalization worked.

- `sample_qc_summary.tsv`: One row per sample. Contains the sample-level
  classification (ok / warning / suspect), a pipe-delimited list of all flag
  components that contributed to the classification, and a plain-English note.
  This file consolidates the key per-sample signals from all secondary tables
  into a single place, so that most diagnostic questions can be answered without
  leaving this file.

### Secondary Outputs

These files provide the underlying metrics that feed into the primary summaries.
They are most useful when a sample or experiment has been flagged and you need
to understand why. All paths below are relative to
`qc_[label]/secondary_outputs/`.

**Tables:**
- `post_basin_anchor_check.tsv`: Pipeline integrity check on canonical scaled
  tracks. One row per sample at the primary drop: the retained scaled anchor
  median (`primary_anchor_trimmed`, expected ~1.0), experiment-group median/CV,
  and nonconvergence flag. Current runs use exact retained BAR identities from
  HARMONY; legacy runs use corrected reconstruction.
- `drop_stability_by_experiment.tsv`: Per experiment group and drop fraction, the
  mean within-genotype CV of depth-normalized (CPM) medBAR (cross-replicate reproducibility of the
  background estimate) and the median per-sample scale-factor step between
  consecutive drops (the elbow signal). Lower CV = more reproducible background
- `drop_stability_by_genotype.tsv`: The granular table behind the above — raw
  medBAR mean, sd, and CV for each genotype group at each drop (CV reported when
  n ≥ 3)
- `drop_recommendation.tsv`: Advisory per-experiment suggested drop (the one
  minimizing cross-replicate medBAR CV), alongside the run's primary drop. Treat
  as advisory; confirm against a positive control before locking quantitative
  cross-sample claims
- `sample_qc_by_threshold.tsv`: Per-sample signal metrics (domain, anchor,
  interpeak medians and ratios) at each drop threshold, from the scaled bigWigs
- `baseline_ratio_table.tsv`: Pre-normalization regional signal ratios from the
  baseline CPM bigWigs
- `interpeak_diversion_metrics.tsv`: Per-sample interpeak diversion component
  flags and robust z-scores, computed on baseline data
- `genotype_spread_qc.tsv`: Within-genotype consistency metrics, computed on
  baseline data
- `retained_anchor_chromosome_fraction_by_threshold.tsv`: Chromosome composition
  of retained BARs across drop thresholds

**Figures:**
- `threshold_domain_interpeak.png`: Domain/interpeak ratio across BASIN drop
  thresholds, by sample and genotype — useful for visualizing how signal
  relationships change as more anchors are trimmed
- `threshold_anchor_interpeak.png`: Anchor/interpeak ratio across BASIN
  thresholds
- `baseline_median_signal_grouped_barplot.png`: Pre-normalization median signal
  by region class and sample
- `baseline_ratio_grouped_barplot.png`: Pre-normalization signal ratios by
  sample
- `primary_drop_median_signal_grouped_barplot.png`: Post-normalization median
  signal by region class and sample at the primary drop threshold
- `primary_drop_ratio_grouped_barplot.png`: Post-normalization signal ratios at
  the primary drop threshold
- `retained_anchor_chromosome_fraction_by_threshold.png`: Chromosome composition
  of retained BARs across thresholds — useful for diagnosing sex chromosome
  anchor skew
