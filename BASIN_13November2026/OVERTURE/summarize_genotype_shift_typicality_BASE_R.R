#!/usr/bin/env Rscript

# ==============================================================================
# OVERTURE secondary QC: medBAR genotype-shift typicality among passing lt50 panels
#
# DEPENDENCIES
# ------------
# Base R only. No tidyverse / ggplot2 installation is required.
#
# PURPOSE
# -------
# This is a POST-EVALUATION diagnostic. It does NOT change OVERTURE pass/fail
# criteria and it does NOT select a panel automatically.
#
# It asks:
#   Among panels that already passed the primary OVERTURE criteria, how closely
#   does each panel follow the central medBAR genotype-shift pattern across the
#   validation marks and genotypes?
#
# medBAR GENOTYPE SHIFT
# ---------------------
# OVERTURE evaluates candidate BAR panels by measuring the trimmed median signal
# across BARs (medBAR) for each sample. medBAR is intended to estimate the
# nonspecific/background component of the assay rather than enrichment at true
# target sites.
#
# For each mark, OVERTURE calculates a medBAR genotype shift:
#
#     medBAR genotype shift =
#       median medBAR(test genotype) - median medBAR(reference genotype)
#
# Thus:
#   negative shift = lower BAR/background signal in the test genotype
#   positive shift = higher BAR/background signal in the test genotype
#
# This quantity should NOT be interpreted simply as a change in biological
# target signal. Instead, it measures whether the background distribution
# sampled by a candidate BAR panel has moved between genotypes.
#
# WHY medBAR CAN SHIFT WHEN BIOLOGY SHIFTS
# ----------------------------------------
# If there is no true global difference in specific signal between genotypes,
# and candidate BARs are sampling nonspecific background appropriately, the
# medBAR distributions should remain broadly comparable between genotypes.
# Individual BARs and individual panels will vary, but the ensemble of passing
# panels should not show a strong coordinated displacement.
#
# If a mark undergoes a genuine global change, however, the background can move
# even when the BARs themselves are not true target sites. Sequencing data are
# compositional: when a larger fraction of reads is devoted to true specific
# signal, the relative representation of nonspecific background can decrease;
# when global specific signal decreases, background can become relatively more
# represented.
#
# Consequently, a genuine global biological change may produce a coordinated
# displacement of medBAR across many independent candidate panels. This is the
# background-titration behavior that BASIN/HARMONY is intended to detect and
# accommodate.
#
# INTERPRETING medBAR GENOTYPE SHIFT
# ----------------------------------
# A medBAR genotype shift can therefore arise from at least two conceptually
# different situations:
#
#   1. PANEL-SPECIFIC BEHAVIOR
#      One candidate panel shifts differently from the other passing panels.
#      This may indicate that the panel contains regions affected by real
#      biological signal or otherwise behaves atypically.
#
#   2. DATASET-WIDE / BIOLOGICAL BEHAVIOR
#      Many passing panels shift together in the same direction for a particular
#      mark/genotype comparison. This can be consistent with a genuine global
#      change in that mark because the whole background distribution is
#      displaced by the changing balance between specific signal and
#      nonspecific background.
#
# For this reason, medBAR genotype shift is a DIAGNOSTIC rather than something
# that should simply be minimized toward zero.
#
# GENOTYPE-SHIFT TYPICALITY
# -------------------------
# Typicality asks whether a candidate panel follows the medBAR genotype-shift
# pattern seen across the ensemble of passing panels.
#
# It does NOT ask whether the shift is close to zero.
#
# A panel can therefore have:
#
#   large medBAR shift + low typicality distance
#
# if many other passing panels show the same shift. Such a pattern is consistent
# with a dataset-wide effect and may reflect real global biology.
#
# Conversely:
#
#   small or moderate medBAR shift + high typicality distance
#
# can occur if one panel behaves differently from the rest of the ensemble,
# which is potentially more concerning for panel robustness.
#
# TYPICALITY CALCULATION
# ----------------------
# Within each mark x genotype comparison:
#   1. Calculate the median signed medBAR genotype shift across passing panels.
#   2. Calculate the MAD (median absolute deviation; constant = 1).
#   3. For each panel calculate:
#
#        | panel_shift - facet_median | / facet_MAD
#
# For each panel, these standardized distances are summarized across all
# mark x genotype comparisons:
#
#   mean_typicality_distance
#     Mean standardized distance from the ensemble center.
#     Lower = more representative overall.
#     Sensitive to occasional large deviations.
#
#   median_typicality_distance
#     Median standardized distance from the ensemble center.
#     Lower = more representative in a typical comparison.
#     Less sensitive to isolated outlying behavior.
#
# IMPORTANT:
# A panel can show a LARGE real medBAR genotype shift and still have a LOW
# typicality distance if the other passing panels show the same biological shift.
#
# Because typicality is defined relative to the evaluation dataset, its value
# depends on which marks, genotypes, and biological perturbations were included.
# It is therefore a secondary robustness diagnostic, not a universal property of
# a BAR panel and not an independent optimization target.
#
# FLEXIBILITY ACROSS EVALUATION DATASETS
# --------------------------------------
# Marks and genotypes are discovered directly from genotype_shift_trimmed.tsv;
# they are not required to match the dm6 example dataset. The script also
# detects the delta_vs_<reference> column automatically, so the configured
# reference genotype may differ between OVERTURE runs.
#
# Plot dimensions scale with the number of marks and genotypes. A few familiar
# dm6 names are given a preferred display order when present; any additional
# marks or genotypes are retained automatically.
#
# INPUT DISCOVERY
# ---------------
# No argument is required when running from the OVERTURE/ directory,
# or when setting OVERTURE_ROOT to that directory.
#
# It automatically searches the OVERTURE tree for evaluation-output directories
# containing BOTH:
#   genotype_shift_trimmed.tsv
#   panel_selection_summary.tsv
#
# It chooses the most recently modified matching directory and PRINTS that path
# prominently before analysis.
#
# Optional explicit override:
#   Rscript summarize_genotype_shift_typicality_BASE_R.R /path/to/panel_eval_out
#
# or:
#   PANEL_EVAL_DIR=/path/to/panel_eval_out \
#     Rscript summarize_genotype_shift_typicality_BASE_R.R
#
# OUTPUTS
# -------
# Written into:
#   <panel_eval_dir>/genotype_shift_typicality_qc/
#
#   panel_genotype_shift_typicality.tsv
#   genotype_shift_typicality_detail.tsv
#   genotype_shift_passed_panels_signed.png
#
# The plot uses ONE COMMON Y AXIS across all mark x genotype panels.
# ==============================================================================

CHOSEN_DROP <- Sys.getenv("EVAL_DROP_FRAC", unset = "median_drop20")
LT50_PATTERN <- "filter_away_lt50bp"

panel_core <- function(x) {
  x <- basename(x)
  x <- sub("\\.chr\\.bed$", "", x)
  x <- sub("\\.nochr\\.bed$", "", x)
  x <- sub("\\.bed$", "", x)
  x
}

valid_eval_dir <- function(d) {
  dir.exists(d) &&
    file.exists(file.path(d, "genotype_shift_trimmed.tsv")) &&
    file.exists(file.path(d, "panel_selection_summary.tsv"))
}

find_overture_root <- function() {
  env_root <- Sys.getenv("OVERTURE_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root))
  }

  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)

  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_path <- tryCatch(normalizePath(script_path), error = function(e) NA_character_)
    if (!is.na(script_path)) {
      maybe_root <- dirname(script_path)
      if (file.exists(file.path(maybe_root, "run_overture.sh"))) {
        return(maybe_root)
      }
      maybe_parent <- dirname(maybe_root)
      if (file.exists(file.path(maybe_parent, "run_overture.sh"))) {
        return(maybe_parent)
      }
    }
  }

  p <- normalizePath(getwd())
  repeat {
    if (file.exists(file.path(p, "run_overture.sh")) &&
        dir.exists(file.path(p, "step4_panel_evaluation"))) {
      return(p)
    }
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }

  stop(
    "Could not identify the OVERTURE root. Place this script in ",
    "OVERTURE/, run from inside OVERTURE, ",
    "or set OVERTURE_ROOT."
  )
}

discover_eval_dir <- function(root) {
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  candidates <- dirs[vapply(dirs, valid_eval_dir, logical(1))]

  if (length(candidates) == 0) {
    stop(
      "No evaluation directory containing both required TSV files was found under:\n",
      root
    )
  }

  scores <- vapply(candidates, function(d) {
    a <- file.info(file.path(d, "genotype_shift_trimmed.tsv"))$mtime
    b <- file.info(file.path(d, "panel_selection_summary.tsv"))$mtime
    as.numeric(max(a, b, na.rm = TRUE))
  }, numeric(1))

  candidates[which.max(scores)]
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  stop("Usage: Rscript summarize_genotype_shift_typicality_BASE_R.R [panel_eval_dir]")
}

env_override <- Sys.getenv("PANEL_EVAL_DIR", unset = "")

if (length(args) == 1) {
  panel_dir <- normalizePath(args[1], mustWork = TRUE)
} else if (nzchar(env_override)) {
  panel_dir <- normalizePath(env_override, mustWork = TRUE)
} else {
  root <- find_overture_root()
  panel_dir <- discover_eval_dir(root)
}

if (!valid_eval_dir(panel_dir)) {
  stop(
    "Selected directory is missing genotype_shift_trimmed.tsv and/or ",
    "panel_selection_summary.tsv:\n", panel_dir
  )
}

cat("\n============================================================\n")
cat("OVERTURE genotype-shift typicality QC\n")
cat("============================================================\n")
cat("Using evaluation directory:\n  ", panel_dir, "\n", sep = "")
cat("Drop fraction:\n  ", CHOSEN_DROP, "\n", sep = "")
cat("Panel subset:\n  ", LT50_PATTERN, "\n", sep = "")
cat("============================================================\n\n")

shift <- read.delim(
  file.path(panel_dir, "genotype_shift_trimmed.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)

panel_summary <- read.delim(
  file.path(panel_dir, "panel_selection_summary.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)

required_shift <- c("drop", "panel", "mark", "genotype")
missing_shift <- setdiff(required_shift, names(shift))
if (length(missing_shift) > 0) {
  stop(
    "genotype_shift_trimmed.tsv is missing required columns: ",
    paste(missing_shift, collapse = ", ")
  )
}

# Detect the genotype-shift column automatically so this works with any
# configured reference genotype (for example delta_vs_yw, delta_vs_WT, etc.).
delta_cols <- grep("^delta_vs_", names(shift), value = TRUE)

if (length(delta_cols) != 1) {
  stop(
    "Expected exactly one genotype-shift column named delta_vs_<reference>, but found: ",
    paste(delta_cols, collapse = ", ")
  )
}

delta_col <- delta_cols[1]
reference_genotype <- sub("^delta_vs_", "", delta_col)

cat("Detected reference genotype: ", reference_genotype, "\n", sep = "")
cat("Detected genotype-shift column: ", delta_col, "\n\n", sep = "")

if (!"panel" %in% names(panel_summary)) {
  stop("panel_selection_summary.tsv is missing required column: panel")
}

shift$panel_key <- panel_core(shift$panel)
panel_summary$panel_key <- panel_core(panel_summary$panel)

passing_keys <- unique(
  panel_summary$panel_key[grepl(LT50_PATTERN, panel_summary$panel_key, fixed = TRUE)]
)

if (length(passing_keys) < 2) {
  stop("Fewer than two passing lt50 panels were found; typicality needs an ensemble.")
}

keep <- shift$drop == CHOSEN_DROP &
        shift$panel_key %in% passing_keys &
        grepl(LT50_PATTERN, shift$panel_key, fixed = TRUE)

shift_pass <- shift[keep, , drop = FALSE]

if (nrow(shift_pass) == 0) {
  stop("No genotype-shift rows matched passing lt50 panels at drop = ", CHOSEN_DROP)
}

cat("Passing lt50 panels in panel_selection_summary.tsv: ",
    length(passing_keys), "\n", sep = "")
cat("Passing lt50 panels represented in genotype-shift data: ",
    length(unique(shift_pass$panel_key)), "\n\n", sep = "")

missing_from_shift <- setdiff(passing_keys, unique(shift_pass$panel_key))
if (length(missing_from_shift) > 0) {
  warning(
    "Some passing lt50 panels have no genotype-shift rows at ", CHOSEN_DROP,
    ":\n  ", paste(missing_from_shift, collapse = "\n  ")
  )
}

# ------------------------------------------------------------------------------
# Calculate facet medians, MADs, and standardized typicality distance
# ------------------------------------------------------------------------------

facet_id <- interaction(
  shift_pass$mark,
  shift_pass$genotype,
  drop = TRUE,
  lex.order = TRUE
)

facet_median <- ave(
  shift_pass[[delta_col]],
  facet_id,
  FUN = function(x) median(x, na.rm = TRUE)
)

facet_mad <- ave(
  shift_pass[[delta_col]],
  facet_id,
  FUN = function(x) {
    m <- mad(x, center = median(x, na.rm = TRUE), constant = 1, na.rm = TRUE)
    if (is.na(m) || m == 0) 1 else m
  }
)

shift_pass$facet_median_shift <- facet_median
shift_pass$facet_mad_shift <- facet_mad
shift_pass$typicality_distance <- abs(
  shift_pass[[delta_col]] - shift_pass$facet_median_shift
) / shift_pass$facet_mad_shift

panels <- unique(shift_pass$panel_key)

summary_rows <- lapply(panels, function(pk) {
  z <- shift_pass[shift_pass$panel_key == pk, , drop = FALSE]
  data.frame(
    panel = z$panel[1],
    panel_key = pk,
    mean_typicality_distance = mean(z$typicality_distance, na.rm = TRUE),
    median_typicality_distance = median(z$typicality_distance, na.rm = TRUE),
    n_mark_genotype_comparisons = sum(!is.na(z$typicality_distance)),
    stringsAsFactors = FALSE
  )
})

typicality_summary <- do.call(rbind, summary_rows)

typicality_summary$rank_mean <- rank(
  typicality_summary$mean_typicality_distance,
  ties.method = "min",
  na.last = "keep"
)

typicality_summary$rank_median <- rank(
  typicality_summary$median_typicality_distance,
  ties.method = "min",
  na.last = "keep"
)

typicality_summary <- typicality_summary[
  order(
    typicality_summary$rank_mean,
    typicality_summary$rank_median
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------------------------
# Write tables
# ------------------------------------------------------------------------------

out_dir <- file.path(panel_dir, "genotype_shift_typicality_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_out <- file.path(
  out_dir,
  "panel_genotype_shift_typicality.tsv"
)

detail_out <- file.path(
  out_dir,
  "genotype_shift_typicality_detail.tsv"
)

write.table(
  typicality_summary,
  summary_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  shift_pass,
  detail_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Signed passing-panel violin plot, BASE R
# ------------------------------------------------------------------------------

marks <- unique(shift_pass$mark)
genotypes <- unique(shift_pass$genotype)

# Keep stable, readable ordering when these expected names are present.
preferred_marks <- c(
  "H3K27ac", "H3K36me2", "H3K36me3",
  "panH4", "RNAPIIS2P", "RNAPIIS5P"
)
marks <- c(
  preferred_marks[preferred_marks %in% marks],
  setdiff(sort(marks), preferred_marks)
)

preferred_genotypes <- c("H33AnullHWT", "Set2null", "TripleK36R")
genotypes <- c(
  preferred_genotypes[preferred_genotypes %in% genotypes],
  setdiff(sort(genotypes), preferred_genotypes)
)

yrange <- range(shift_pass[[delta_col]], na.rm = TRUE)
pad <- diff(yrange) * 0.06
if (!is.finite(pad) || pad == 0) pad <- 0.1
ylim <- c(yrange[1] - pad, yrange[2] + pad)

plot_out <- file.path(
  out_dir,
  "genotype_shift_passed_panels_signed.png"
)

# Scale the raster dimensions with the number of facets so runs containing
# different numbers of marks or genotypes remain legible.
plot_width_in <- max(8, 3.2 * length(genotypes))
plot_height_in <- max(6, 1.65 * length(marks))

png(
  filename = plot_out,
  width = plot_width_in,
  height = plot_height_in,
  units = "in",
  res = 300
)

old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

par(
  mfrow = c(length(marks), length(genotypes)),
  mar = c(1.7, 3.2, 2.0, 0.8),
  oma = c(3.0, 4.0, 4.0, 1.0)
)

draw_violin <- function(vals, ylim, show_y = TRUE) {
  vals <- vals[is.finite(vals)]

  plot(
    NA,
    xlim = c(0.55, 1.45),
    ylim = ylim,
    xaxt = "n",
    yaxt = if (show_y) "s" else "n",
    xlab = "",
    ylab = "",
    bty = "l"
  )

  abline(h = 0, lty = 2, col = "grey55")

  if (length(vals) >= 2 && length(unique(vals)) > 1) {
    d <- density(vals, na.rm = TRUE)
    w <- d$y / max(d$y) * 0.32
    polygon(
      c(1 - w, rev(1 + w)),
      c(d$x, rev(d$x)),
      col = "grey90",
      border = "grey55"
    )
  }

  if (length(vals) > 0) {
    q <- quantile(vals, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    segments(0.88, q[1], 1.12, q[1], lwd = 1)
    segments(0.88, q[3], 1.12, q[3], lwd = 1)
    rect(0.90, q[1], 1.10, q[3], col = "white", border = "black")
    segments(0.90, q[2], 1.10, q[2], lwd = 2)

    # Deterministic horizontal jitter
    offsets <- seq(-0.10, 0.10, length.out = length(vals))
    points(1 + offsets, vals, pch = 16, cex = 0.75, col = grDevices::adjustcolor("black", alpha.f = 0.55))
  }
}

for (i in seq_along(marks)) {
  for (j in seq_along(genotypes)) {
    vals <- shift_pass[[delta_col]][
      shift_pass$mark == marks[i] &
      shift_pass$genotype == genotypes[j]
    ]

    draw_violin(vals, ylim, show_y = (j == 1))

    if (i == 1) {
      title(main = genotypes[j], cex.main = 0.85)
    }
    if (j == length(genotypes)) {
      mtext(marks[i], side = 4, line = 0.1, cex = 0.75)
    }
  }
}

mtext(
  "medBAR genotype-shift distributions across passing OVERTURE panels",
  side = 3, outer = TRUE, line = 2.2, cex = 1.15, font = 2
)

mtext(
  paste0(
    "Signed medBAR shift relative to ",
    reference_genotype,
    "; ",
    CHOSEN_DROP,
    "; lt50-filtered panels"
  ),
  side = 3, outer = TRUE, line = 0.8, cex = 0.85
)

mtext(
  "medBAR genotype shift",
  side = 2, outer = TRUE, line = 2.2, cex = 0.95
)

dev.off()

# ------------------------------------------------------------------------------
# Report
# ------------------------------------------------------------------------------

best_mean <- typicality_summary[
  which.min(typicality_summary$mean_typicality_distance),
  ,
  drop = FALSE
]

best_median <- typicality_summary[
  which.min(typicality_summary$median_typicality_distance),
  ,
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("Typicality summary\n")
cat("------------------------------------------------------------\n")

cat("\nMost representative by MEAN typicality distance:\n")
cat("  (lowest mean standardized deviation from the passing-panel ensemble)\n")
print(best_mean, row.names = FALSE)

cat("\nMost representative by MEDIAN typicality distance:\n")
cat("  (lowest median standardized deviation in a typical mark x genotype comparison)\n")
print(best_median, row.names = FALSE)

cat("\nInterpretation reminder:\n")
cat("  Typicality is relative to the marks/genotypes used for evaluation.\n")
cat("  A real global biological shift can be highly typical if it is reproduced\n")
cat("  across the passing-panel ensemble. This QC does not alter panel pass/fail\n")
cat("  status and should not be treated as an independent optimization target.\n")

cat("\nOutputs written to:\n  ", out_dir, "\n", sep = "")
cat("  ", basename(summary_out), "\n", sep = "")
cat("  ", basename(detail_out), "\n", sep = "")
cat("  ", basename(plot_out), "\n", sep = "")
cat("\nDone.\n")
