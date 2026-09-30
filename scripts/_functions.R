# Shared loaders, tests and figure style for the analysis scripts in scripts/
# (01_four_groups.Rmd to 05_sensitivity.Rmd). The loaders read the harmonised
# inputs in ../data/, which 00_harmonise_cohorts.Rmd writes from ../raw/.
# Everything here is written plainly, one step per line, so it can be read
# top to bottom. Each script sources this file and then does its own work.

# ------------------------------------------------------------- constants ----

CELL_LEVELS = c("EAS Never", "EAS Ever", "Western Never", "Western Ever")

CONTRAST_PAIRS = list(
  c("EAS Never", "Western Never"),
  c("EAS Ever", "Western Ever"),
  c("EAS Never", "EAS Ever"),
  c("Western Never", "Western Ever"),
  c("EAS Ever", "Western Never"))

CONTRAST_LABELS = c(
  "ancestry within never-smokers",
  "ancestry within ever-smokers",
  "smoking within EAS",
  "smoking within Western",
  "EAS ever vs Western never")

# Labelling sets only. Every analysis runs genome-wide or over all signatures;
# these are the features that get a text label on a figure.
LABEL_GENES = c("EGFR", "KRAS", "KEAP1", "TP53", "STK11", "RBM10", "BRAF",
                "PIK3CA", "SMARCA4", "NF1", "PTEN", "RB1", "ERBB2", "MET", "ALK")
LABEL_SIGNATURES = c("SBS1", "SBS2", "SBS4", "SBS5", "SBS13", "SBS18", "SBS40a")

SIGNATURE_COHORTS = c("Chen2020", "Xu2020", "Guo2026", "TRACERx2023", "TCGA2014")
ALL_COHORTS = c("Chen2019", "Chen2020", "Xu2020", "Guo2026", "TCGA2014", "TRACERx2023")

# Cohorts with enough tumours in both smoking strata for within-cohort work.
PER_COHORT_SET = c("Chen2020", "Xu2020", "TCGA2014", "TRACERx2023")

MIN_GENE_MUTATED = 10
MIN_SIGNATURE_ACTIVE = 20

# The 39 autosomal arms, in genome order. Acrocentric p arms are excluded.
ARM_ORDER = c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q",
              "6p", "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q",
              "11p", "11q", "12p", "12q", "13q", "14q", "15q", "16p", "16q",
              "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21q", "22q")

# ----------------------------------------------------------- figure style ----

CELL_COLS = c("EAS Never" = "#3B6FB6", "EAS Ever" = "#9FC0E4",
              "Western Never" = "#C0392B", "Western Ever" = "#EFA9A2")
GROUP_COLS = c("EAS" = "#3B6FB6", "Western" = "#C0392B")
SMOKE_COLS = c("Never" = "#7FB3D5", "Ever" = "#B03A2E")
CLASS_COLS = c("cross-ancestry" = "#C0392B", "within-EAS" = "#3B6FB6",
               "within-Western" = "#7F8C8D")
DIR_COLS = c("EAS-enriched" = "#3B6FB6", "Western-enriched" = "#C0392B",
             "ns" = "grey78")
CNV_COLS = c("Gain" = "#C0392B", "Loss" = "#3B6FB6")

theme_pub = function(base = 9) {
  out = ggplot2::theme_bw(base_size = base)
  out = out + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  out = out + ggplot2::theme(panel.grid.major =
    ggplot2::element_line(linewidth = 0.25, colour = "grey90"))
  out = out + ggplot2::theme(panel.border =
    ggplot2::element_rect(linewidth = 0.4, colour = "grey30"))
  out = out + ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"))
  out = out + ggplot2::theme(strip.background =
    ggplot2::element_rect(fill = "grey95", linewidth = 0.3))
  out = out + ggplot2::theme(strip.text =
    ggplot2::element_text(face = "bold", size = base - 0.5))
  out = out + ggplot2::theme(legend.key.size = grid::unit(0.35, "cm"))
  out = out + ggplot2::theme(plot.title =
    ggplot2::element_text(face = "bold", size = base + 1))
  out
}

# ---------------------------------------------------------------- loading ----

load_metadata = function(drop_ancestry = TRUE, drop_hyper = FALSE) {
  meta = read.csv("../data/analysis_metadata.csv", stringsAsFactors = FALSE)
  # The in-house Guo2026 variant calls are restricted and ship separately.
  # Without them Guo2026 tumours cannot enter the mutation-frequency analyses.
  if (!file.exists("../data/Guo2026_snv.csv")) {
    warning("../data/Guo2026_snv.csv not found. The in-house Guo2026 variant ",
            "calls are restricted and not distributed with this release, so ",
            "Guo2026 tumours are left out of the mutation-frequency analyses. ",
            "Request the file from the corresponding author and place it in ",
            "../data/ to reproduce the published numbers exactly.",
            call. = FALSE)
    meta$HasFullMAF[meta$Cohort == "Guo2026"] = FALSE
  }
  if (drop_ancestry) {
    meta = meta[!meta$AncestryExcluded, ]
  }
  if (drop_hyper) {
    meta = meta[!(meta$HypermutatedMechanism %in% TRUE), ]
  }
  meta = meta[!is.na(meta$Smoking), ]
  meta$Cell = factor(paste(meta$CohortGroup, meta$Smoking), levels = CELL_LEVELS)
  meta
}

# Unique sample-by-gene table of protein-altering variants only (R2-4).
load_protein_altering = function(cohorts) {
  out = NULL
  for (i in seq_along(cohorts)) {
    d = data.table::fread(paste0("../data/", cohorts[i], "_snv.csv"),
                          data.table = FALSE, showProgress = FALSE)
    d = d[d$VariantClass == "protein_altering", ]
    d = d[, c("SampleID", "Gene")]
    out = rbind(out, d)
  }
  out = out[!duplicated(out), ]
  out
}

# Per-sample relative signature contributions, one row per sample.
load_signature_relative = function(cohorts = SIGNATURE_COHORTS) {
  out = NULL
  common = NULL
  parts = list()
  for (i in seq_along(cohorts)) {
    w = read.csv(paste0("../data/signature/", cohorts[i],
                        "_signature_weights.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
    rownames(w) = w$Signature
    m = as.matrix(w[, setdiff(names(w), "Signature"), drop = FALSE])
    total = colSums(m)
    for (j in seq_len(ncol(m))) {
      if (total[j] > 0) {
        m[, j] = m[, j] / total[j]
      } else {
        m[, j] = NA
      }
    }
    d = as.data.frame(t(m))
    d$SampleID = rownames(d)
    parts[[i]] = d
    sigs = setdiff(names(d), "SampleID")
    if (is.null(common)) {
      common = sigs
    } else {
      common = intersect(common, sigs)
    }
  }
  for (i in seq_along(parts)) {
    out = rbind(out, parts[[i]][, c("SampleID", common)])
  }
  out
}

# Signatures active above 1 per cent in at least min_active tumours.
active_signatures = function(sig_wide, min_active = MIN_SIGNATURE_ACTIVE) {
  sigs = setdiff(names(sig_wide), "SampleID")
  keep = c()
  for (i in seq_along(sigs)) {
    v = sig_wide[[sigs[i]]]
    if (sum(v > 0.01, na.rm = TRUE) >= min_active) {
      keep = c(keep, sigs[i])
    }
  }
  keep
}

# Per-sample counts of the six substitution classes, plus their fractions.
load_point_mutation_classes = function(cohorts = SIGNATURE_COHORTS) {
  classes = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  out = NULL
  for (i in seq_along(cohorts)) {
    sp = read.csv(paste0("../data/signature/", cohorts[i],
                         "_spectrum96.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
    samples = setdiff(names(sp), "Type")
    m = as.matrix(sp[, samples, drop = FALSE])
    d = data.frame(SampleID = samples, stringsAsFactors = FALSE)
    d$TotalSBS = colSums(m)
    for (k in seq_along(classes)) {
      rows = grepl(paste0("\\[", classes[k], "\\]"), sp$Type)
      d[[classes[k]]] = colSums(m[rows, , drop = FALSE]) / d$TotalSBS
    }
    out = rbind(out, d)
  }
  out
}

# Aggregated 96-context spectrum, one row per context per cohort.
load_sbs96_long = function(cohorts = SIGNATURE_COHORTS) {
  out = NULL
  for (i in seq_along(cohorts)) {
    sp = read.csv(paste0("../data/signature/", cohorts[i],
                         "_spectrum96.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
    samples = setdiff(names(sp), "Type")
    m = as.matrix(sp[, samples, drop = FALSE])
    d = data.frame(Cohort = cohorts[i], Type = sp$Type,
                   Count = rowSums(m), stringsAsFactors = FALSE)
    out = rbind(out, d)
  }
  out$MutationClass = substr(out$Type, 3, 5)
  out
}

# ------------------------------------------------------------------ tests ----

# Two-group Wilcoxon, returned as one row.
test_two_group = function(value, group, level_a, level_b, label) {
  a = value[group == level_a]
  b = value[group == level_b]
  a = a[!is.na(a)]
  b = b[!is.na(b)]
  if (length(a) < 3 || length(b) < 3) {
    return(NULL)
  }
  tst = wilcox.test(a, b)
  out = data.frame(feature = label, group_a = level_a, group_b = level_b,
                   n_a = length(a), n_b = length(b),
                   median_a = median(a), median_b = median(b),
                   mean_a = mean(a), mean_b = mean(b),
                   stringsAsFactors = FALSE)
  out$log2_ratio = log2((out$median_a + 1e-3) / (out$median_b + 1e-3))
  out$p_value = tst$p.value
  out
}

# Two-group Fisher on a binary feature, returned as one row.
test_two_group_binary = function(hit, group, level_a, level_b, label) {
  a = hit[group == level_a]
  b = hit[group == level_b]
  a = a[!is.na(a)]
  b = b[!is.na(b)]
  if (length(a) < 3 || length(b) < 3) {
    return(NULL)
  }
  tab = matrix(c(sum(a), length(a) - sum(a), sum(b), length(b) - sum(b)), nrow = 2)
  tst = fisher.test(tab)
  odds = ((sum(a) + 0.5) * (length(b) - sum(b) + 0.5))
  odds = odds / ((length(a) - sum(a) + 0.5) * (sum(b) + 0.5))
  out = data.frame(feature = label, group_a = level_a, group_b = level_b,
                   n_a = length(a), mutated_a = sum(a), freq_a = mean(a),
                   n_b = length(b), mutated_b = sum(b), freq_b = mean(b),
                   odds_ratio = odds, log2_odds_ratio = log2(odds),
                   p_value = tst$p.value, stringsAsFactors = FALSE)
  out
}

# The five planned contrasts for one continuous feature.
contrasts_continuous = function(value, group, smoking, label) {
  cell = paste(group, smoking)
  out = NULL
  for (i in seq_along(CONTRAST_PAIRS)) {
    r = test_two_group(value, cell, CONTRAST_PAIRS[[i]][1],
                       CONTRAST_PAIRS[[i]][2], label)
    if (is.null(r)) {
      next
    }
    r$contrast = CONTRAST_LABELS[i]
    out = rbind(out, r)
  }
  out
}

# The five planned contrasts for one binary feature.
contrasts_binary = function(hit, group, smoking, label) {
  cell = paste(group, smoking)
  out = NULL
  for (i in seq_along(CONTRAST_PAIRS)) {
    r = test_two_group_binary(hit, cell, CONTRAST_PAIRS[[i]][1],
                              CONTRAST_PAIRS[[i]][2], label)
    if (is.null(r)) {
      next
    }
    r$contrast = CONTRAST_LABELS[i]
    out = rbind(out, r)
  }
  out
}

# Benjamini-Hochberg applied inside each contrast family.
adjust_within_contrast = function(df) {
  df$p_adjusted = NA
  families = unique(df$contrast)
  for (i in seq_along(families)) {
    rows = df$contrast == families[i]
    df$p_adjusted[rows] = p.adjust(df$p_value[rows], method = "BH")
  }
  df
}

# Difference in differences of medians, with a percentile bootstrap interval.
# (Western ever minus Western never) minus (EAS ever minus EAS never).
did_medians = function(value, group, smoking, label, n_boot = 2000) {
  keep = !is.na(value) & !is.na(group) & !is.na(smoking)
  value = value[keep]
  group = group[keep]
  smoking = smoking[keep]
  en = value[group == "EAS" & smoking == "Never"]
  ee = value[group == "EAS" & smoking == "Ever"]
  wn = value[group == "Western" & smoking == "Never"]
  we = value[group == "Western" & smoking == "Ever"]
  if (min(length(en), length(ee), length(wn), length(we)) < 3) {
    return(NULL)
  }
  delta_eas = median(ee) - median(en)
  delta_west = median(we) - median(wn)
  boot = rep(NA, n_boot)
  for (i in seq_len(n_boot)) {
    b_en = median(sample(en, length(en), replace = TRUE))
    b_ee = median(sample(ee, length(ee), replace = TRUE))
    b_wn = median(sample(wn, length(wn), replace = TRUE))
    b_we = median(sample(we, length(we), replace = TRUE))
    boot[i] = (b_we - b_wn) - (b_ee - b_en)
  }
  p_value = 2 * min(mean(boot <= 0), mean(boot >= 0))
  if (p_value > 1) {
    p_value = 1
  }
  if (p_value < 1 / (n_boot + 1)) {
    p_value = 1 / (n_boot + 1)
  }
  out = data.frame(feature = label,
                   n_eas_never = length(en), n_eas_ever = length(ee),
                   n_west_never = length(wn), n_west_ever = length(we),
                   median_eas_never = median(en), median_eas_ever = median(ee),
                   median_west_never = median(wn), median_west_ever = median(we),
                   delta_eas = delta_eas, delta_western = delta_west,
                   did = delta_west - delta_eas,
                   ci_low = unname(quantile(boot, 0.025)),
                   ci_high = unname(quantile(boot, 0.975)),
                   p_bootstrap = p_value, stringsAsFactors = FALSE)
  out
}

# Model-based interaction on the same 2x2, reported alongside the bootstrap.
did_models = function(value, group, smoking, cohort, label) {
  keep = !is.na(value) & !is.na(group) & !is.na(smoking)
  d = data.frame(y = value[keep], anc = factor(group[keep], c("EAS", "Western")),
                 smoke = factor(smoking[keep], c("Never", "Ever")),
                 cohort = factor(cohort[keep]), stringsAsFactors = FALSE)
  fit_ols = lm(y ~ anc * smoke, data = d)
  co = summary(fit_ols)$coefficients
  term = grep(":", rownames(co), value = TRUE)[1]
  d$rank_y = rank(d$y)
  fit_rank = lm(rank_y ~ anc * smoke, data = d)
  co_rank = summary(fit_rank)$coefficients
  out = data.frame(feature = label,
                   ols_beta = co[term, 1], ols_p = co[term, 4],
                   rank_p = co_rank[term, 4],
                   mixed_p = NA, stringsAsFactors = FALSE)
  if (length(unique(d$cohort)) > 2) {
    fit_mix = try(nlme::lme(y ~ anc * smoke, random = ~ 1 | cohort, data = d,
                            control = nlme::lmeControl(opt = "optim")),
                  silent = TRUE)
    if (!inherits(fit_mix, "try-error")) {
      tt = summary(fit_mix)$tTable
      term2 = grep(":", rownames(tt), value = TRUE)[1]
      out$mixed_p = tt[term2, "p-value"]
    }
  }
  out
}

# Change from never to ever smoking inside one cohort, with a bootstrap CI.
within_cohort_delta = function(value, smoking, label, n_boot = 2000) {
  keep = !is.na(value) & !is.na(smoking)
  value = value[keep]
  smoking = smoking[keep]
  nv = value[smoking == "Never"]
  ev = value[smoking == "Ever"]
  if (length(nv) < 3 || length(ev) < 3) {
    return(NULL)
  }
  boot = rep(NA, n_boot)
  for (i in seq_len(n_boot)) {
    boot[i] = median(sample(ev, length(ev), replace = TRUE)) -
              median(sample(nv, length(nv), replace = TRUE))
  }
  p_value = 2 * min(mean(boot <= 0), mean(boot >= 0))
  if (p_value > 1) {
    p_value = 1
  }
  if (p_value < 1 / (n_boot + 1)) {
    p_value = 1 / (n_boot + 1)
  }
  out = data.frame(feature = label, n_never = length(nv), n_ever = length(ev),
                   median_never = median(nv), median_ever = median(ev),
                   delta = median(ev) - median(nv),
                   ci_low = unname(quantile(boot, 0.025)),
                   ci_high = unname(quantile(boot, 0.975)),
                   p_bootstrap = p_value, stringsAsFactors = FALSE)
  out
}

# Difference in differences between two named cohorts.
pairwise_did = function(value, cohort, smoking, cohort_a, cohort_b, label,
                        n_boot = 2000) {
  an = value[cohort == cohort_a & smoking == "Never"]
  ae = value[cohort == cohort_a & smoking == "Ever"]
  bn = value[cohort == cohort_b & smoking == "Never"]
  be = value[cohort == cohort_b & smoking == "Ever"]
  an = an[!is.na(an)]
  ae = ae[!is.na(ae)]
  bn = bn[!is.na(bn)]
  be = be[!is.na(be)]
  if (min(length(an), length(ae), length(bn), length(be)) < 5) {
    return(NULL)
  }
  boot = rep(NA, n_boot)
  for (i in seq_len(n_boot)) {
    b_an = median(sample(an, length(an), replace = TRUE))
    b_ae = median(sample(ae, length(ae), replace = TRUE))
    b_bn = median(sample(bn, length(bn), replace = TRUE))
    b_be = median(sample(be, length(be), replace = TRUE))
    boot[i] = (b_be - b_bn) - (b_ae - b_an)
  }
  p_value = 2 * min(mean(boot <= 0), mean(boot >= 0))
  if (p_value > 1) {
    p_value = 1
  }
  if (p_value < 1 / (n_boot + 1)) {
    p_value = 1 / (n_boot + 1)
  }
  out = data.frame(feature = label, cohort_a = cohort_a, cohort_b = cohort_b,
                   n_a_never = length(an), n_a_ever = length(ae),
                   n_b_never = length(bn), n_b_ever = length(be),
                   delta_a = median(ae) - median(an),
                   delta_b = median(be) - median(bn),
                   did = (median(be) - median(bn)) - (median(ae) - median(an)),
                   ci_low = unname(quantile(boot, 0.025)),
                   ci_high = unname(quantile(boot, 0.975)),
                   p_bootstrap = p_value, stringsAsFactors = FALSE)
  out
}

# ---------------------------------------------------------------- figures ----

# Volcano used for every genome-wide or all-signature scan, so that the driver
# and signature volcanoes look like one family.
volcano_plot = function(df, x_col, p_col, label_col, title, x_lab,
                        label_set = NULL, n_label = 12) {
  d = df
  d$x = d[[x_col]]
  d$p = d[[p_col]]
  d$name = d[[label_col]]
  d = d[is.finite(d$x) & !is.na(d$p), ]
  d$y = -log10(pmax(d$p, 1e-300))
  d$direction = "ns"
  d$direction[d$p < 0.05 & d$x > 0] = "EAS-enriched"
  d$direction[d$p < 0.05 & d$x < 0] = "Western-enriched"
  sig = d[d$direction != "ns", ]
  up = sig[sig$x > 0, ]
  down = sig[sig$x < 0, ]
  up = up[order(-up$y), ]
  down = down[order(-down$y), ]
  lab = rbind(head(up, n_label), head(down, n_label))
  if (!is.null(label_set)) {
    forced = d[d$name %in% label_set, ]
    lab = rbind(lab, forced)
    lab = lab[!duplicated(lab$name), ]
  }
  p = ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = direction))
  p = p + ggplot2::geom_point(size = 0.8, alpha = 0.65, stroke = 0)
  p = p + ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                              linewidth = 0.3, colour = "grey40")
  p = p + ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40")
  p = p + ggrepel::geom_text_repel(data = lab,
            ggplot2::aes(label = name), size = 2.4, max.overlaps = 30,
            segment.size = 0.2, min.segment.length = 0, show.legend = FALSE)
  p = p + ggplot2::scale_colour_manual(values = DIR_COLS, name = NULL)
  p = p + ggplot2::labs(x = x_lab, y = expression(-log[10] * " adjusted P"),
                        title = title)
  p = p + theme_pub()
  p = p + ggplot2::theme(legend.position = "top")
  p
}

# Boxplot with every tumour shown as a jittered point.
box_jitter_plot = function(df, x_col, y_col, fill_col, palette, ylab, title,
                           logscale = FALSE) {
  d = df
  d$x = d[[x_col]]
  d$y = d[[y_col]]
  d$fill = d[[fill_col]]
  d = d[!is.na(d$y), ]
  p = ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, fill = fill))
  p = p + ggplot2::geom_jitter(width = 0.26, size = 0.45, alpha = 0.28,
                               colour = "grey25", stroke = 0, show.legend = FALSE)
  p = p + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6,
                                linewidth = 0.35, width = 0.55)
  p = p + ggplot2::scale_fill_manual(values = palette, name = NULL)
  p = p + ggplot2::labs(x = NULL, y = ylab, title = title)
  p = p + theme_pub()
  if (logscale) {
    p = p + ggplot2::scale_y_log10()
  }
  p
}

# Forest plot for any estimate that carries a confidence interval.
forest_plot = function(df, estimate_col, low_col, high_col, y_col, colour_col,
                       palette, xlab, title) {
  d = df
  d$est = d[[estimate_col]]
  d$lo = d[[low_col]]
  d$hi = d[[high_col]]
  d$y = d[[y_col]]
  d$col = d[[colour_col]]
  p = ggplot2::ggplot(d, ggplot2::aes(x = est, y = y, colour = col))
  p = p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                              linewidth = 0.3, colour = "grey40")
  p = p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
                                  height = 0.2, linewidth = 0.5)
  p = p + ggplot2::geom_point(size = 2.1)
  p = p + ggplot2::scale_colour_manual(values = palette, name = NULL)
  p = p + ggplot2::labs(x = xlab, y = NULL, title = title)
  p = p + theme_pub()
  p = p + ggplot2::theme(legend.position = "top")
  p
}
