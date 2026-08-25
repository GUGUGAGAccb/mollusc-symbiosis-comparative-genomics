#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: script PROJECT STEP14_RUN STEP15_RUN OUT")
project <- args[[1]]
step14 <- args[[2]]
step15 <- args[[3]]
out <- args[[4]]
if (dir.exists(out)) stop("Output directory already exists: ", out)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

paths <- list(
  busco = file.path(project, "primary_proteome_busco_summary.tsv"),
  orthogroups = file.path(project, "results_or_reports/orthogroup_distribution/run_18269662/step04_summary.tsv"),
  candidates = file.path(project, "results_or_reports/candidate_gene_families/base_model_v1/run_18271359/step06_summary.tsv"),
  annotation = file.path(project, "results_or_reports/functional_annotation/combined_summary/run_18272016/step07_annotation_coverage_summary.tsv"),
  enrichment = file.path(project, "results_or_reports/functional_enrichment/orthogroup_level_v1/run_18272100/step08_summary.tsv"),
  validation = file.path(project, "results_or_reports/candidate_sequence_validation/combined_summary/run_18273530/step09_summary.tsv"),
  trees = file.path(project, "results_or_reports/candidate_family_phylogenetics/validated_top30_v1/run_18274798/step10_summary.tsv"),
  reconciliation = file.path(project, "results_or_reports/gene_tree_species_tree_reconciliation/combined_summary/dl_lca_v1/run_18286960/step11_summary.tsv"),
  evolution = file.path(project, "results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/run_18288232/step12_summary.tsv"),
  expression = file.path(project, "results_or_reports/transcriptomic_integration/single_library_salmon_v1/run_18420304/step13_summary.tsv"),
  priority_summary = file.path(step14, "step14_summary.tsv"),
  priority_families = file.path(step14, "ranked_candidate_families.tsv"),
  sensitivity_summary = file.path(step15, "step15_summary.tsv"),
  sensitivity_scenarios = file.path(step15, "scenario_comparison_summary.tsv"),
  stability = file.path(step15, "candidate_rank_stability.tsv")
)

busco <- read_tsv(paths$busco)
og <- read_tsv(paths$orthogroups)
cand <- read_tsv(paths$candidates)
ann <- read_tsv(paths$annotation)
enrich <- read_tsv(paths$enrichment)
valid <- read_tsv(paths$validation)
trees <- read_tsv(paths$trees)
recon <- read_tsv(paths$reconciliation)
evol <- read_tsv(paths$evolution)
expr <- read_tsv(paths$expression)
prio <- read_tsv(paths$priority_summary)
families <- read_tsv(paths$priority_families)
sens <- read_tsv(paths$sensitivity_summary)
scenarios <- read_tsv(paths$sensitivity_scenarios)
stability <- read_tsv(paths$stability)

bivalves <- c("Acanthocardia_echinata", "Americardia_media", "Cerastoderma_edule",
              "Fragum_sueziense", "Tridacna_gigas", "Tridacna_maxima")
gastropods <- c("Aplysia_californica", "Berghia_stephanieae", "Bullacta_exarata",
                "Elysia_chlorotica", "Elysia_crispata", "Elysia_marginata",
                "Onchidella_celtica", "Plakobranchus_ocellatus")
busco$group <- ifelse(busco$species %in% bivalves, "Bivalvia",
                      ifelse(busco$species %in% gastropods, "Gastropoda", "Supplementary"))
busco_core <- busco[busco$group != "Supplementary", ]
busco_core$species_label <- gsub("_", " ", busco_core$species)

palette <- c(Bivalvia = "#2C7FB8", Gastropoda = "#D95F0E")
theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "#F0F0F0"),
        plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot <- function(plot, stem, width = 8.5, height = 6.0) {
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
}

# Figure 1: proteome quality.
p1 <- ggplot(busco_core, aes(x = complete_pct, y = duplicated_pct, color = group,
                             size = proteins, label = species_label)) +
  geom_hline(yintercept = 10, linetype = 3, color = "grey60") +
  geom_vline(xintercept = 80, linetype = 3, color = "grey60") +
  geom_point(alpha = 0.85) +
  geom_text(size = 2.5, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  scale_color_manual(values = palette) + scale_size_continuous(range = c(3, 10)) +
  labs(title = "Proteome completeness and redundancy", x = "Complete BUSCOs (%)",
       y = "Duplicated BUSCOs (%)", color = "Class", size = "Proteins") + theme_pub
save_plot(p1, "figure01_proteome_quality", 9, 6.5)

# Figure 2: orthogroup composition.
og2 <- og[og$group %in% c("bivalvia", "gastropoda"), ]
og_long <- do.call(rbind, lapply(seq_len(nrow(og2)), function(i) {
  data.frame(group = tools::toTitleCase(og2$group[i]),
             category = c("Universal single-copy", "Universal multicopy", "Near-core", "Shell", "Species-specific"),
             count = as.numeric(og2[i, c("universal_single_copy", "universal_multicopy", "near_core", "shell", "species_specific")]))
}))
og_long$fraction <- ave(og_long$count, og_long$group, FUN = function(x) x / sum(x))
og_long$category <- factor(og_long$category,
  levels = c("Species-specific", "Shell", "Near-core", "Universal multicopy", "Universal single-copy"))
p2 <- ggplot(og_long, aes(x = group, y = fraction, fill = category)) + geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("#E34A33", "#FDBB84", "#FEE8C8", "#A6BDDB", "#2B8CBE")) +
  labs(title = "Orthogroup distribution by molluscan class", x = NULL,
       y = "Fraction of classified orthogroups", fill = "Category") + theme_pub
save_plot(p2, "figure02_orthogroup_distribution", 8.5, 6)

# Figure 3: annotation coverage.
ann$group <- tools::toTitleCase(ann$group)
ann_keep <- ann[ann$metric %in% c("seed_ortholog", "functional_annotation", "preferred_name",
                                  "description", "GO", "EC", "KEGG_KO", "Pfam"), ]
ann_keep$metric_label <- factor(ann_keep$metric,
  levels = rev(c("seed_ortholog", "functional_annotation", "description", "preferred_name", "GO", "KEGG_KO", "EC", "Pfam")),
  labels = rev(c("Seed ortholog", "Functional annotation", "Description", "Preferred name", "GO", "KEGG KO", "EC", "Pfam")))
p3 <- ggplot(ann_keep, aes(x = fraction, y = metric_label, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = palette) +
  labs(title = "Functional annotation coverage", x = "Annotated candidate proteins", y = NULL, fill = "Class") + theme_pub
save_plot(p3, "figure03_functional_annotation_coverage", 8.5, 6)

# Figure 4: validation outcome and high-risk flags.
valid$group <- tools::toTitleCase(valid$group)
val_long <- do.call(rbind, lapply(seq_len(nrow(valid)), function(i) {
  data.frame(group = valid$group[i], outcome = c("Pass", "Review", "Review high", "Exclude recommended"),
             count = as.numeric(valid[i, c("pass", "review", "review_high", "exclude_recommended")]),
             total = valid$total_sequences[i])
}))
val_long$fraction <- val_long$count / val_long$total
val_long$outcome <- factor(val_long$outcome, levels = c("Exclude recommended", "Review high", "Review", "Pass"))
p4 <- ggplot(val_long, aes(x = group, y = fraction, fill = outcome)) + geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("#B2182B", "#EF8A62", "#FDDBC7", "#4D9221")) +
  labs(title = "Candidate-sequence validation outcomes", x = NULL, y = "Fraction of sequences", fill = "Outcome") + theme_pub
save_plot(p4, "figure04_sequence_validation", 8.5, 6)

# Figure 5: evidence funnel from CAFE to robust shortlist.
funnel <- rbind(
  data.frame(group = tools::toTitleCase(cand$group), stage = "CAFE FDR candidates", count = cand$bh_fdr_le_0.05_families),
  data.frame(group = tools::toTitleCase(trees$group), stage = "Gene trees", count = trees$completed_trees),
  data.frame(group = tools::toTitleCase(evol$group), stage = "HyPhy eligible", count = evol$families_eligible),
  data.frame(group = tools::toTitleCase(prio$group), stage = "Tier 1 integrated", count = prio$tier1_families),
  data.frame(group = tools::toTitleCase(sens$group), stage = "Stable high", count = sens$stable_high_families)
)
funnel$stage <- factor(funnel$stage, levels = rev(c("CAFE FDR candidates", "Gene trees", "HyPhy eligible", "Tier 1 integrated", "Stable high")))
p5 <- ggplot(funnel, aes(x = count, y = stage, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000)) + scale_fill_manual(values = palette) +
  labs(title = "Candidate-evidence funnel", subtitle = "Logarithmic count scale",
       x = "Families", y = NULL, fill = "Class") + theme_pub
save_plot(p5, "figure05_candidate_evidence_funnel", 8.5, 6)

# Figure 6: reconciliation and molecular evolution evidence.
event_long <- rbind(
  data.frame(group = tools::toTitleCase(recon$group), metric = "Strong duplications", count = recon$strong_support_duplications),
  data.frame(group = tools::toTitleCase(recon$group), metric = "Inferred losses", count = recon$loss_events),
  data.frame(group = tools::toTitleCase(evol$group), metric = "BUSTED FDR families", count = evol$busted_fdr_significant),
  data.frame(group = tools::toTitleCase(evol$group), metric = "aBSREL FDR branches", count = evol$absrel_global_fdr_significant_branches),
  data.frame(group = tools::toTitleCase(evol$group), metric = "RELAX FDR families", count = evol$relax_fdr_significant)
)
event_long$metric <- factor(event_long$metric, levels = rev(c("Strong duplications", "Inferred losses", "BUSTED FDR families", "aBSREL FDR branches", "RELAX FDR families")))
p6 <- ggplot(event_long, aes(x = count, y = metric, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) + scale_fill_manual(values = palette) +
  labs(title = "Phylogenetic and molecular-evolution signals", x = "Events or significant tests", y = NULL, fill = "Class") + theme_pub
save_plot(p6, "figure06_phylogenetic_evolution_signals", 8.5, 6)

# Figure 7: top integrated families.
families$group_label <- tools::toTitleCase(families$group)
top <- do.call(rbind, lapply(split(families, families$group), function(x) x[order(x$group_rank_integrated), ][1:10, ]))
top$label <- paste0(top$orthogroup, "  ", ifelse(is.na(top$consensus_preferred_name) | top$consensus_preferred_name == "", "unannotated", top$consensus_preferred_name))
top$label_id <- paste(top$group_label, top$label, sep = "::")
top$label_id <- factor(top$label_id, levels = rev(unique(top$label_id[order(top$group_label, top$group_rank_integrated)])))
p7 <- ggplot(top, aes(x = integrated_score_normalized, y = label_id, color = group_label, shape = priority_tier)) +
  geom_point(size = 3.2) + facet_wrap(~group_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = palette) + scale_y_discrete(labels = function(x) sub("^[^:]+::", "", x)) +
  labs(title = "Top integrated candidate families", x = "Evidence-normalised score", y = NULL,
       color = "Class", shape = "Priority tier") + theme_pub + theme(legend.position = "bottom")
save_plot(p7, "figure07_top_integrated_families", 11, 7)

# Figure 8: sensitivity heatmap.
sc <- scenarios[scenarios$scenario != "baseline", ]
sc$group <- tools::toTitleCase(sc$group)
sc$scenario_label <- gsub("_", " ", sc$scenario)
scenario_order <- unique(sc$scenario_label)
sc$scenario_label <- factor(sc$scenario_label, levels = rev(scenario_order))
sc$label <- paste0(sprintf("%.2f", sc$spearman_rho_vs_baseline), " / ", sc$top10_overlap)
p8 <- ggplot(sc, aes(x = group, y = scenario_label, fill = spearman_rho_vs_baseline)) +
  geom_tile(color = "white", linewidth = 0.5) + geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "#D73027", mid = "#FEE08B", high = "#1A9850", midpoint = 0.75,
                       limits = c(0, 1), name = "Spearman rho") +
  labs(title = "Candidate-ranking sensitivity", subtitle = "Cell label: Spearman rho / baseline Top-10 overlap",
       x = NULL, y = NULL) + theme_pub + theme(axis.text.y = element_text(size = 9))
save_plot(p8, "figure08_ranking_sensitivity", 9.5, 7.5)

# Figure 9: transcriptomic detection thresholds.
expr$group <- tools::toTitleCase(expr$group)
expr_long <- rbind(
  data.frame(group = expr$group, threshold = "TPM >= 0.1", count = expr$candidate_detected_tpm_ge_0.1, total = expr$candidate_transcripts),
  data.frame(group = expr$group, threshold = "TPM >= 1", count = expr$candidate_detected_tpm_ge_1, total = expr$candidate_transcripts),
  data.frame(group = expr$group, threshold = "TPM >= 10", count = expr$candidate_detected_tpm_ge_10, total = expr$candidate_transcripts)
)
expr_long$fraction <- expr_long$count / expr_long$total
expr_long$threshold <- factor(expr_long$threshold, levels = c("TPM >= 0.1", "TPM >= 1", "TPM >= 10"))
p9 <- ggplot(expr_long, aes(x = threshold, y = fraction, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = palette) +
  labs(title = "Candidate transcript detection", subtitle = "Single selected library per represented species; not differential expression",
       x = NULL, y = "Detected candidate transcripts", fill = "Class") + theme_pub
save_plot(p9, "figure09_transcriptomic_detection", 8.5, 6)

# Multi-page PDF booklet.
pdf(file.path(out, "step16_figure_booklet.pdf"), width = 10, height = 7.5, onefile = TRUE)
for (p in list(p1, p2, p3, p4, p5, p6, p7, p8, p9)) print(p)
dev.off()

# Statistical tests: selected contrasts only; biological units and caveats are explicit.
test_rows <- list()
add_test <- function(section, metric, method, statistic, p, effect, unit, note) {
  test_rows[[length(test_rows) + 1]] <<- data.frame(section, metric, method,
    statistic = as.character(statistic), p_value = as.numeric(p), effect = as.character(effect),
    biological_unit = unit, interpretation_note = note, stringsAsFactors = FALSE)
}

for (metric in c("complete_pct", "duplicated_pct", "proteins")) {
  x <- busco_core[busco_core$group == "Bivalvia", metric]
  y <- busco_core[busco_core$group == "Gastropoda", metric]
  wt <- wilcox.test(x, y, exact = FALSE, conf.int = FALSE)
  add_test("proteome_quality", metric, "Wilcoxon rank-sum", unname(wt$statistic), wt$p.value,
           paste0("median_B=", median(x), ";median_G=", median(y)), "species",
           "Exploratory two-class contrast; small and taxonomically non-independent samples.")
}
ct <- cor.test(busco_core$complete_pct, log10(busco_core$proteins), method = "spearman", exact = FALSE)
add_test("proteome_quality", "BUSCO completeness vs log10 protein count", "Spearman correlation",
         unname(ct$estimate), ct$p.value, paste0("rho=", round(unname(ct$estimate), 4)), "species",
         "Descriptive association; does not establish annotation causality.")

og_b <- og2[og2$group == "bivalvia", ]; og_g <- og2[og2$group == "gastropoda", ]
tab <- matrix(c(og_b$assigned_proteins, og_b$proteins - og_b$assigned_proteins,
                og_g$assigned_proteins, og_g$proteins - og_g$assigned_proteins), nrow = 2, byrow = TRUE)
ft <- fisher.test(tab)
add_test("orthogroup_assignment", "assigned vs unassigned proteins", "Fisher exact", unname(ft$estimate), ft$p.value,
         paste0("OR=", signif(unname(ft$estimate), 5)), "protein",
         "Protein-level counts are nested within species; p-value is descriptive and may be anti-conservative.")

for (i in seq_len(nrow(ann_keep) / 2)) {
  metric <- unique(ann_keep$metric)[i]
  z <- ann_keep[ann_keep$metric == metric, ]
  if (nrow(z) == 2) {
    tab <- cbind(z$count, z$total_queries - z$count)
    ft <- fisher.test(tab)
    add_test("functional_annotation", metric, "Fisher exact", unname(ft$estimate), ft$p.value,
             paste0("B=", signif(z$fraction[z$group == "Bivalvia"], 4), ";G=", signif(z$fraction[z$group == "Gastropoda"], 4)),
             "protein", "Protein-level counts are nested within families and species; descriptive comparison.")
  }
}

for (group in c("Bivalvia", "Gastropoda")) {
  vv <- valid[valid$group == group, ]
  add_test("sequence_validation", paste0(group, " pass fraction"), "Binomial Wilson-compatible descriptive proportion",
           vv$pass, NA, paste0(vv$pass, "/", vv$total_sequences, "=", signif(vv$pass/vv$total_sequences, 4)),
           "protein", "No independence claim across homologous proteins.")
}

stats <- do.call(rbind, test_rows)
stats$p_adjust_bh <- NA_real_
valid_p <- !is.na(stats$p_value)
stats$p_adjust_bh[valid_p] <- p.adjust(stats$p_value[valid_p], method = "BH")
write_tsv(stats, "statistical_tests.tsv")

# Unified tables for downstream reporting.
write_tsv(og_long, "orthogroup_composition_long.tsv")
write_tsv(funnel, "candidate_evidence_funnel.tsv")
write_tsv(event_long, "phylogenetic_evolution_summary.tsv")
write_tsv(top[, c("global_rank", "group", "group_rank_integrated", "orthogroup", "priority_tier",
                  "integrated_score_normalized", "evidence_domain_count", "critical_risk_flag",
                  "consensus_preferred_name", "consensus_description", "best_event_node",
                  "best_event_direction", "best_event_change")], "top20_integrated_candidates.tsv")

workflow <- data.frame(
  step = 1:16,
  analysis = c("Proteome quality", "Orthogroup inference", "Species tree", "Orthogroup distribution",
               "Family expansion/contraction", "Candidate-family identification", "Functional annotation",
               "Functional enrichment", "Candidate sequence validation", "Candidate-family phylogenetics",
               "Gene-tree/species-tree reconciliation", "Molecular evolution", "Transcriptomic integration",
               "Integrated prioritisation", "Sensitivity and bias", "Statistics and visualisation"),
  status = rep("COMPLETE", 16),
  note = c("BUSCO and redundancy metrics", "OrthoFinder", "IQ-TREE class analyses and combined baseline",
           "Conserved/shell/species-specific summaries", "CAFE5 exploratory relative-time screen",
           "CAFE Base-model candidate extraction", "eggNOG-mapper", "Candidate-internal orthogroup enrichment",
           "Pfam/HMMER and risk flags", "60 supported gene trees", "Duplication-loss LCA reconciliation",
           "BUSTED, aBSREL and RELAX", "Six single-library Salmon datasets; no DE", "Evidence-weighted family/gene ranking",
           "13 ranking scenarios", "Unified tables, selected statistical contrasts and figures"),
  stringsAsFactors = FALSE
)
write_tsv(workflow, "pipeline_status.tsv")

source_manifest <- data.frame(source = names(paths), path = unlist(paths), stringsAsFactors = FALSE)
write_tsv(source_manifest, "source_manifest.tsv")

summary_lines <- c(
  "# Step 16 statistical analysis and visualisation", "",
  "## Scope", "",
  "This report consolidates completed Steps 1-15 into selected statistical contrasts, unified summary tables and publication-ready figures. It does not re-run upstream inference.", "",
  "## Principal descriptive results", "",
  sprintf("- Core dataset: %d Bivalvia and %d Gastropoda proteomes.", sum(busco_core$group == "Bivalvia"), sum(busco_core$group == "Gastropoda")),
  sprintf("- Orthogroup assignment rates: %.1f%% Bivalvia and %.1f%% Gastropoda.", og_b$assignment_rate_percent, og_g$assignment_rate_percent),
  sprintf("- CAFE/Step 6 FDR candidates: %d Bivalvia and %d Gastropoda families.", cand$bh_fdr_le_0.05_families[cand$group == "bivalvia"], cand$bh_fdr_le_0.05_families[cand$group == "gastropoda"]),
  sprintf("- Integrated Tier 1 families: %d Bivalvia and %d Gastropoda.", prio$tier1_families[prio$group == "bivalvia"], prio$tier1_families[prio$group == "gastropoda"]),
  sprintf("- Stable-high families across 13 scenarios: %d Bivalvia and %d Gastropoda.", sens$stable_high_families[sens$group == "bivalvia"], sens$stable_high_families[sens$group == "gastropoda"]),
  "", "## Statistical interpretation", "",
  "Species-level Wilcoxon tests are exploratory because each class has few species and observations are phylogenetically non-independent. Protein- and family-level Fisher tests use nested counts and therefore provide descriptive evidence rather than fully independent biological replication. BH-adjusted values are supplied for the selected formal contrasts.",
  "", "## Transcriptomic limitation", "",
  "Step 13 has one selected RNA-seq library per represented species and no biological replicates. Figure 9 reports detection thresholds only; it is not differential-expression evidence, and TPM is not compared quantitatively across species.",
  "", "## CAFE limitation", "",
  "CAFE5 analyses use relative root age 100 as exploratory screens. They do not supply absolute divergence-time conclusions. Gastropoda Gamma-model significance was unreliable; candidate work used the Base-model results and sensitivity analysis.",
  "", "## Files", "",
  "Individual figures are supplied as 300-dpi PNG files, with all figures also in a multi-page PDF booklet. Data behind principal panels are supplied as TSV files."
)
writeLines(summary_lines, file.path(out, "STEP16_REPORT.md"), useBytes = TRUE)

captions <- data.frame(
  figure = sprintf("Figure %d", 1:9),
  file = c("figure01_proteome_quality.png", "figure02_orthogroup_distribution.png",
           "figure03_functional_annotation_coverage.png", "figure04_sequence_validation.png",
           "figure05_candidate_evidence_funnel.png", "figure06_phylogenetic_evolution_signals.png",
           "figure07_top_integrated_families.png", "figure08_ranking_sensitivity.png",
           "figure09_transcriptomic_detection.png"),
  caption = c(
    "BUSCO completeness versus duplication for the 14 core proteomes; point area represents primary-proteome size.",
    "Relative composition of conserved, near-core, shell and species-specific orthogroups in class-specific OrthoFinder analyses.",
    "eggNOG-derived annotation coverage among Step 6 high-confidence candidate proteins.",
    "Sequence-validation categories from direct Pfam/HMMER scanning and rule-based quality assessment.",
    "Reduction from all CAFE FDR candidates to 30 gene trees per class, HyPhy-eligible families, integrated Tier 1 and cross-scenario stable-high candidates.",
    "Duplication/loss events and significant molecular-evolution results. Counts have different biological units and should not be compared as rates.",
    "Top ten Step 14 integrated families per class, including evidence-normalised score and priority tier.",
    "Sensitivity of class-specific rankings to 12 alternative scenarios. Labels show Spearman correlation and Top-10 overlap relative to corrected Step 14 baseline.",
    "Candidate-transcript detection at three within-library TPM thresholds. No differential-expression inference was performed."
  ), stringsAsFactors = FALSE
)
write_tsv(captions, "figure_captions.tsv")

capture.output(sessionInfo(), file = file.path(out, "R_session_info.txt"))

files_now <- list.files(out, full.names = TRUE)
manifest <- data.frame(file = basename(files_now), bytes = file.info(files_now)$size,
                       stringsAsFactors = FALSE)
write_tsv(manifest[order(manifest$file), ], "result_manifest.tsv")
writeLines("Step 16 statistical analysis and visualisation complete", file.path(out, "COMPLETE"))
cat("Step 16 completed:", out, "\n")
