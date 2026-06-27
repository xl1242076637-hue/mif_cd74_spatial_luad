#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scTenifoldKnk))

input_dir <- file.path("results", "tables", "sctenifoldknk_smoke")
output_dir <- file.path("results", "tables", "sctenifoldknk_original_expanded")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

runs <- data.frame(
  source_dataset = c(rep("GSE308103", 1), rep("GSE308103", 6)),
  broad_class = c("epithelial", rep("macrophage", 6)),
  target_gene = c("MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7"),
  input_csv = c(
    file.path(input_dir, "gse308103_epithelial_mif_counts_for_sctenifoldknk.csv"),
    rep(file.path(input_dir, "gse308103_macrophage_cd74_counts_for_sctenifoldknk.csv"), 6)
  ),
  stringsAsFactors = FALSE
)

safe_name <- function(text) {
  tolower(gsub("[^A-Za-z0-9]+", "_", text))
}

run_one <- function(source_dataset, broad_class, target_gene, input_csv) {
  counts <- read.csv(input_csv, row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  if (!(target_gene %in% rownames(counts))) {
    stop(sprintf("Target gene %s is absent from %s", target_gene, input_csv))
  }

  prefix <- paste(
    tolower(source_dataset),
    safe_name(broad_class),
    tolower(target_gene),
    "original_sctenifoldknk",
    sep = "_"
  )
  diffreg_csv <- file.path(output_dir, paste0(prefix, "_diffregulation.csv"))
  top_csv <- file.path(output_dir, paste0(prefix, "_top_non_target_genes.csv"))

  set.seed(7)
  result <- scTenifoldKnk(
    countMatrix = counts,
    qc = FALSE,
    gKO = target_gene,
    nc_nNet = 3,
    nc_nCells = min(75, ncol(counts)),
    nc_nComp = 2,
    td_K = 2,
    td_maxIter = 120,
    ma_nDim = 2,
    nCores = 1
  )

  diff_reg <- as.data.frame(result$diffRegulation)
  if (!("gene" %in% colnames(diff_reg))) {
    diff_reg <- cbind(gene = rownames(diff_reg), diff_reg)
  }
  rownames(diff_reg) <- NULL
  diff_reg <- cbind(
    source_dataset = source_dataset,
    broad_class = broad_class,
    target_gene = target_gene,
    n_input_genes = nrow(counts),
    n_input_cells = ncol(counts),
    run_mode = "original_scTenifoldKnk_expanded_sensitivity",
    diff_reg
  )
  write.csv(diff_reg, diffreg_csv, row.names = FALSE)

  ordered <- diff_reg
  if ("p.adj" %in% colnames(ordered)) {
    ordered <- ordered[order(ordered$p.adj, decreasing = FALSE), , drop = FALSE]
  }
  ordered_non_target <- ordered[ordered$gene != target_gene, , drop = FALSE]
  write.csv(head(ordered_non_target, 20), top_csv, row.names = FALSE)

  top_genes <- head(ordered_non_target$gene, 10)
  top_p <- if (nrow(ordered_non_target) > 0 && "p.adj" %in% colnames(ordered_non_target)) {
    as.numeric(ordered_non_target$p.adj[[1]])
  } else {
    NA_real_
  }

  data.frame(
    source_dataset = source_dataset,
    broad_class = broad_class,
    target_gene = target_gene,
    input_csv = input_csv,
    diffregulation_csv = diffreg_csv,
    top_non_target_csv = top_csv,
    n_input_genes = nrow(counts),
    n_input_cells = ncol(counts),
    n_diffregulation_rows = nrow(diff_reg),
    top_non_target_gene_by_p_adj = if (length(top_genes) > 0) top_genes[[1]] else NA_character_,
    top_non_target_p_adj = top_p,
    top_10_non_target_genes = paste(top_genes, collapse = ","),
    stringsAsFactors = FALSE
  )
}

summary_rows <- do.call(
  rbind,
  Map(run_one, runs$source_dataset, runs$broad_class, runs$target_gene, runs$input_csv)
)

summary_rows$r_version <- R.version.string
summary_rows$sctenifoldknk_version <- as.character(utils::packageVersion("scTenifoldKnk"))
summary_rows$seed <- 7
summary_rows$qc <- FALSE
summary_rows$nc_nNet <- 3
summary_rows$nc_nCells <- 75
summary_rows$nc_nComp <- 2
summary_rows$td_K <- 2
summary_rows$td_maxIter <- 120
summary_rows$ma_nDim <- 2
summary_rows$nCores <- 1

write.csv(
  summary_rows,
  file.path(output_dir, "sctenifoldknk_original_expanded_summary.csv"),
  row.names = FALSE
)

print(summary_rows)
