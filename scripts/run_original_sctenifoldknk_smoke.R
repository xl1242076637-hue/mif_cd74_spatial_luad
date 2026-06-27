#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scTenifoldKnk))

output_dir <- file.path("results", "tables", "sctenifoldknk_smoke")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

runs <- data.frame(
  source_dataset = c("GSE308103", "GSE308103"),
  broad_class = c("epithelial", "macrophage"),
  target_gene = c("MIF", "CD74"),
  input_csv = c(
    file.path(output_dir, "gse308103_epithelial_mif_counts_for_sctenifoldknk.csv"),
    file.path(output_dir, "gse308103_macrophage_cd74_counts_for_sctenifoldknk.csv")
  ),
  output_csv = c(
    file.path(output_dir, "gse308103_epithelial_mif_sctenifoldknk_diffregulation.csv"),
    file.path(output_dir, "gse308103_macrophage_cd74_sctenifoldknk_diffregulation.csv")
  ),
  stringsAsFactors = FALSE
)

run_one <- function(source_dataset, broad_class, target_gene, input_csv, output_csv) {
  counts <- read.csv(input_csv, row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  if (!(target_gene %in% rownames(counts))) {
    stop(sprintf("Target gene %s is absent from %s", target_gene, input_csv))
  }

  set.seed(7)
  result <- scTenifoldKnk(
    countMatrix = counts,
    qc = FALSE,
    gKO = target_gene,
    nc_nNet = 2,
    nc_nCells = min(60, ncol(counts)),
    nc_nComp = 2,
    td_K = 2,
    td_maxIter = 80,
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
    run_mode = "original_scTenifoldKnk_smoke",
    diff_reg
  )

  write.csv(diff_reg, output_csv, row.names = FALSE)

  top_gene <- NA_character_
  top_p_adj <- NA_real_
  top_non_target_gene <- NA_character_
  top_non_target_p_adj <- NA_real_
  top_non_target_output_csv <- sub(
    "_diffregulation[.]csv$",
    "_top_non_target_genes.csv",
    output_csv
  )
  if ("p.adj" %in% colnames(diff_reg)) {
    ordered <- diff_reg[order(diff_reg$p.adj, decreasing = FALSE), , drop = FALSE]
    if (nrow(ordered) > 0) {
      top_gene <- as.character(ordered$gene[[1]])
      top_p_adj <- as.numeric(ordered$p.adj[[1]])
    }

    ordered_non_target <- ordered[ordered$gene != target_gene, , drop = FALSE]
    if (nrow(ordered_non_target) > 0) {
      top_non_target_gene <- as.character(ordered_non_target$gene[[1]])
      top_non_target_p_adj <- as.numeric(ordered_non_target$p.adj[[1]])
    }
    write.csv(head(ordered_non_target, 20), top_non_target_output_csv, row.names = FALSE)
  }

  data.frame(
    source_dataset = source_dataset,
    broad_class = broad_class,
    target_gene = target_gene,
    input_csv = input_csv,
    output_csv = output_csv,
    n_input_genes = nrow(counts),
    n_input_cells = ncol(counts),
    n_diffregulation_rows = nrow(diff_reg),
    top_gene_by_p_adj = top_gene,
    top_p_adj = top_p_adj,
    top_non_target_gene_by_p_adj = top_non_target_gene,
    top_non_target_p_adj = top_non_target_p_adj,
    top_non_target_output_csv = top_non_target_output_csv,
    stringsAsFactors = FALSE
  )
}

summary_rows <- do.call(
  rbind,
  Map(
    run_one,
    runs$source_dataset,
    runs$broad_class,
    runs$target_gene,
    runs$input_csv,
    runs$output_csv
  )
)

summary_rows$r_version <- R.version.string
summary_rows$sctenifoldknk_version <- as.character(utils::packageVersion("scTenifoldKnk"))
summary_rows$seed <- 7
summary_rows$qc <- FALSE
summary_rows$nc_nNet <- 2
summary_rows$nc_nCells <- 60
summary_rows$nc_nComp <- 2
summary_rows$td_K <- 2
summary_rows$td_maxIter <- 80
summary_rows$ma_nDim <- 2
summary_rows$nCores <- 1

write.csv(
  summary_rows,
  file.path(output_dir, "sctenifoldknk_original_smoke_summary.csv"),
  row.names = FALSE
)

print(summary_rows)
