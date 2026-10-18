# single cell analysis

# Parameters
WORK_DIR <- Sys.getenv("STUDY_WORK_DIR", unset = "~/rawData")
OUTPUT_DIR <- file.path(WORK_DIR, "publication_outputs", "hub_v3")

study.samples <- c("133_525_als", "147_521_als", "189_control", "888_control")
int.nfeatures <- 3000
int.npcs <- 50
int.chosen_dims <- 1:30
int.max_iter <- 30
int.resolutions <- seq(0.2, 1.0, by = 0.2)
int.final_res <- "SCT_snn_res.0.4"
de.test_use <- "wilcox"
de.logfc_threshold <- 0.25
de.sig_padj <- 0.05
de.sig_logfc <- 0.5
de.min_pct <- 0.1
met.bsg_subset <- "Activated (Late/Core)"
met.glycolysis_core <- c(
  "SLC2A1", "SLC2A3", "HK2", "PFKP", "PFKFB3", "GAPDH", "PGK1", "ENO1",
  "PKM", "LDHA", "SLC16A3", "PDK1"
)
en.kegg_organism <- "hsa"
en.kegg_pv_cutoff <- 0.05
en.gsea_pv_cutoff <- 0.05

# Helpers
group_from_sample <- function(x) {
  ifelse(grepl("525", x), "ALS_525",
         ifelse(grepl("521", x), "ALS_521", "Control"))
}
state_order <- c(
  "Homeostatic",
  "Activated (Early/Transitional)",
  "Activated (Late/Core)",
  "Activated (Lipid/Phagocytic)",
  "Activated (Interferon)",
  "Activated (Terminal/Stress)",
  "Proliferative",
  "Myeloid progenitors"
)

.hub_write_matrix <- function(x, file, id = "gene") {
  tab <- data.frame(id_value = rownames(x), as.data.frame(x, check.names = FALSE),
                    check.names = FALSE, row.names = NULL)
  names(tab)[1] <- id
  write.csv(tab, file, row.names = FALSE)
}

.hub_read_seurat <- function(path) {
  .hub_require_string(path, "input_rds")
  if (!file.exists(path) || dir.exists(path)) stop("Required input file is missing: ", path,
                             ". Supply input_rds explicitly if stored elsewhere.")
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Expected a cell-level Seurat object: ", path)
  cells <- colnames(obj)
  if (!nrow(obj) || !ncol(obj)) stop("The input object has no features or cells.")
  if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
    stop("Cell names must be nonempty and unique.")
  }
  if (!identical(rownames(obj@meta.data), cells)) {
    stop("Metadata rows must match the expression cell names in order.")
  }
  obj
}

.hub_object_record <- function(obj, file = "input_object_summary.csv") {
  write.csv(data.frame(
    cells = ncol(obj), features = nrow(obj),
    default_assay = DefaultAssay(obj),
    assays = paste(Assays(obj), collapse = ";"),
    reductions = paste(Reductions(obj), collapse = ";")
  ), file, row.names = FALSE)
}

.hub_require_states <- function(obj) {
  missing_cols <- setdiff(c("orig.ident", "Group", "Trajectory_States"),
                          colnames(obj@meta.data))
  if (length(missing_cols)) stop("Run annotation first; missing metadata: ",
                                  paste(missing_cols, collapse = ", "))
  .hub_require_study_samples(obj)
  states <- as.character(obj$Trajectory_States)
  if (anyNA(states) || any(!states %in% state_order)) {
    stop("Unassigned or unknown states; verify the original cluster-to-state map.")
  }
  groups <- as.character(obj$Group)
  if (anyNA(groups) || any(groups != group_from_sample(as.character(obj$orig.ident)))) {
    stop("Group labels must match the study's sample-to-group mapping. Run annotation first.")
  }
  invisible(obj)
}

.hub_require_study_samples <- function(obj) {
  if (!"orig.ident" %in% colnames(obj@meta.data)) stop("Required metadata orig.ident is absent.")
  samples <- as.character(obj$orig.ident)
  if (anyNA(samples) || any(!nzchar(samples))) stop("orig.ident contains missing/empty sample labels.")
  unknown <- setdiff(unique(samples), study.samples)
  if (length(unknown)) stop("Unexpected sample IDs for this study: ", paste(unknown, collapse = ", "),
                            ". The original grouping rules apply only to study.samples.")
  invisible(obj)
}

.hub_process_de <- function(de_df, comparison_label) {
  de_df$gene        <- rownames(de_df)
  de_df$p_val_adj   <- p.adjust(de_df$p_val, method = "BH")
  de_df$direction   <- ifelse(de_df$avg_log2FC > 0, "+", "-")
  de_df$comparison  <- comparison_label
  de_df$cell_type   <- "ALL"
  de_df$significance <- ifelse(abs(de_df$avg_log2FC) > de.sig_logfc &
                                 de_df$p_val_adj < de.sig_padj,
                               "logFC > 0.5 & FDR < 0.05", "ns")
  de_df[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj",
            "cell_type", "comparison", "direction", "significance")]
}

.hub_classify_de <- function(de, label) {
  data.frame(gene = rownames(de), comparison = label,
             avg_log2FC = de$avg_log2FC, p_val_adj_source = de$p_val_adj,
             change = dplyr::case_when(
               de$avg_log2FC > 0.25 & de$p_val_adj < 0.05 ~ "Up",
               de$avg_log2FC < -0.25 & de$p_val_adj < 0.05 ~ "Down",
               TRUE ~ "Stable"))
}

.hub_require_string <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(label, " must be one nonempty string.", call. = FALSE)
  }
  invisible(x)
}

.hub_require_groups <- function(x, required, label) {
  missing_groups <- setdiff(required, unique(as.character(x)))
  if (length(missing_groups)) {
    stop("Required ", label, " groups are absent: ", paste(missing_groups, collapse = ", "),
         call. = FALSE)
  }
  invisible(x)
}

.hub_require_assay <- function(obj, assay) {
  if (!assay %in% Assays(obj)) stop("Required assay is absent: ", assay, call. = FALSE)
  invisible(obj)
}

.hub_require_feature <- function(obj, feature) {
  if (!feature %in% rownames(obj)) {
    stop("Required feature is absent from the active assay: ", feature,
         ". Verify the original assay and feature identifiers.", call. = FALSE)
  }
  invisible(obj)
}

.hub_require_cell_alignment <- function(values, obj) {
  if (!identical(rownames(values), colnames(obj))) {
    stop("Fetched expression rows do not match the input cells in order. ",
         "Check missing expression/metadata and the original assay; no cells were relabeled.",
         call. = FALSE)
  }
  invisible(values)
}

# Run one stage
run_hub <- function(stage, work_dir = WORK_DIR,
                    output_dir = file.path(work_dir, "publication_outputs", "hub_v3"),
                    input_rds = NULL) {
  allowed <- c("integration", "annotation", "de", "bsg", "enrichment", "metabolism")
  if (missing(stage) || !is.character(stage) || length(stage) != 1L ||
      is.na(stage) || !stage %in% allowed) {
    stop("Choose one stage: ", paste(allowed, collapse = ", "))
  }
  .hub_require_string(work_dir, "work_dir")
  .hub_require_string(output_dir, "output_dir")
  work_dir <- normalizePath(path.expand(work_dir), mustWork = TRUE)
  if (!dir.exists(work_dir)) stop("work_dir must be an existing directory.")
  output_dir <- path.expand(output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  if (is.null(input_rds)) {
    input_rds <- switch(stage,
      integration = file.path(work_dir, "ALS_Microglia_QC_Annotated.rds"),
      annotation = file.path(work_dir, "sc_hd_wgcna_completed.rds"),
      de = file.path(work_dir, "sc_hd_wgcna_completed.rds"),
      file.path(output_dir, "annotation", "ALS_Microglia_State_Annotated.rds"))
  }
  if (!is.null(input_rds)) {
    .hub_require_string(input_rds, "input_rds")
    input_rds <- path.expand(input_rds)
    if (!file.exists(input_rds) || dir.exists(input_rds)) stop("Required input is missing: ", input_rds, ". Supply input_rds explicitly if stored elsewhere.")
    input_rds <- normalizePath(input_rds, mustWork = TRUE)
  }
  deps <- switch(stage,
    integration = c("Seurat", "harmony", "future"),
    annotation = c("Seurat", "dplyr", "Matrix"),
    de = c("Seurat", "dplyr"),
    bsg = c("Seurat", "dplyr"),
    enrichment = c("Seurat", "dplyr", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "msigdbr"),
    metabolism = c("Seurat", "dplyr"))
  deps <- unique(c(deps, "future"))
  missing_pkgs <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install the documented dependencies: ", paste(missing_pkgs, collapse = ", "))
  suppressPackageStartupMessages(for (pkg in deps) library(pkg, character.only = TRUE))
  stage_dir <- file.path(output_dir, stage)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  previous_wd <- getwd()
  setwd(stage_dir)
  completed <- FALSE
  on.exit({
    tryCatch({
      writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
      saveRDS(if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL,
              "rng_state_after.rds")
      write.csv(data.frame(stage = stage, completed = completed,
                           time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                "stage_status.csv", row.names = FALSE)
    }, finally = setwd(previous_wd))
  }, add = TRUE)
  saveRDS(if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL,
          "rng_state_before.rds")
  write.csv(data.frame(stage = stage, work_dir = work_dir,
                       input_rds = if (is.null(input_rds)) NA_character_ else input_rds),
            "run_input.csv", row.names = FALSE)
  parameter_names <- c("study.samples", "int.nfeatures",
                       "int.npcs", "int.chosen_dims", "int.max_iter", "int.resolutions", "int.final_res",
                       "de.test_use", "de.logfc_threshold", "de.sig_padj", "de.sig_logfc", "de.min_pct",
                       "met.bsg_subset", "met.glycolysis_core", "en.kegg_organism", "en.kegg_pv_cutoff", "en.gsea_pv_cutoff")
  parameters_used <- mget(parameter_names, envir = environment(run_hub), inherits = TRUE)
  saveRDS(parameters_used, "analysis_parameters.rds")
  dput(parameters_used, file = "analysis_parameters.txt")
  future::plan("sequential")
  options(future.globals.maxSize = 500 * 1024^2)
  fun <- get(paste0(".hub_", stage), mode = "function")
  fun(work_dir = work_dir, input_rds = input_rds)
  completed <- TRUE
  invisible(stage_dir)
}

# Integration
.hub_integration <- function(work_dir, input_rds) {
  sc_data <- .hub_read_seurat(input_rds)
  .hub_require_study_samples(sc_data)
  .hub_require_assay(sc_data, "RNA")
  .hub_object_record(sc_data)
  sc_data <- NormalizeData(sc_data, assay = "RNA", verbose = FALSE)
  sc_data <- SCTransform(sc_data, assay = "RNA", vst.flavor = "v2", verbose = TRUE)

  sc_list <- SplitObject(sc_data, split.by = "orig.ident")
  integ_features <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = int.nfeatures)
  rm(sc_list); gc()

  sc_data <- PrepSCTFindMarkers(sc_data)
  VariableFeatures(sc_data, assay = "SCT") <- integ_features

  sc_data <- RunPCA(sc_data, assay = "SCT",
                    features = VariableFeatures(sc_data, assay = "SCT"),
                    npcs = int.npcs, verbose = FALSE)
  chosen_dims <- int.chosen_dims

  sc_data <- RunHarmony(sc_data,
                        dims.use         = chosen_dims,
                        group.by.vars    = "orig.ident",
                        reduction        = "pca",
                        reduction.save   = "harmony",
                        max.iter.harmony = int.max_iter,
                        verbose          = FALSE)

  sc_data <- FindNeighbors(sc_data, reduction = "harmony", dims = chosen_dims, verbose = FALSE)
  sc_data <- FindClusters(sc_data, resolution = int.resolutions, verbose = FALSE)
  final_res <- int.final_res
  Idents(sc_data) <- final_res
  sc_data <- RunUMAP(sc_data, reduction = "harmony", dims = chosen_dims,
                     reduction.name = "umap", verbose = FALSE)
  if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
  saveRDS(sc_data, "data/processed/ALS_Microglia_Harmony_Clustered.rds")

  for (reduction in c("pca", "harmony", "umap")) {
    .hub_write_matrix(Embeddings(sc_data, reduction),
                      paste0(reduction, "_embeddings.csv"), id = "cell")
  }
  .hub_write_matrix(sc_data@meta.data, "cluster_metadata.csv", id = "cell")
  write.csv(data.frame(PC = seq_along(Stdev(sc_data, reduction = "pca")),
                       standard_deviation = Stdev(sc_data, reduction = "pca")),
            "pca_standard_deviations.csv", row.names = FALSE)
  write.csv(data.frame(feature = integ_features), "integration_features.csv", row.names = FALSE)
  .hub_object_record(sc_data, "output_object_summary.csv")
}

# Cell annotation
.hub_annotation <- function(work_dir, input_rds) {
  sc_hd <- .hub_read_seurat(input_rds)
  .hub_require_study_samples(sc_hd)
  .hub_object_record(sc_hd)
  cluster_col <- intersect(c("RNA_snn_res.0.7", int.final_res, "seurat_clusters"),
                           colnames(sc_hd@meta.data))[1]
  if (is.na(cluster_col)) {
    stop("No supported cluster column found. Expected RNA_snn_res.0.7, ",
         int.final_res, ", or seurat_clusters.")
  }
  sc_hd$Cluster_final <- factor(sc_hd[[cluster_col, drop = TRUE]])

  trajectory_ids <- c(
    "2"  = "Homeostatic",
    "3"  = "Homeostatic",
    "1"  = "Activated (Early/Transitional)",
    "0"  = "Activated (Late/Core)",
    "4"  = "Activated (Lipid/Phagocytic)",
    "10" = "Activated (Interferon)",
    "9"  = "Activated (Terminal/Stress)",
    "5"  = "Proliferative",
    "7"  = "Proliferative",
    "8"  = "Proliferative",
    "11" = "Proliferative",
    "6"  = "Myeloid progenitors"
  )

  unmapped <- setdiff(unique(as.character(sc_hd$Cluster_final)), names(trajectory_ids))
  if (length(unmapped)) stop("Clusters absent from the author's fixed map: ",
                             paste(unmapped, collapse = ", "),
                             ". Confirm the original clustering; do not guess state labels.")
  write.csv(data.frame(cluster_column = cluster_col, default_assay = DefaultAssay(sc_hd),
                       mapping_origin = "Original HUB script, clusters 0-11"),
            "annotation_selection.csv", row.names = FALSE)
  write.csv(data.frame(cluster = names(trajectory_ids), state = unname(trajectory_ids)),
            "cluster_to_state_map.csv", row.names = FALSE)
  Idents(sc_hd) <- "Cluster_final"
  sc_hd <- RenameIdents(sc_hd, trajectory_ids)
  sc_hd$Trajectory_States <- Idents(sc_hd)
  sc_hd$Trajectory_States <- factor(sc_hd$Trajectory_States, levels = state_order)
  Idents(sc_hd) <- "Trajectory_States"
  trajectory_markers <- list(
    "Homeostatic"          = c("ADGRG1", "PREX2", "CDH20", "CLEC10A"),
    "Activated_Early"      = c("G0S2", "EDNRA", "PGA5", "KIF26B"),
    "Activated_Late_Core"  = c("MRC1", "SELENOP", "CD163L1"),
    "Activated_Lipid_DAM"  = c("FUCA1", "STAB1", "C15orf48"),
    "Activated_Interferon" = c("IFIT1", "IFIT3", "MX2", "IFI44L"),
    "Activated_Terminal"   = c("SAMD4A", "SULT1C2", "KCNMA1", "APBA1"),
    "Proliferative"        = c("MKI67", "TOP2A", "DTL", "CDC20", "H1-5"),
    "Myeloid_Progenitors"  = c("FCN1", "MCEMP1", "CHRM3")
  )
  valid_markers <- lapply(trajectory_markers, function(x) intersect(x, rownames(sc_hd)))
  sc_hd$Group    <- factor(group_from_sample(sc_hd$orig.ident),
                           levels = c("Control", "ALS_521", "ALS_525"))
  sc_hd$SampleID <- factor(sc_hd$orig.ident,
                           levels = c("189_control", "888_control",
                                      "147_521_als", "133_525_als"))
  meta <- sc_hd@meta.data
  meta$Trajectory_States <- factor(meta$Trajectory_States, levels = state_order)
  meta$Group    <- factor(meta$Group, levels = c("Control", "ALS_521", "ALS_525"))
  meta$SampleID <- factor(meta$SampleID,
                          levels = c("189_control", "888_control",
                                     "147_521_als", "133_525_als"))

  cell_counts <- meta %>%
    group_by(Group, SampleID, Trajectory_States) %>%
    summarise(count = n(), .groups = "drop")
  cell_prop <- cell_counts %>%
    group_by(Group) %>%
    mutate(proportion = count / sum(count))

  # Divide sample-state counts by the total cell count within each group.
  write.csv(cell_counts, "cell_counts_by_sample_state.csv", row.names = FALSE)
  write.csv(cell_prop, "cell_proportions_group_denominator.csv", row.names = FALSE)
  cell_group_totals <- cell_counts %>%
    group_by(Group, Trajectory_States) %>%
    summarise(count = sum(count), .groups = "drop")
  write.csv(cell_group_totals, "cell_counts_by_group_state.csv", row.names = FALSE)

  marker_panel <- do.call(rbind, lapply(names(trajectory_markers), function(panel) {
    data.frame(panel = panel, gene = trajectory_markers[[panel]],
               present = trajectory_markers[[panel]] %in% valid_markers[[panel]])
  }))
  write.csv(marker_panel, "state_marker_panel.csv", row.names = FALSE)
  features <- unique(unlist(valid_markers, use.names = FALSE))
  if (length(features)) {
    marker_values <- FetchData(sc_hd, vars = features)
    marker_states <- as.character(sc_hd$Trajectory_States[match(rownames(marker_values), colnames(sc_hd))])
    marker_stats <- list()
    for (state in state_order) {
      idx <- which(marker_states == state)
      if (!length(idx)) next
      for (gene in features) {
        x <- marker_values[idx, gene]
        marker_stats[[length(marker_stats) + 1L]] <- data.frame(
          state = state, gene = gene, n_cells = length(x),
          avg_exp = mean(expm1(x)), pct_exp = 100 * sum(x > 0) / length(x),
          mean_fetched_expression = mean(x),
          assay = DefaultAssay(sc_hd), requested_layer = "data (FetchData default)"
        )
      }
    }
    write.csv(do.call(rbind, marker_stats), "state_marker_expression.csv", row.names = FALSE)
  }
  .hub_write_matrix(sc_hd@meta.data, "annotated_cell_metadata.csv", id = "cell")
  umap_available <- "umap" %in% Reductions(sc_hd)
  write.csv(data.frame(reduction = "umap", available = umap_available,
                       status = if (umap_available) "coordinates_exported" else "absent_in_input_not_recomputed"),
            "UMAP_availability.csv", row.names = FALSE)
  if (umap_available) {
    .hub_write_matrix(Embeddings(sc_hd, "umap"), "annotated_umap_coordinates.csv", id = "cell")
  }
  saveRDS(sc_hd, "ALS_Microglia_State_Annotated.rds")
}

# Differential expression
.hub_de <- function(work_dir, input_rds) {
  sc_de <- .hub_read_seurat(input_rds)
  .hub_require_study_samples(sc_de)
  .hub_object_record(sc_de)
  sc_de$Group <- ifelse(
    sc_de$orig.ident == "133_525_als", "ALS_525",
    ifelse(sc_de$orig.ident == "147_521_als", "ALS_521", "Control"))
  Idents(sc_de) <- "Group"

  contrasts <- list(
    ALS521_vs_Control = c("ALS_521", "Control", "ALS_521 vs Control"),
    ALS525_vs_Control = c("ALS_525", "Control", "ALS_525 vs Control"),
    ALS525_vs_521 = c("ALS_525", "ALS_521", "ALS_525 vs ALS_521")
  )
  .hub_require_groups(sc_de$Group, c("ALS_521", "ALS_525", "Control"), "Group")
  de_results <- lapply(contrasts, function(contrast) {
    FindMarkers(sc_de, ident.1 = contrast[1], ident.2 = contrast[2],
                test.use = de.test_use,
                logfc.threshold = de.logfc_threshold, min.pct = de.min_pct)
  })
  for (key in names(contrasts)) {
    .hub_write_matrix(de_results[[key]], paste0("DE_", key, "_raw.csv"))
  }
  for (key in names(contrasts)) {
    write.csv(.hub_process_de(de_results[[key]], contrasts[[key]][3]),
              paste0("DE_", key, "_final.csv"), row.names = FALSE)
  }
  Idents(sc_de) <- "Group"
  sc_de$Group_AllvsCtrl <- ifelse(sc_de$Group %in% c("ALS_521", "ALS_525"),
                                  "ALS_All", "Control")
  Idents(sc_de) <- "Group_AllvsCtrl"
  de_all <- FindMarkers(sc_de, ident.1 = "ALS_All", ident.2 = "Control",
                        test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1)

  .hub_write_matrix(de_all, "DE_ALSAll_vs_Control_raw.csv")
  # Classify DE using the FindMarkers-adjusted P values.

  classes <- rbind(.hub_classify_de(de_results$ALS525_vs_Control, "ALS_525 vs Control"),
                   .hub_classify_de(de_results$ALS521_vs_Control, "ALS_521 vs Control"),
                   .hub_classify_de(de_all, "ALS_All vs Control"))
  write.csv(classes, "DE_original_significance_classes.csv", row.names = FALSE)
  write.csv(classes[classes$gene == "BSG", , drop = FALSE], "BSG_DE_source_statistics.csv", row.names = FALSE)
  write.csv(data.frame(
    output = c("*_raw.csv and DE_original_significance_classes.csv", "*_final.csv"),
    adjusted_p_source = c("FindMarkers return value, original p_val_adj",
                          "BH recalculated over p_val rows returned by each FindMarkers call"),
    classification_abs_log2FC = c(0.25, de.sig_logfc),
    classification_p_cutoff = c(0.05, de.sig_padj)),
    "DE_pvalue_definitions.csv", row.names = FALSE)
  .hub_write_matrix(sc_de@meta.data, "DE_cell_groups.csv", id = "cell")
}

# BSG expression
.hub_bsg <- function(work_dir, input_rds) {
  sc_hd <- .hub_read_seurat(input_rds)
  .hub_require_states(sc_hd)
  .hub_require_feature(sc_hd, "BSG")
  .hub_object_record(sc_hd)
  bsg_data <- FetchData(sc_hd, vars = c("BSG", "Group", "Trajectory_States"))
  bsg_data$cell <- rownames(bsg_data)
  bsg_data <- bsg_data %>%
    mutate(Study_Group = case_when(
      Group == "Control" ~ "WT iMG",
      Group == "ALS_521" ~ "R521C iMG",
      Group == "ALS_525" ~ "P525L iMG",
      TRUE ~ NA_character_)) %>%
    filter(!is.na(Study_Group)) %>%
    mutate(Study_Group = factor(Study_Group, levels = c("WT iMG", "R521C iMG", "P525L iMG")),
           Trajectory_States = factor(Trajectory_States, levels = state_order))

  my_comparisons <- list(
    c("WT iMG", "R521C iMG"),
    c("WT iMG", "P525L iMG")
  )

  # Use finite, nonnegative BSG values for summaries and tests.
  bsg_data$included_in_source_statistics <- is.finite(bsg_data$BSG) & bsg_data$BSG >= 0
  write.csv(bsg_data, "BSG_expression_by_cell.csv", row.names = FALSE)
  summary_rows <- list()
  test_rows <- list()
  for (state in state_order) {
    sub_data <- bsg_data[bsg_data$Trajectory_States == state, , drop = FALSE]
    for (group in levels(bsg_data$Study_Group)) {
      g <- sub_data[sub_data$Study_Group == group, , drop = FALSE]
      x <- g$BSG[g$included_in_source_statistics]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        state = state, group = group, n_cells = nrow(g), n_tested = length(x),
        n_excluded_by_source_data_boundary = nrow(g) - length(x),
        mean = if (length(x)) mean(x) else NA_real_,
        sd = if (length(x) > 1) sd(x) else NA_real_,
        sem = if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
      )
    }
    for (pair in my_comparisons) {
      d <- sub_data[sub_data$included_in_source_statistics, , drop = FALSE]
      x <- d$BSG[d$Study_Group == pair[1]]
      y <- d$BSG[d$Study_Group == pair[2]]
      test_warning <- character()
      test_error <- NA_character_
      result <- if (length(x) && length(y)) {
        tryCatch(withCallingHandlers(
          stats::wilcox.test(x, y, paired = FALSE),
          warning = function(w) {
            test_warning <<- c(test_warning, conditionMessage(w))
            invokeRestart("muffleWarning")
          }), error = function(e) {
            test_error <<- conditionMessage(e)
            NULL
          })
      } else {
        test_error <- "A comparison group has no eligible cells"
        NULL
      }
      test_rows[[length(test_rows) + 1L]] <- data.frame(
        state = state, group1 = pair[1], group2 = pair[2], n1 = length(x), n2 = length(y),
        statistic = if (is.null(result)) NA_real_ else unname(result$statistic),
        p_value_unadjusted = if (is.null(result)) NA_real_ else result$p.value,
        method = if (is.null(result)) "wilcox.test" else result$method,
        alternative = "two.sided", paired = FALSE, multiple_testing_adjustment = "none",
        exact = "R default; no explicit override in source", warning = paste(unique(test_warning), collapse = "; "),
        status = if (is.null(result)) test_error else if (is.finite(result$p.value)) "completed" else "undefined_p_value"
      )
    }
  }
  write.csv(do.call(rbind, summary_rows), "BSG_mean_SEM_by_state_group.csv", row.names = FALSE)
  write.csv(do.call(rbind, test_rows), "BSG_Wilcoxon_by_state.csv", row.names = FALSE)
}

# Enrichment
.hub_enrichment <- function(work_dir, input_rds) {
  sc_hd <- .hub_read_seurat(input_rds)
  .hub_require_states(sc_hd)
  .hub_require_feature(sc_hd, "BSG")
  .hub_object_record(sc_hd)
  options(timeout = 3600)

  .hub_require_groups(sc_hd$Trajectory_States, met.bsg_subset, "Trajectory_States")
  sc_late_core <- subset(sc_hd, subset = Trajectory_States == met.bsg_subset)

  if (!"Group" %in% colnames(sc_late_core@meta.data)) {
    sc_late_core$Group <- group_from_sample(sc_late_core$orig.ident)
  }
  sc_late_core$Group <- factor(sc_late_core$Group, levels = c("Control", "ALS_521", "ALS_525"))
  sc_late_core$ALS_vs_Ctrl <- ifelse(sc_late_core$Group == "Control", "Control", "ALS")
  Idents(sc_late_core) <- "ALS_vs_Ctrl"

  .hub_require_groups(sc_late_core$ALS_vs_Ctrl, c("ALS", "Control"), "ALS_vs_Ctrl")
  deg_late_core <- FindMarkers(sc_late_core, ident.1 = "ALS", ident.2 = "Control",
                               test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  deg_late_core$gene       <- rownames(deg_late_core)
  deg_late_core$delta.pct  <- deg_late_core$pct.1 - deg_late_core$pct.2
  write.csv(deg_late_core, "DE_ActivatedLateCore_ALS_vs_Control.csv")

  sig_degs <- deg_late_core %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25 & abs(delta.pct) > 0.1)
  up_genes   <- sig_degs %>% filter(avg_log2FC > 0) %>% pull(gene)
  down_genes <- sig_degs %>% filter(avg_log2FC < 0) %>% pull(gene)

  dir.create("Enrichment_Results", showWarnings = FALSE)
  write.csv(sig_degs, "Enrichment_Results/DE_significant_LateCore.csv", row.names = FALSE)
  write.csv(data.frame(gene = up_genes), "Enrichment_Results/up_genes.csv", row.names = FALSE)
  write.csv(data.frame(gene = down_genes), "Enrichment_Results/down_genes.csv", row.names = FALSE)
  gene_ids <- if (length(up_genes)) {
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  } else data.frame(SYMBOL = character(), ENTREZID = character())
  entrez_ids <- gene_ids$ENTREZID
  write.csv(gene_ids, "Enrichment_Results/up_gene_entrez_mapping.csv", row.names = FALSE)
  kegg <- if (length(entrez_ids) > 0) {
    enrichKEGG(gene          = entrez_ids,
               organism      = en.kegg_organism,
               pAdjustMethod = "BH",
               pvalueCutoff  = en.kegg_pv_cutoff,
               qvalueCutoff  = 0.2)
  } else NULL

  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    dir.create("Enrichment_Results", showWarnings = FALSE)
    write.csv(as.data.frame(kegg), "Enrichment_Results/KEGG_up_enrichment.csv", row.names = FALSE)
   }
  if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
    write.csv(data.frame(ID = character(), Description = character(), pvalue = numeric(),
                         p.adjust = numeric(), qvalue = numeric()),
              "Enrichment_Results/KEGG_up_enrichment.csv", row.names = FALSE)
  }
  saveRDS(kegg, "Enrichment_Results/KEGG_up_enrichment.rds")
  write.csv(data.frame(foreground_symbols = length(up_genes),
                       mapped_ids_before_internal_deduplication = length(entrez_ids),
                       reported_terms = if (is.null(kegg)) 0L else nrow(as.data.frame(kegg)),
                       background = "enrichKEGG default; no custom universe in source"),
            "Enrichment_Results/KEGG_status.csv", row.names = FALSE)

  bsg_values <- FetchData(sc_late_core, vars = "BSG")
  .hub_require_cell_alignment(bsg_values, sc_late_core)
  bsg_expr <- bsg_values$BSG
  sc_late_core$BSG_Express <- ifelse(bsg_expr == 0, "No_express", "BSG_express")
  Idents(sc_late_core) <- "BSG_Express"

  .hub_require_groups(sc_late_core$BSG_Express, c("BSG_express", "No_express"), "BSG_Express")
  deg_gsea <- FindMarkers(sc_late_core, ident.1 = "BSG_express", ident.2 = "No_express",
                          test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  deg_gsea$gene <- rownames(deg_gsea)
  deg_gsea <- deg_gsea %>% filter(gene != "BSG")

  deg_gsea$rank_score <- sign(deg_gsea$avg_log2FC) * (-log10(deg_gsea$p_val + 1e-300))
  deg_gsea <- deg_gsea[order(-deg_gsea$rank_score), ]
  gene_list <- setNames(deg_gsea$rank_score, deg_gsea$gene)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)

  dir.create("GSEA_MultiDatabase_Results", showWarnings = FALSE)
  write.csv(data.frame(cell = colnames(sc_late_core), BSG = bsg_expr,
                       BSG_Express = sc_late_core$BSG_Express,
                       Group = sc_late_core$Group, SampleID = sc_late_core$orig.ident),
            "GSEA_MultiDatabase_Results/BSG_binary_cell_groups.csv", row.names = FALSE)
  write.csv(deg_gsea, "GSEA_MultiDatabase_Results/BSG_DE_ranked_without_BSG.csv", row.names = FALSE)
  write.csv(data.frame(gene = names(gene_list), rank_score = unname(gene_list)),
            "GSEA_MultiDatabase_Results/GSEA_ranked_gene_list.csv", row.names = FALSE)
  # Support both msigdbr collection argument names.
  hallmark <- tryCatch(
    msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr(species = "Homo sapiens", category = "H")
  )
  t2g <- hallmark %>% dplyr::select(gs_name, gene_symbol)

  gsea_h <- GSEA(geneList     = gene_list,
                 TERM2GENE    = t2g,
                 pvalueCutoff = en.gsea_pv_cutoff,
                 pAdjustMethod = "BH")

  saveRDS(hallmark, "GSEA_MultiDatabase_Results/MSigDB_Hallmark_snapshot.rds")
  write.csv(t2g, "GSEA_MultiDatabase_Results/MSigDB_Hallmark_TERM2GENE.csv", row.names = FALSE)
  saveRDS(gsea_h, "GSEA_MultiDatabase_Results/GSEA_MSigDB_Hallmark_Result.rds")
  gsea_table <- if (is.null(gsea_h)) {
    data.frame(ID = character(), Description = character(), NES = numeric(), pvalue = numeric(), p.adjust = numeric())
  } else {
    as.data.frame(gsea_h)
  }
  write.csv(gsea_table, "GSEA_MultiDatabase_Results/GSEA_MSigDB_Hallmark_Result.csv", row.names = FALSE)
  focus_pathways <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                      "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                      "HALLMARK_GLYCOLYSIS",
                      "HALLMARK_PI3K_AKT_MTOR_SIGNALING")

  focus_table <- gsea_table[gsea_table$Description %in% focus_pathways, , drop = FALSE]
  write.csv(focus_table, "GSEA_MultiDatabase_Results/GSEA_focus_results.csv", row.names = FALSE)
  write.csv(data.frame(pathway = focus_pathways,
                       returned_by_source_GSEA_cutoff = focus_pathways %in% gsea_table$Description),
            "GSEA_MultiDatabase_Results/GSEA_focus_availability.csv", row.names = FALSE)
  sets <- split(hallmark$gene_symbol, hallmark$gs_name)
  ranked <- gene_list
  N <- length(ranked)
  curves <- list()
  hits <- list()
  for (pathway in focus_pathways[focus_pathways %in% focus_table$Description]) {
    gs <- intersect(sets[[pathway]], names(ranked))
    hit <- names(ranked) %in% gs
    NH <- sum(hit)
    if (NH == 0 || NH == N) next
    P_hit <- cumsum(hit * abs(ranked)) / sum(hit * abs(ranked))
    P_miss <- cumsum(1 - hit) / (N - NH)
    es <- P_hit - P_miss
    curves[[pathway]] <- data.frame(pathway = pathway, Rank = seq_len(N),
                                    P_hit = P_hit, P_miss = P_miss, ES = es)
    hits[[pathway]] <- data.frame(pathway = pathway, Rank = which(hit), gene = names(ranked)[hit])
  }
  curve_table <- if (length(curves)) do.call(rbind, curves) else data.frame(
    pathway = character(), Rank = integer(), P_hit = numeric(), P_miss = numeric(), ES = numeric())
  hit_table <- if (length(hits)) do.call(rbind, hits) else data.frame(
    pathway = character(), Rank = integer(), gene = character())
  write.csv(curve_table, "GSEA_MultiDatabase_Results/GSEA_running_ES.csv", row.names = FALSE)
  write.csv(hit_table, "GSEA_MultiDatabase_Results/GSEA_hit_positions.csv", row.names = FALSE)
}

# Metabolism
.hub_metabolism <- function(work_dir, input_rds) {
  sc_hd <- .hub_read_seurat(input_rds)
  .hub_require_states(sc_hd)
  .hub_require_feature(sc_hd, "BSG")
  .hub_object_record(sc_hd)
  glycolysis_pool <- c(
    "AKR1A1", "ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "GALM", "ADH7",
    "LDHAL6A", "DLAT", "DLD", "ENO1", "ENO2", "ENO3", "ALDH2", "ALDH3A1", "ALDH1B1",
    "ALDH1A3", "FBP1", "ALDH3B1", "ALDH3B2", "ALDH9A1", "ALDH3A2", "ALDOA", "ALDOB",
    "ALDOC", "G6PC1", "GAPDH", "GCK", "GPI", "HK1", "HK2", "HK3", "LDHA", "LDHB",
    "LDHC", "PGAM4", "ALDH7A1", "PCK1", "PCK2", "PDHA1", "PDHA2", "PDHB", "PFKL",
    "PFKM", "PFKP", "PGAM1", "PGAM2", "PGK1", "PGK2", "PGM1", "PKLR", "PKM", "PGM2",
    "ACSS2", "G6PC2", "BPGM", "TPI1", "ACSS1", "FBP2", "LDHAL6B",
    "SLC2A1", "SLC2A3", "SLC16A3", "PDK1", "PFKFB3", "G6PD", "TXN", "SOD1",
    "IER3", "GALE", "MDH1", "BSG", "SRD5A3", "FAM162A", "PYGL", "CYB5A", "SAP30",
    "HK1", "PGAM1", "GAPDH"
  )
  glycolysis_pool <- unique(glycolysis_pool)

  .hub_require_assay(sc_hd, "RNA")
  .hub_require_groups(sc_hd$Trajectory_States, met.bsg_subset, "Trajectory_States")
  sc_met <- subset(sc_hd, subset = Trajectory_States == met.bsg_subset)
  bsg_values <- FetchData(sc_met, vars = "BSG")
  .hub_require_cell_alignment(bsg_values, sc_met)
  bsg_expr <- bsg_values$BSG
  pos_bsg  <- bsg_expr[bsg_expr > 0]
  bsg_median <- median(pos_bsg)

  sc_met$BSG_Group3 <- factor(
    ifelse(bsg_expr == 0, "No_express",
           ifelse(bsg_expr >= bsg_median, "BSG_High", "BSG_Low")),
    levels = c("No_express", "BSG_Low", "BSG_High"))
  Idents(sc_met) <- "BSG_Group3"

  write.csv(data.frame(cell = colnames(sc_met), BSG = bsg_expr, BSG_Group3 = sc_met$BSG_Group3,
                       Group = sc_met$Group, SampleID = sc_met$orig.ident),
            "BSG_three_group_cell_assignments.csv", row.names = FALSE)
  write.csv(data.frame(positive_BSG_median = bsg_median, positive_cells = length(pos_bsg),
                       total_cells = ncol(sc_met), threshold_ties_assigned_to = "BSG_High",
                       fetched_default_assay = DefaultAssay(sc_met)),
            "BSG_three_group_threshold.csv", row.names = FALSE)
  write.csv(as.data.frame(table(BSG_Group3 = sc_met$BSG_Group3)), "BSG_three_group_counts.csv", row.names = FALSE)
  missing_groups <- names(which(table(sc_met$BSG_Group3) == 0))
  if (length(missing_groups)) stop("Source analysis requires all three BSG groups. Missing: ",
                                 paste(missing_groups, collapse = ", "))
  deg_high_vs_none <- FindMarkers(sc_met, ident.1 = "BSG_High", ident.2 = "No_express",
                                  test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  deg_high_vs_none$gene <- rownames(deg_high_vs_none)
  deg_high_vs_none <- deg_high_vs_none[order(-deg_high_vs_none$avg_log2FC), ]

  core_valid <- intersect(met.glycolysis_core, rownames(sc_met))
  candidate_pool <- setdiff(intersect(glycolysis_pool, rownames(sc_met)), core_valid)
  top13 <- head(deg_high_vs_none$gene[deg_high_vs_none$gene %in% candidate_pool],
                25 - length(core_valid))
  final_glycolysis_genes <- c(core_valid, top13)
  cat("glycolysis genes: core", length(core_valid),
      "+ log2FC filled", length(top13), "\n")

  write.csv(deg_high_vs_none, "DE_BSG_High_vs_No_express.csv", row.names = FALSE)
  write.csv(data.frame(gene = glycolysis_pool, available_in_default_assay = glycolysis_pool %in% rownames(sc_met),
                       core_member = glycolysis_pool %in% core_valid,
                       in_remaining_candidate_pool = glycolysis_pool %in% candidate_pool),
            "glycolysis_candidate_pool.csv", row.names = FALSE)
  write.csv(data.frame(gene = final_glycolysis_genes,
                       selection = c(rep("predefined_core", length(core_valid)),
                                     rep("descending_log2FC_fill", length(top13))),
                       order_before_expression_aggregation = seq_along(final_glycolysis_genes)),
            "glycolysis_selected_genes.csv", row.names = FALSE)
  avg_expr <- AverageExpression(sc_met, assays = "RNA", slot = "data",
                                features = final_glycolysis_genes)$RNA

  # Match group labels with underscores or hyphens.
  expected_groups <- c("No_express", "BSG_Low", "BSG_High")
  observed_groups <- colnames(avg_expr)
  column_index <- vapply(expected_groups, function(target) {
    candidates <- unique(c(target, gsub("_", "-", target, fixed = TRUE)))
    idx <- which(observed_groups %in% candidates)
    if (length(idx) != 1L) stop("Cannot unambiguously match RNA mean-expression group: ", target,
                                "; returned columns: ", paste(observed_groups, collapse = ", "))
    idx
  }, integer(1))
  if (anyDuplicated(column_index)) stop("Mean-expression group mapping is not one-to-one.")
  write.csv(data.frame(expected_group = expected_groups, returned_group = observed_groups[column_index]),
            "RNA_mean_expression_column_mapping.csv", row.names = FALSE)
  avg_expr <- avg_expr[, column_index, drop = FALSE]
  colnames(avg_expr) <- expected_groups
  mat_z <- t(scale(t(as.matrix(avg_expr))))
  mat_z[is.nan(mat_z)] <- 0

  .hub_write_matrix(avg_expr, "glycolysis_RNA_group_mean_expression.csv")
  .hub_write_matrix(mat_z, "glycolysis_group_mean_Zscores.csv")
  write.csv(data.frame(assay = "RNA", slot = "data",
                       aggregation = "Seurat AverageExpression; default identity is BSG_Group3",
                       zscore = "Row-wise scale across the three group means; NaN replaced by zero"),
            "glycolysis_expression_definitions.csv", row.names = FALSE)
}

