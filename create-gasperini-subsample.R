library(sceptre)
library(sceptreIGVF)
library(MuData)
library(SummarizedExperiment)

# Analysis options ------------------------------------------------------------
ENHANCERS_TO_KEEP <- c("chr8.847_top_two", "chr9.1594_top_two", 
                       "chr9.2869_top_two", "chr9.3633_top_two", 
                       "chr9.871_top_two")
NUM_PC_PAIRS <- 10
NUM_NT_GRNAS <- 25
NUM_GENES_TO_KEEP <- 289
NUM_MT_GENES_TO_KEEP <- 2
NUM_CELLS_TO_KEEP <- 10000
CIS_DISTANCE_THRESHOLD <- 5e6
NUM_DISCOVERY_PAIRS <- 100

# Load and preprocess data -----------------------------------------------------
gene_table <- gene_position_data_frame_grch38
gasperini_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/processed/")
gasperini_dir_raw <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/raw/")
gasperini_dir_intermediate <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/intermediate/")

# 0. set the Gasperini directory load the group, gene, and grna data
gasp_fp <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/processed/")
pairs_grouped <- readRDS(paste0(gasp_fp, "/pairs_grouped.rds")) |> dplyr::distinct()
covariate_info <- readRDS(paste0(gasperini_dir_intermediate, "/cell_covariates.rds")) |>
  dplyr::rename(cell_id = cell) |>
  dplyr::select(cell_id, prep_batch, within_batch_chip, within_chip_lane)
all_deg_results <- readr::read_tsv(paste0(gasperini_dir_raw, "/GSE120861_all_deg_results.at_scale.txt"))
gene_info <- all_deg_results |> 
  dplyr::select(ENSG, gene_short_name, target_gene.chr, 
                target_gene.start, target_gene.stop) |> 
  dplyr::group_by(ENSG) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

gene_odm_fp <- paste0(gasp_fp, "gene/matrix.odm")
gene_metadata_fp <- paste0(gasp_fp, "gene/metadata.rds")
gene_odm <- ondisc::read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)

grna_odm_fp <- paste0(gasp_fp, "grna_expression/matrix.odm")
grna_metadata_fp <-  paste0(gasp_fp, "grna_expression/metadata.rds")
grna_odm <- ondisc::read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)
grna_feature_df <- grna_odm |>
  ondisc::get_feature_covariates()
grna_feature_df <- grna_feature_df |>
  dplyr::mutate(grna_id = rownames(grna_feature_df)) |>
  `rownames<-`(NULL)

# 1. Construct the discovery pairs data frame, which contains both cis and positive control pairs
set.seed(6) 
cis_pairs <- pairs_grouped |>
  dplyr::filter(grna_group %in% ENHANCERS_TO_KEEP, gene_id %in% gene_table$response_id)
my_pc_grna_groups <- pairs_grouped |>
  dplyr::filter(site_type == "pos_cntrl") |>
  dplyr::pull(grna_group) |>
  sample(NUM_PC_PAIRS)
pc_pairs <- pairs_grouped |>
  dplyr::filter(grna_group %in% my_pc_grna_groups, site_type == "pos_cntrl",
                gene_id %in% gene_table$response_id)
targeting_pairs <- rbind(cis_pairs, pc_pairs) |>
  dplyr::rename(type = site_type, response_id = gene_id)

# 2. sample some set of negative control pairs
nt_grnas <- grna_feature_df |>
  dplyr::filter(target_type == "non-targeting") |>
  dplyr::sample_n(NUM_NT_GRNAS) |> 
  dplyr::pull(grna_id)

# 3. construct the grna group table
grna_group_data_frame_highmoi <- rbind(grna_feature_df |>
                                         dplyr::filter(grna_group %in% targeting_pairs$grna_group) |>
                                         dplyr::select(grna_id, grna_group),
                                       data.frame(grna_id = nt_grnas, grna_group = "non-targeting")) |>
  dplyr::arrange(grna_group)

# 4. add chromosomal locations to the grna group table
grna_loc_info <- readr::read_tsv(
  file = paste0(gasperini_dir_raw, "/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"),
  col_names = c("chr", "start", "end", "grna_group", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"),
  skip = 1) |>
  dplyr::select(chr, start, end, grna_group) |>
  dplyr::distinct()
targeting_grna_grps <- grna_group_data_frame_highmoi |> dplyr::filter(grna_group != "non-targeting") |> dplyr::pull()
target_group_loc_info <- grna_loc_info |> dplyr::filter(grna_group %in% targeting_grna_grps) |>
  dplyr::mutate(start = as.integer(start), end = as.integer(end))
grna_group_data_frame_highmoi <- dplyr::left_join(x = grna_group_data_frame_highmoi,
                                                  y = target_group_loc_info, by = c("grna_group"))

# 5. rename some of the grna groups
renamed_grna_group <- pc_pairs$gene_id[match(x = grna_group_data_frame_highmoi$grna_group, table = pc_pairs$grna_group)]
grna_group_data_frame_highmoi$grna_group[!is.na(renamed_grna_group)] <- renamed_grna_group[!is.na(renamed_grna_group)]
enh_group_idxs <- grep(pattern = "chr", x = grna_group_data_frame_highmoi$grna_group)
enh_groups <- grna_group_data_frame_highmoi$grna_group[enh_group_idxs]
grna_group_data_frame_highmoi$grna_group[enh_group_idxs] <- factor(
  x = enh_groups, 
  levels = unique(enh_groups),
  labels = paste0("candidate_enh_", seq_along(unique(enh_groups)))
  ) |>
  as.character()

# 6. reset the start and end positions
grna_group_data_frame_highmoi <- dplyr::group_by(grna_group_data_frame_highmoi, grna_group) |>
  dplyr::group_modify(.f = function(tbl, key) {
    if (!all(is.na(tbl$start))) {
      group_start <- min(tbl$start)
      group_end <- max(tbl$end)
      if (group_end - group_start == 1L) group_end <- group_end + 30
      posits <- as.integer(floor(seq(group_start, group_end, length.out = nrow(tbl) + 1)))
      tbl$start <- posits[seq(1, nrow(tbl))]
      tbl$end <- posits[seq(2, nrow(tbl) + 1)]
    }
    return(tbl)
  }) |> dplyr::relocate(grna_id)

# 7. determine cells that contain at least one gRNA
my_grna_ids <- unique(grna_group_data_frame_highmoi$grna_id)
grna_matrix <- grna_odm[[my_grna_ids,]]
cell_ids <- which(apply(X = as.matrix(grna_matrix >= 5), 2, any)) |> 
  sample(size = NUM_CELLS_TO_KEEP)

# 8. construct the response and grna matrices; downsample cells
multimodal_odm <- ondisc::multimodal_ondisc_matrix(
  covariate_ondisc_matrix_list = list(grna = grna_odm, gene = gene_odm)
)
multimodal_odm_downsample <- multimodal_odm[,cell_ids]

# 9. get the gene matrix, adding a couple MT genes for good measure
gene_sub <- ondisc::get_modality(multimodal_odm_downsample, "gene")
all_gene_ids <- ondisc::get_feature_ids(gene_sub)
mt_genes <- gene_table |> 
  dplyr::filter(chr == "chrM", response_id %in% all_gene_ids) |> 
  dplyr::pull(response_id)
my_gene_ids <- c(unique(targeting_pairs$response_id), 
                 sample(mt_genes, NUM_MT_GENES_TO_KEEP))

gene_matrix <- gene_sub[[my_gene_ids,]]
rownames(gene_matrix) <- my_gene_ids
colnames(gene_matrix) <- gene_sub |> ondisc::get_cell_barcodes()

# 10. get the grna matrix
grna_sub <- ondisc::get_modality(multimodal_odm_downsample, "grna")
my_grna_ids <- unique(grna_group_data_frame_highmoi$grna_id)
grna_matrix <- grna_sub[[my_grna_ids,]]
rownames(grna_matrix) <- my_grna_ids
colnames(grna_matrix) <- grna_sub |> ondisc::get_cell_barcodes()

# 11. Compute the covariate matrix
covariate_matrix <- multimodal_odm_downsample |>
  ondisc::get_cell_covariates() |>
  dplyr::select(prep_batch = gene_batch) |>
  tibble::rownames_to_column(var = "cell_id") |>
  dplyr::left_join(covariate_info, by = c("cell_id", "prep_batch")) |>
  tibble::column_to_rownames(var = "cell_id")

# 12. sort according to batch; also remove cells containing no gRNAs
cell_order <- order(covariate_matrix$prep_batch)
response_matrix_highmoi <- gene_matrix[,cell_order]
grna_matrix_highmoi <- grna_matrix[,cell_order]
covariate_data_frame_highmoi <- covariate_matrix[cell_order,,drop = FALSE]

# 13. rename the data objects
gene_info_subset <- dplyr::tibble(ENSG = rownames(response_matrix_highmoi)) |>
  dplyr::left_join(dplyr::tibble(ENSG = rownames(response_matrix_highmoi)) |>
  dplyr::left_join(gene_info, by = "ENSG"), by = "ENSG") |>
  tibble::column_to_rownames(var = "ENSG") |>
  dplyr::rename(symbol = gene_short_name, 
                gene_chr = target_gene.chr, 
                gene_start = target_gene.start,
                gene_end = target_gene.stop)
grna_target_data_frame_highmoi <- grna_group_data_frame_highmoi |> 
  dplyr::rename("grna_target" = "grna_group") |> 
  as.data.frame()

# Import and process via sceptre -----------------------------------------------

# Import data
sceptre_object_highmoi <- import_data(
  response_matrix = response_matrix_highmoi,
  grna_matrix = grna_matrix_highmoi,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = covariate_data_frame_highmoi,
)

# Construct positive control pairs and discovery pairs
positive_control_pairs <- construct_positive_control_pairs(sceptre_object_highmoi)
discovery_pairs <- construct_cis_pairs(
  sceptre_object_highmoi,
  positive_control_pairs = positive_control_pairs,
  distance_threshold = CIS_DISTANCE_THRESHOLD
)

# Analyze via sceptre
sceptre_object_highmoi <- sceptre_object_highmoi |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs,
    side = "left"
  ) |>
  assign_grnas(parallel = TRUE) |>
  run_qc() |>
  run_power_check(parallel = TRUE) |>
  run_discovery_analysis(parallel = TRUE)
  
# Convert to MuData -----------------------------------------------------------

# Convert to MuData inputs and outputs
mudata_list_highmoi <- sceptre_object_to_mudata_inputs_outputs(
  sceptre_object = sceptre_object_highmoi,
  num_discovery_pairs = NUM_DISCOVERY_PAIRS, 
  gene_info = gene_info_subset,
  guide_capture_method = "CROP-seq"
)

# Save to disk
save_mudata_list(mudata_list = mudata_list_highmoi, path = "data", prefix = "gasperini_")
