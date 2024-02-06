library(sceptre)
library(sceptreIGVF)
library(SeuratData)

# Analysis options ------------------------------------------------------------
NUM_GENES_TO_KEEP <- 289
NUM_MT_GENES_TO_KEEP <- 1
NUM_CELLS_TO_KEEP <- 10000
NUM_DISCOVERY_PAIRS <- 100
PC_GENES <- c("STAT1", "JAK2", "CMTM6", "STAT2", 
              "UBE2L6", "STAT3", "TNFRSF14", 
              "IFNGR2", "NFKBIA")

# Load and preprocess data -----------------------------------------------------
eccite <- LoadData(ds = "thp1.eccite")

set.seed(1)

# Extract gene and gRNA expressions
gene_expressions <- eccite@assays$RNA@counts
grna_expressions <- as.matrix(eccite@assays$GDO@counts)-1 # it appears that 1 was added to all counts

# Extra gRNA target data frame
grna_target_data_frame_lowmoi <- eccite@meta.data |> 
  dplyr::select(guide_ID, gene) |> 
  unique() |> 
  dplyr::arrange(gene) |>
  dplyr::rename(grna_id = guide_ID, grna_target = gene) |>
  dplyr::mutate(grna_target = as.character(grna_target),
                grna_target = ifelse(grna_target == "NT", "non-targeting", grna_target)) |>
  `rownames<-`(NULL)

# Extra gene and gRNA IDs
gene_ids <- rownames(gene_expressions)
grna_ids <- unique(grna_target_data_frame_lowmoi$grna_id) |> as.character()

# Construct PC pairs
pc_pairs_lowmoi <- tibble::tibble(
  response_id = PC_GENES,
  grna_group = response_id) |> 
  as.data.frame()

# Sample genes to keep
mt_genes <- grep(pattern = "^MT-", x = gene_ids, value = TRUE)
my_genes <- unique(c(sample(gene_ids, NUM_GENES_TO_KEEP), 
                     sample(mt_genes, NUM_MT_GENES_TO_KEEP), 
                     pc_pairs_lowmoi$response_id))

# Sample cells to keep
cell_ids <- colnames(gene_expressions)
cells_to_keep <- sample(cell_ids, size = NUM_CELLS_TO_KEEP)

# Subset data
response_matrix_lowmoi <- gene_expressions[my_genes,cells_to_keep]
grna_matrix_lowmoi <- grna_expressions[grna_ids,cells_to_keep]
covariate_data_frame_lowmoi <- eccite@meta.data[cells_to_keep,] |>
  dplyr::select(replicate, orig.ident) |>
  dplyr::rename(lane = orig.ident)

# Import and process via sceptre -----------------------------------------------

# Import into sceptre
sceptre_object_lowmoi <- import_data(
  response_matrix = response_matrix_lowmoi,
  grna_matrix = grna_matrix_lowmoi,
  extra_covariates = covariate_data_frame_lowmoi,
  grna_target_data_frame = grna_target_data_frame_lowmoi,
  moi = "low"
)

# Define positive control pairs and discovery pairs
positive_control_pairs_lowmoi <- construct_positive_control_pairs(
  sceptre_object = sceptre_object_lowmoi
)
discovery_pairs_lowmoi <- construct_trans_pairs(
  sceptre_object = sceptre_object_lowmoi,
  positive_control_pairs = positive_control_pairs_lowmoi,
  pairs_to_exclude = "pc_pairs"
)

# Run sceptre analysis
sceptre_object_lowmoi <- sceptre_object_lowmoi |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs_lowmoi,
    positive_control_pairs = positive_control_pairs_lowmoi,
    formula = formula(~ log(response_n_nonzero) + log(response_n_umis) + replicate)
  ) |>
  assign_grnas() |>
  run_qc() |>
  run_power_check(parallel = TRUE) |>
  run_discovery_analysis(parallel = TRUE)

# Convert to MuData -----------------------------------------------------------

# Convert to MuData inputs and outputs
mudata_list_lowmoi <- sceptre_object_to_mudata_inputs_outputs(
  sceptre_object = sceptre_object_lowmoi,
  num_discovery_pairs = NUM_DISCOVERY_PAIRS,
  guide_capture_method = "direct capture"
)

# Save to disk
save_mudata_list(mudata_list = mudata_list_lowmoi, path = "data", prefix = "papalexi_")