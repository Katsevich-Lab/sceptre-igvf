library(sceptre)
library(sceptreIGVF)
library(MuData)

############### low-MOI ###############

cat("Creating MuData for Papalexi subset...\n")

# 1. import data
sceptre_object_lowmoi <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  extra_covariates = lowmoi_example_data$extra_covariates,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  moi = "low"
)
positive_control_pairs <- construct_positive_control_pairs(sceptre_object_lowmoi)
discovery_pairs <- construct_trans_pairs(
  sceptre_object = sceptre_object_lowmoi,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude = "pc_pairs"
)

# 2-4. set analysis parameters, assign gRNAs, run qc
sceptre_object_lowmoi <- sceptre_object_lowmoi |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs
  ) |>
  assign_grnas() |>
  run_qc(p_mito_threshold = 0.075)

mae_lowmoi <- sceptre_object_to_mudata(sceptre_object_lowmoi)
writeH5MU(mae_lowmoi, "data/papalexi_subset.h5mu")

############### high-MOI ###############

cat("Creating MuData for Gasperini subset...\n")

# 1. import data
sceptre_object_highmoi <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = highmoi_example_data$grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)
positive_control_pairs <- construct_positive_control_pairs(sceptre_object_highmoi)
discovery_pairs <- construct_cis_pairs(sceptre_object_highmoi,
                                       positive_control_pairs = positive_control_pairs,
                                       distance_threshold = 5e6
)

# 2-4. set analysis parameters, assign gRNAs, run qc
sceptre_object_highmoi <- sceptre_object_highmoi |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs,
    side = "left"
  ) |>
  assign_grnas(method = "thresholding", threshold = 1, parallel = TRUE) |>
  run_qc(p_mito_threshold = 0.075)

mae_highmoi <- sceptre_object_to_mudata(sceptre_object_highmoi)

writeH5MU(mae_highmoi, "data/gasperini_subset.h5mu")

cat("Done.\n")