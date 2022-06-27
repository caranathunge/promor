## code to prepare `ecoli_norm_df` dataset goes here
raw <- create_df(
prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/PXD000279_proteinGroups.txt",
exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/PXD000279_expDesign.txt",
)
raw_filtered <- filterbygroup_na(raw)
imp_df_mp <- impute_na(raw_filtered)
ecoli_norm_df <- normalize_data(imp_df_mp)
usethis::use_data(ecoli_norm_df, overwrite = TRUE)
