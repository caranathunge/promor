## code to prepare `ecoli_fit_df` dataset goes here
raw <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
)
raw_filtered <- filterbygroup_na(raw)
imp_df_mp <- impute_na(raw_filtered)
ecoli_norm_df <- normalize_data(imp_df_mp)
ecoli_fit_df <- find_dep(ecoli_norm_df)
usethis::use_data(ecoli_fit_df, overwrite = TRUE)
