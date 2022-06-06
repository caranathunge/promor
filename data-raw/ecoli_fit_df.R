## code to prepare `ecoli_fit_df` dataset goes here
raw <- create_df(prot_groups = system.file("extdata", "ecoli_proteinGroups.txt",
                                           package = "promor"),
                 exp_design = system.file("extdata", "expDesign.txt",
                                          package = "promor"))
raw_filtered <- filterbygroup_na(raw)
imp_df_mp <- impute_na(raw_filtered)
ecoli_norm_df <- normalize_data(imp_df_mp)
ecoli_fit_df <- find_dep(ecoli_norm_df)
usethis::use_data(ecoli_fit_df, overwrite = TRUE)
