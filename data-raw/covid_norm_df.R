## code to prepare `covid_norm_df` dataset goes here
covid_raw <- create_df(
prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg3.txt",
exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed3.txt",
)
covid_filt <- filterbygroup_na(covid_raw)
covid_imp_df <- impute_na(covid_filt, method = "kNN")
covid_norm_df <- normalize_data(covid_imp_df)
usethis::use_data(covid_norm_df, overwrite = TRUE)
