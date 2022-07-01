#' Cox et al 2014 LFQ data (normalized)
#'
#' A dataframe containing normalized LFQ protein intensity data for 4360
#' proteins in 6 samples
#'
#' @references \url{https://europepmc.org/article/MED/24942700#id609082}
#' @docType data
#' @keywords datasets
#' @name ecoli_norm_df
#' @usage data(ecoli_norm_df)
#' @format A data frame with 4360 rows (proteins) and 6 columns (samples)
NULL

#' Cox et al 2014 LFQ data (fit object)
#'
#' An object of class "MArrayLM" from running find_dep on ecoli_norm_df
#'
#' @references \url{https://europepmc.org/article/MED/24942700#id609082}
#' @docType data
#' @keywords datasets
#' @name ecoli_fit_df
#' @usage data(ecoli_fit_df)
#' @format An object of class "MArrayLM"
NULL

#' Suvarna et al 2021 LFQ data (normalized)
#'
#' A dataframe containing normalized LFQ protein intensity data for 230
#' proteins in 35 samples (a subset of the original data set)
#'
#' @references \url{https://www.frontiersin.org/articles/10.3389/fphys.2021.652799/full#h3}
#' @docType data
#' @keywords datasets
#' @name covid_norm_df
#' @usage data(covid_norm_df)
#' @format A data frame with 230 rows (proteins) and 35 columns (samples)
NULL

#' Suvarna et al 2021 LFQ data (fit object)
#'
#' An object of class "MArrayLM" from running find_dep on covid_norm_df
#'
#' @references \url{https://www.frontiersin.org/articles/10.3389/fphys.2021.652799/full#h3}
#' @docType data
#' @keywords datasets
#' @name covid_fit_df
#' @usage data(covid_fit_df)
#' @format An object of class "MArrayLM"
NULL

