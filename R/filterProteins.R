# Filter proteins by group level missing data -----------------------------
#' Filter proteins by group level missing data
#' @description This function filters out proteins based on missing data
#' at the group level.
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom stats aggregate
#'
#' @param raw_df A \code{raw_df} object (output of \code{\link{create_df}})
#' @param set_na The proportion of missing data allowed.
#' Default is 0.34 (one third of the samples in the group).
#' @param filter_condition If set to \code{"each"}, proteins that exceed
#' the missing value proportion threshold set by \code{set_na} in each group
#' will be removed (lenient).
#' If set to \code{"either"}(default), proteins that exceed the missing value
#' proportion threshold set by \code{set_na} in at least one group will be
#' removed (stringent).
#'
#' @details
#'  \itemize{\item This function first
#'  extracts group or condition information from the \code{raw_df} object and
#'  assigns samples to their groups.
#'  \item If \code{filter_condition = "each"}, it then removes proteins (rows)
#'  from the data frame if the proportion of NAs in **each** group exceeds the
#'  threshold indicated by \code{set_na} (default is 0.34). This option is
#'  more lenient in comparison to \code{filter_condition = "either"}, where
#'  proteins that exceeds the missing data threshold in **either** group gets
#'  removed from the data frame.}
#'
#' @return A \code{raw_df} object.
#' @seealso \code{\link{create_df}}
#'
#' @examples
#' \dontrun{
#'
#' # Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/PXD000279_proteinGroups.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/PXD000279_expDesign.txt",
#' )
#'
#' ## Remove proteins that exceed 34% NAs in either group (default)
#' rawdf_filt1 <- filterbygroup_na(raw_df)
#'
#' ## Remove proteins that exceed 34% NAs in each group
#' rawdf_filt2 <- filterbygroup_na(raw_df, filter_condition = "each")
#'
#' ## Proportion of samples with NAs allowed in each group = 0.5
#' rawdf_filt3 <- filterbygroup_na(raw_df, set_na = 0.5, filter_condition = "each")
#' }
#'
#' @export

filterbygroup_na <- function(raw_df,
                             set_na = 0.34,
                             filter_condition = "either") {

  # Extract number of row
  orig_rows <- nrow(raw_df)

  # Extract group information from sample names in the dataframe x
  group <- factor(c(sapply(strsplit(colnames(raw_df), "_"), getElement, 1)))

  # Transpose the data frame. Columns are now proteins and rows are samples.
  transdf <- as.data.frame(t(raw_df))

  # Add a new column with the group information.
  # Group column is the rightmost column in the data frame.
  transdf$Group <- group

  # Get the number of proteins
  n_proteins <- ncol(transdf) - 1

  # First calculate the mean  number of missing values for each group.
  # This outputs a data frame of mean NA: groups as rows and proteins as
  # columns.
  df_na <- aggregate(
    is.na(transdf[, 1:n_proteins]),
    list(transdf$Group), mean
  )

  # Remove the first column that contains the group name in the data frame
  df_na[, 1] <- NULL

  #Remove proteins that don't meet the set_na in at least one group
  if (filter_condition == "either"){
    # Make a list of proteins where at least one group NA > set_na
    rem_prot <- as.list(colnames(Filter(
      function(y) any(as.numeric(y) > set_na),
      df_na
  )))

  # Subset the dataframe by removing proteins in the list with >33% NA
  raw_3 <- raw_df[!rownames(raw_df) %in% rem_prot, ]
  message(paste0(
    orig_rows - nrow(raw_3),
    " proteins with higher than ",
    set_na * 100,
    "% NAs in at least one group removed."))
  }
  if (filter_condition == "each"){
    # Make a list of proteins where at least one group NA < set_na
    keep_prot <- as.list(colnames(Filter(
      function(y) any(as.numeric(y) < set_na),
      df_na
    )))

    # Subset the dataframe by keeping proteins in the list
    raw_3 <- raw_df[rownames(raw_df) %in% keep_prot, ]
    message(paste0(
      orig_rows - nrow(raw_3),
      " proteins with higher than ",
      set_na * 100,
      "% NAs in each group removed."))
  }
  return(raw_3)
}


# Output proteins only present in one group -------------------------------

#' Proteins that are only expressed in a given group
#' @description This function outputs a list of proteins that are only
#' expressed (present) in one user-specified group while not expressed
#' (completely absent) in another user-specified group.
#'
#' @author Chathurani Ranathunge
#'
#' @param raw_df A \code{raw_df} object (output of \code{\link{create_df}})
#' @param abs_group Name of the group in which proteins are not expressed.
#' @param pres_group Name of the group in which proteins are expressed.
#' @param set_na The percentage of missing data allowed in \code{pres_group}.
#' Default is 0.34 (one thrid of the samples in the group).
#' @param save Logical. If \code{TRUE} (default), it saves the output in a text
#' file named "Group_\code{pres_group}_only.txt."
#'
#' @details Note: \code{onegroup_only} function assumes that column names in
#' the \code{raw_df} object provided as \code{df} follow "Group_UniqueSampleID"
#' notation. (Use \code{head(raw_df)} to check the structure of your
#' \code{raw_df} object.)
#' \itemize{\item Given a pair of groups, \code{onegroup_only}
#' function finds proteins that are only expressed in \code{pres_group} while
#' completely absent or not expressed in \code{abs_group}.}
#'
#' @return A list of majority protein IDs.
#'
#' @examples
#' \dontrun{
#'
#' # Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/PXD000279_proteinGroups.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/PXD000279_expDesign.txt",
#' )
#'
#' ## Save a list of proteins only expressed in group L, but absent in group H.
#' onegroup_only(raw_df,abs_group = "H", pres_group = "L")
#'
#' }
#'
#' @export
onegroup_only <- function(raw_df,
                          abs_group,
                          pres_group,
                          set_na = 0.34,
                          save = TRUE) {

  # Extract group information from sample names in the data frame x
  group <- factor(c(sapply(strsplit(colnames(raw_df), "_"), getElement, 1)))

  # Transpose the data frame: columns are now proteins and rows are samples.
  transdf <- as.data.frame(t(raw_df))

  # Add a new column with the group information.
  # Group column is the rightmost column in the data frame.
  transdf$Group <- group

  # Number of proteins
  n_proteins <- ncol(transdf) - 1

  # Check if there are proteins with 100% missing data in one group and
  #>1-set_na% valid data in the other group
  group_only <- as.data.frame(
    ifelse(
      colSums(is.na(transdf[transdf$Group == abs_group, 1:n_proteins]))
      / nrow(transdf[transdf$Group == abs_group, ]) == 1 &
        colSums(is.na(transdf[transdf$Group == pres_group, 1:n_proteins]))
        / nrow(transdf[transdf$Group == pres_group, ]) <= set_na,
      print("TRUE"), print("FALSE")
    )
  )

  group_only$proteins <- rownames(group_only)

  # Print out a list of proteins only present in user specified group
  prot_list <- group_only[group_only[1] == "TRUE", 2]

  if (identical(prot_list, character(0))) {
    message(paste0(
      "None of the proteins are expressed only in ",
      pres_group
    ))
  } else {
    if (save == TRUE) {
      cat(prot_list,
        file = paste0("Group_", pres_group, "_only.txt"), sep = "\n"
      )
    }
    return(prot_list)
  }
}
