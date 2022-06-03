# Filter proteins by group level missing data -----------------------------
#' Filter proteins by group level missing data
#' @description This function filters out proteins that exceed a given
#' missing data percentage in each group.
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom stats aggregate
#'
#' @param df A \code{raw_df} object (output of \code{\link{create_df}})
#' @param set_na The percentage of missing data allowed in any group.
#' Default is 0.33.
#'
#' @details This function assumes that column headers in the \code{raw_df}
#'  object provided as \code{df} follow "Group_UniqueSampleID" notation.
#'  \itemize{\item It first
#'  extracts group or condition information from column headers and assigns
#'  samples to different groups.
#'  \item It then removes proteins (rows) from the data frame if the percentage
#'  of NAs in any one of the given groups exceeds the threshold indicated by
#'  \code{set_na} (default is 0.33)}
#'
#' @return A \code{raw_df} object or a data frame of filtered proteins
#'  as rows and sample LFQ intensities as columns.
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt")
#'
#' ## Missing data percentage allowed in each group = 0.33
#' raw_filtered <- filterbygroup_na(raw)
#'
#' ## Missing data percentage allowed in each group = 0.5
#' raw_filtered1 <- filterbygroup_na(raw, set_na = 0.5)
#' }
#' @export

filterbygroup_na <- function(df,
                             set_na = 0.33) {

  # Extract group information from sample names in the dataframe x
  group <- factor(c(sapply(strsplit(colnames(df), "_"), getElement, 1)))

  # Transpose the data frame. Columns are now proteins and rows are samples.
  transdf <- as.data.frame(t(df))

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

  # Make a list of proteins with >33% NA or > set_na
  rem_prot <- as.list(colnames(Filter(
    function(y) any(as.numeric(y) > set_na),
    df_na
  )))

  # Subset the dataframe by removing proteins in the list with >33% NA
  raw_3 <- df[!rownames(df) %in% rem_prot, ]
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
#' @param df A \code{raw_df} object (output of \code{\link{create_df}})
#' @param abs_group Name of the group in which proteins are not expressed.
#' @param pres_group Name of the group in which proteins are expressed.
#' @param set_na The percentage of missing data allowed in \code{pres_group}.
#' Default is 0.33.
#' @param save Logical. If \code{TRUE} (default), it saves the output in a text
#' file named "Group_\code{pres_group}_only.txt."
#'
#' @details Note: \code{onegroup_only} function assumes that column headers in
#' the \code{raw_df} object provided as \code{df} follow "Group_UniqueSampleID"
#' notation. \itemize{\item Given a pair of groups, \code{onegroup_only}
#' function finds proteins that are only expressed in \code{pres_group} while
#' completely absent or not expressed in \code{abs_group}.}
#'
#' @return A list of majority protein IDs.
#'
#' @examples
#' \dontrun{
#' ## Create a raw_df object from a proteinGroups.txt file.
#' raw <- create_df(file.path = "./proteinGroups.txt")
#'
#' ## Save a list of proteins only expressed in group B, but absent in group A.
#' onegroup_only(raw, abs_group = "A", pres_group = "B")
#'
#' ## Save the above list in a variable.
#' protein_list <- onegroup_only(raw, abs_group = "A", pres_group = "B")
#' }
#'
#' @export
onegroup_only <- function(df,
                          abs_group,
                          pres_group,
                          set_na = 0.33,
                          save = TRUE) {

  # Extract group information from sample names in the data frame x
  group <- factor(c(sapply(strsplit(colnames(df), "_"), getElement, 1)))

  # Transpose the data frame: columns are now proteins and rows are samples.
  transdf <- as.data.frame(t(df))

  # Add a new column with the group information.
  # Group column is the rightmost column in the data frame.
  transdf$Group <- group

  # Number of proteins
  n_proteins <- ncol(transdf) - 1

  # Check if there are proteins with 100% missing data in one group and
  #>66% valid data in the other group
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
    message(paste0("None of the proteins are expressed only in ",
                   pres_group))
  } else {
    if (save == TRUE) {
      cat(prot_list,
        file = paste0("Group_", pres_group, "_only.txt"), sep = "\n"
      )
    }
    return(prot_list)
  }
}
