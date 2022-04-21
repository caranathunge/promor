# Filter missing data -----------------------------------------------------
#' Filter out proteins and/ or samples that do not have valid data
#' @description This function removes proteins (rows) that have missing data
#'  (NAs) across all samples and samples (columns) with missing data across all
#'  proteins
#' @param df A \code{raw.df} object (output of \code{\link {create.df}}).
#' @details \itemize{\item This function first removes proteins (rows) with
#' missing values (NAs) across all samples.
#' \item Then if present, it removes samples (columns) with missing data
#' across all proteins }
#' @examples
#' \dontrun{
#' raw <- create.df(file.path = "./proteinGroups.txt")
#' raw_1 <- filter.NA(raw)
#' }
#' @return A \code{raw.df} object or a data frame with no empty (all NAs) rows
#' or columns.
#' @export

filter.NA <- function(df){
    #Remove proteins (rows) with missing values (NA) across all samples
    raw_2 <- df[rowSums(is.na(df)) != ncol(df), ]
    #Remove samples (columns) with missing values (NA) across all proteins
    raw_2<- raw_2[, colSums(is.na(raw_2)) != nrow(raw_2)]
    return(raw_2)
}


# Filter proteins by group level missing data -----------------------------
#' Filter proteins by group level missing data
#' @description This function filters out proteins that exceed a given
#' missing data percentage in each group.
#'
#' @param df A \code{raw.df} object (output of \code{\link{create.df}})
#' @param set.na The percentage of missing data allowed in any group.
#' Default is 0.33.
#'
#' @details This function assumes that column headers in the \code{raw.df}
#'  object provided as \code{df} follow "Group_UniqueSampleID" notation.
#'  \itemize{\item It first
#'  extracts group or condition information from column headers and assigns
#'  samples to different groups.
#'  \item It then removes proteins (rows) from the data frame if the percentage
#'  of NAs in any one of the given groups exceeds the threshold indicated by
#'  \code{set.na} (default is 0.33)}
#'
#' @return A \code{raw.df} object or a data frame of filtered proteins
#'  as rows and sample LFQ intensities as columns.
#'
#' @examples
#'  \dontrun{
#'  ## Create a raw.df object from a proteinGroups.txt file.
#'  raw <- create.df(file.path = "./proteinGroups.txt")
#'
#'  ## Missing data percentage allowed in each group = 0.33
#'  raw_filtered <- filterby.groupNA(raw)
#'
#'  ## Missing data percentage allowed in each group = 0.5
#'  raw_filtered1 <- filterby.groupNA(raw, set.na = 0.5)
#'  }
#' @export

filterby.groupNA <- function(df,
                             set.na = 0.33){

  #Extract group information from sample names in the dataframe x
  group <- factor(c(sapply(strsplit(colnames(df), "_"), getElement, 1)))

  #Transpose the data frame. Columns are now proteins and rows are samples.
  df_trans <- as.data.frame(t(df))

  #Add a new column with the group information.
  #Group column is the rightmost column in the data frame.
  df_trans$Group <- group

  #Get the number of proteins
  n_proteins <- ncol(df_trans)-1

  #First calculate the mean  number of missing values for each group.
  #This outputs a data frame of mean NA: groups as rows and proteins as columns.
  df_na <- aggregate(is.na(df_trans[, 1:n_proteins]),
                     list(df_trans$Group), mean)

  #Remove the first column that contains the group name in the data frame
  df_na[,1]<- NULL

  #Make a list of proteins with >33% NA
  rem_prot <- as.list(colnames(Filter(function(y) any(as.numeric(y) > set.na),
                                    df_na)))

  #Subset the dataframe by removing proteins in the list with >33% NA
  raw_3 <- df[!rownames(df) %in% rem_prot, ]
  return(raw_3)
}


# Output proteins only present in one group -------------------------------

#' Proteins that are only expressed in a given group
#' @description This function outputs a list of proteins that are only
#' expressed (present) in one user-specified group while not expressed
#' (completely absent) in another user-specified group.
#'
#' @param df A \code{raw.df} object (output of \code{\link{create.df}})
#' @param abs.group Name of the group in which proteins are not expressed.
#' @param pres.group Name of the group in which proteins are expressed.
#' @param set.na The percentage of missing data allowed in \code{pres.group}.
#' Default is 0.33.
#' @param save Logical. If \code{TRUE} (default), it saves the output in a text
#' file named "Group_\code{pres.group}_only.txt."
#'
#' @details Note: \code{onegroup.only} function assumes that column headers in
#' the \code{raw.df} object provided as \code{df} follow "Group_UniqueSampleID"
#' notation. \itemize{\item Given a pair of groups, \code{onegroup.only}
#' function finds proteins that are only expressed in \code{pres.group} while
#' completely absent or not expressed in \code{abs.group}.}
#'
#' @return A list of majority protein IDs.
#'
#' @examples
#' \dontrun{
#' ## Create a raw.df object from a proteinGroups.txt file.
#' raw <- create.df(file.path = "./proteinGroups.txt")
#'
#' ## Save a list of proteins only expressed in group B, but absent in group A.
#' onegroup.only(raw, abs.group = "A", pres.group = "B")
#'
#' ## Save the above list in a variable.
#' protein_list <- onegroup.only(raw, abs.group = "A", pres.group = "B")
#' }

#' @export
onegroup.only<- function(df,
                         abs.group,
                         pres.group,
                         set.na= 0.33,
                         save = TRUE){

  #Extract group information from sample names in the data frame x
  group <- factor(c(sapply(strsplit(colnames(df), "_"), getElement,1)))

  #Transpose the data frame: columns are now proteins and rows are samples.
  df_trans <- as.data.frame(t(df))

  #Add a new column with the group information.
  #Group column is the rightmost column in the data frame.
  df_trans$Group <- group

  #Number of proteins
  n_proteins <- ncol(df_trans)-1

  #Check if there are proteins with 100% missing data in one group and
  #>66% valid data in the other group
  group_only <- as.data.frame(
    ifelse(
      colSums(is.na(df_trans[df_trans$Group == abs.group, 1:n_proteins]))
      /nrow(df_trans[df_trans$Group == abs.group,]) == 1 &
      colSums(is.na(df_trans[df_trans$Group == pres.group, 1:n_proteins]))
      /nrow(df_trans[df_trans$Group == pres.group,]) <= set.na,
      print("TRUE"), print("FALSE")))

  group_only$proteins <- rownames(group_only)

  # Print out a list of proteins only present in user specified group
  prot.list <- group_only[group_only[1] == "TRUE",2]
  return(prot.list)

  if(save == TRUE){
    cat(prot.list,
        file = paste0("Group_",pres.group,"_only.txt"), sep="\n")
  }

}
