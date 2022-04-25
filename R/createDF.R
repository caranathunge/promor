# Create data frame with LFQ intensities ----------------------------------
#' Create a data frame of protein intensities
#' @description This function creates a data frame of label-free quantitative
#' (LFQ) protein intensities from MaxQuant's proteinGroups.txt file.
#'
#' @importFrom utils read.csv
#'
#' @param file.path File path to proteinGroups.txt file produced by MaxQuant.
#' @param filter.prot Logical. If \code{TRUE} (default), filters out
#' reverse proteins, proteins only identified by site, potential contaminants,
#' and proteins identified with less than the minimum number of unique peptides
#' indicated by \code{uniq.pep}.
#' @param uniq.pep Numerical. The minimum number of unique peptides required to
#' identify a protein (default is 2). Proteins that are identified by less than
#' this number of unique peptides are filtered out.
#' @param zero.NA Logical. If \code{TRUE} (default), zeros are considered missing
#' values and replaced with NAs.
#' @param log.tr Logical. If \code{TRUE} (default), intensity values are log
#' transformed to the base indicated by \code{base}.
#' @param base Numerical. Logarithm base. Default is 2.
#'
#' @details
#' \itemize{\item This function first reads in the proteinGroups.txt file
#' produced by MaxQuant.
#' \item It filters out reverse proteins, proteins that were
#' only identified by site, and potential contaminants.
#' \item Then it removes proteins identified with less than
#' the number of unique peptides indicated by \code{uniq.pep} from the the
#' data frame.
#' \item Next, it extracts the LFQ intensity columns and the selected protein
#' rows from the data  frame.
#' \item Converts missing values (zeros) to NAs.
#' \item Finally, the function log transforms the LFQ intensity values.}
#'
#' @return A \code{raw.df} object which is a data frame containing selected
#' proteins as rows and sample LFQ intensities as columns.
#'
#' @examples
#' \dontrun{
#' raw <- create.df(file.path = "./proteinGroups.txt")
#' raw_1 <- create.df(file.path = "./proteinGroups.txt", uniq.pep = 4)
#' }
#' @export

create.df <- function(file.path,
                      filter.prot = TRUE,
                      uniq.pep = 2,
                      zero.NA = TRUE,
                      log.tr = TRUE,
                      base = 2){
  #Load the data
  df = read.csv(file.path,
                sep = "\t",
                stringsAsFactors = FALSE)


  #Filter out some proteins
  if (filter.prot == TRUE){
    df = subset(df,
                df$Only.identified.by.site !="+" &
                  df$Reverse !="+" & df$Potential.contaminant !="+" &
                  df$Unique.peptides > uniq.pep)
  }else{
    warning("Proteins have not been filtered")
  }

  #Extract majority protein group names
  maj_proteins <- df$Majority.protein.IDs
  samples <- df[ , grepl( "LFQ.intensity" , names(df))]
  df <- as.matrix(samples)
  #remove LFQ intensity part from the column name
  colnames(df) <- gsub("LFQ.intensity.", "",colnames(df))
  #add majority protein names to the matrix
  rownames(df) <- c(maj_proteins)

  #Convert zeros to NA
  if (zero.NA == TRUE){
    #Convert matrix to dataframe and convert zeros to NAs
    df <- as.data.frame(df)
    df[df == 0] <- NA
  }else{
    warning("Zeros have not been converted to NAs in the data frame")
  }

  #log2 transform the data
  if(log.tr == TRUE){
    df<- log(df, base)
  }else{
    warning("Intensities have not been log transformed")
  }
  return(df)
}


