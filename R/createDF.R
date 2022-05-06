# Create data frame with LFQ intensities ----------------------------------
#' Create a data frame of protein intensities
#' @description This function creates a data frame of label-free quantitative
#' (LFQ) protein intensities from MaxQuant's proteinGroups.txt file.
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom utils read.csv
#'
#' @param prot.groups File path to proteinGroups.txt file produced by MaxQuant.
#' @param exp.design File path to a text file containing the experimental
#' design.
#' @param filter.n Logical. If \code{TRUE}(default), filters out empty rows and columns
#' from the data frame.
#' @param filter.prot Logical. If \code{TRUE} (default), filters out
#' reverse proteins, proteins only identified by site, potential contaminants,
#' and proteins identified with less than the minimum number of unique peptides
#' indicated by \code{uniq.pep}.
#' @param uniq.pep Numerical. The minimum number of unique peptides required to
#' identify a protein (default is 2). Proteins that are identified by less than
#' this number of unique peptides are filtered out.
#' @param tech.reps Logical. Indicate as \code{TRUE} if technical replicates
#' are present in the data. Defualt is \code{FALSE}.
#' @param zero.na Logical. If \code{TRUE} (default), zeros are considered
#' missing values and replaced with NAs.
#' @param log.tr Logical. If \code{TRUE} (default), intensity values are log
#' transformed to the base indicated by \code{base}.
#' @param base Numerical. Logarithm base. Default is 2.
#'
#' @details
#' \itemize{\item This function first reads in the proteinGroups.txt file
#' produced by MaxQuant.
#' \item It then reads in the expDesign.txt file provided as
#' \code{exp.design} and extracts relevant information from it to add to the
#' data frame.
#' \item First, empty rows and columns are removed from the data frame.
#' \item Next, it filters out reverse proteins, proteins that were
#' only identified by site, and potential contaminants.
#' \item Then it removes proteins identified with less than
#' the number of unique peptides indicated by \code{uniq.pep} from the
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
#' #Generate a raw.df object with default settings.
#' raw <- create_df(file.path = "./proteinGroups.txt",
#'                  exp.design = "./experiment_design.txt")
#' }
#' @export

create_df <- function(prot.groups,
                      exp.design,
                      filter.na = TRUE,
                      filter.prot = TRUE,
                      uniq.pep = 2,
                      tech.reps = FALSE,
                      zero.na = TRUE,
                      log.tr = TRUE,
                      base = 2){
  #Load the data
  df = read.csv(prot.groups,
                sep = "\t",
                stringsAsFactors = FALSE)

  #Load the design file, which is a tab-delimited file containing the
  #experimental design.
  design <- read.csv(exp.design,
                     stringsAsFactors = FALSE,
                     sep ="\t")

  #If tech.reps == FALSE, remove the tech_rep column from the design file and
  #combine remaining columns to make a new sample label
  if(tech.reps == FALSE){
    design$tech_rep <- NULL
    design$new_label <- paste(design$condition,
                              design$sample_ID,
                              sep ="_")
  }else{
    #if tech.reps == TRUE, combine all columns to make new sample label
    design$new_label <- paste(design$condition,
                              design$sample_ID,
                              design$tech_rep,
                              sep ="_")
  }

  #Filter out empty rows and columns if they exist in the dataframe.
  if (filter.na == TRUE){
    #Remove proteins (rows) with missing values (NA) across all samples
    df <- df[rowSums(is.na(df)) != ncol(df), ]
    #Remove samples (columns) with missing values (NA) across all proteins
    df <- df[, colSums(is.na(df)) != nrow(df)]
  }else{
    warning("Data frame contains empty rows and/or columns.")
  }
  #Filter out some proteins based on some columns. First check if the columns
  #are present in the data frame before removing rows based on the presence of
  #"+" signs.

  if (filter.prot == TRUE){
    if("Only.identified.by.site" %in% colnames(df)){
      df = subset(df,
                  df$Only.identified.by.site !="+")
    }
    if("Reverse" %in% colnames(df)){
      df = subset(df,
                  df$Reverse !="+")
    }
    if("Potential.contaminant" %in% colnames(df)){
      df = subset(df,
                  df$Potential.contaminant !="+")
    }
    if("Unique.peptides" %in% colnames(df)){
      df = subset(df,
                  df$Unique.peptides > uniq.pep)
    }
  }else{
    warning("Proteins have not been filtered")
  }

  #Extract majority protein group names
  maj_proteins <- df$Majority.protein.IDs

  #Subset the data frame to only include the LFQ.intensity columns
  samples <- df[ , grepl( "LFQ.intensity" , names(df))]
  df <- as.matrix(samples)

  #remove LFQ intensity part from the column name
  raw_col <- gsub("LFQ.intensity.", "",colnames(df))

  #Compare the mq_label column in the design file with raw_col and replace
  #raw_col with the appropriate new_label.
  raw_col_edited <- mapply(gsub,
                           design$mq_label,
                           design$new_label,
                           raw_col,
                           USE.NAMES = FALSE)

  #add the newly edited column names to the matrix
  colnames(df) <- raw_col_edited

  #add majority protein names to the matrix
  rownames(df) <- c(maj_proteins)

  #Convert zeros to NA
  if (zero.na == TRUE){
    #Convert matrix to dataframe and convert zeros to NAs
    df <- as.data.frame(df)
    df[df == 0] <- NA
  }else{
    warning("Zeros have not been converted to NAs in the data frame")
  }

  #log2 transform the data
  if(log.tr == TRUE){
    df <- log(df, base)
  }else{
    warning("Intensities have not been log transformed")
  }
  return(df)
}


