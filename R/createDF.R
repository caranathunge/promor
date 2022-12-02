# Create data frame with protein intensities ----------------------------------
#' Create a data frame of protein intensities
#' @description This function creates a data frame of protein intensities
#'
#' @author Chathurani Ranathunge
#'
#' @importFrom utils read.csv
#'
#' @param prot_groups File path to a proteinGroups.txt file produced by MaxQuant
#' or a standard input file containing a quantitative matrix
#' where the proteins or protein groups are indicated by rows and the
#' samples by columns.
#' @param exp_design File path to a text file containing the experimental
#' design.
#' @param input_type Type of input file indicated by \code{prot_groups}.
#' Available options are: "MaxQuant", if a proteinGroups.txt file is used, or
#' "standard" if a standard input file is used. Default is "MaxQuant."
#' @param data_type Type of sample protein intensity data columns to use from
#' the proteinGroups.txt file. Some available options are "LFQ", "iBAQ",
#' "Intensity". Default is "LFQ." User-defined prefixes in the proteinGroups.txt
#' file are also allowed. The \code{data_type} argument is case-sensitive, and
#' only applies when \code{input_type = "MaxQuant"}.
#' @param filter_na Logical. If \code{TRUE}(default), filters out empty rows and
#' columns from the data frame.
#' @param filter_prot Logical. If \code{TRUE} (default), filters out
#' reverse proteins, proteins only identified by site, potential contaminants,
#' and proteins identified with less than the minimum number of unique peptides
#' indicated by \code{uniq_pep}. Only applies when
#' \code{input_type = "MaxQuant"}.
#' @param uniq_pep Numerical. The minimum number of unique peptides required to
#' identify a protein (default is 2). Proteins that are identified by less than
#' this number of unique peptides are filtered out. only applies when
#' \code{input_type = "MaxQuant"}.
#' @param tech_reps Logical. Indicate as \code{TRUE} if technical replicates
#' are present in the data. Default is \code{FALSE}.
#' @param zero_na Logical. If \code{TRUE} (default), zeros are considered
#' missing values and replaced with NAs.
#' @param log_tr Logical. If \code{TRUE} (default), intensity values are log
#' transformed to the base indicated by \code{base}.
#' @param base Numerical. Logarithm base. Default is 2.
#'
#' @details
#' \itemize{\item This function first reads in the proteinGroups.txt file
#' produced by MaxQuant or a standard input file containing a quantitative
#' matrix where the proteins or protein groups are indicated by rows and the
#' samples by columns.
#' \item It then reads in the expDesign.txt file provided as
#' \code{exp_design} and extracts relevant information from it to add to the
#' data frame. an example of the expDesign.txt is provided
#' \link[https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt]{here}.
#' \item First, empty rows and columns are removed from the data frame.
#' \item Next, if a proteinGroups.txt file is used, it filters out reverse
#' proteins, proteins that were only identified by site, and potential
#' contaminants.Then it removes proteins identified with less than
#' the number of unique peptides indicated by \code{uniq_pep} from the
#' data frame.
#' \item Next, it extracts the intensity columns indicated by \code{data type}
#' and the selected protein rows from the data frame.
#' \item Converts missing values (zeros) to NAs.
#' \item Finally, the function log transforms the intensity values.}
#'
#' @return A \code{raw_df} object which is a data frame containing protein
#' intensities. Proteins or protein groups are indicated by rows and samples
#' by columns.
#'
#' @examples
#' \donttest{
#'
#' ### Using a proteinGroups.txt file produced by MaxQuant as input.
#' ## Generate a raw_df object with default settings. No technical replicates.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
#'   input_type = "MaxQuant"
#' )
#'
#' ## Data containing technical replicates
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
#'   input_type = "MaxQuant",
#'   tech_reps = TRUE
#' )
#'
#' ## Alter the number of unique peptides needed to retain a protein
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
#'   input_type = "MaxQuant",
#'   uniq_pep = 1
#' )
#'
#' ## Use "iBAQ" values instead of "LFQ" values
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
#'   input_type = "MaxQuant",
#'   data_type = "iBAQ"
#' )
#'
#' ### Using a universal standard input file instead of MaxQuant output.
#' raw_df <- create_df(
#'   prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/st.txt",
#'   exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
#'   input_type = "standard"
#' )
#' }
#' @export

create_df <- function(prot_groups,
                      exp_design,
                      input_type = "MaxQuant",
                      data_type = "LFQ",
                      filter_na = TRUE,
                      filter_prot = TRUE,
                      uniq_pep = 2,
                      tech_reps = FALSE,
                      zero_na = TRUE,
                      log_tr = TRUE,
                      base = 2) {
  # Load the data
  df <- read.csv(prot_groups,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  # Load the design file, which is a tab-delimited file containing the
  # experimental design.
  design <- read.csv(exp_design,
    stringsAsFactors = FALSE,
    sep = "\t"
  )

  #convert data type to lowercase
  input_type <- tolower(input_type)

  #check if the correct input type was entered
  stopifnot("input_type not recognized." = input_type == "maxquant" || input_type == "standard")


  if (tech_reps == TRUE) {
    # if tech_reps == TRUE, combine all columns to make new sample label
    design$new_label <- paste(design$condition,
      design$sample_ID,
      design$tech_rep,
      sep = "_"
    )

    # If tech_reps == FALSE, remove the tech_rep column from the design file and
    # combine remaining columns to make a new sample label
  } else {
    design$tech_rep <- NULL
    design$new_label <- paste(design$condition,
      design$sample_ID,
      sep = "_"
    )
  }
  # extract number of rows
  orig_rows <- nrow(df)
  # extract number of columns
  orig_col <- ncol(df)

  # Filter out empty rows and columns if they exist in the dataframe.
  if (filter_na == TRUE) {
    # Remove proteins (rows) with missing values (NA) across all samples
    df <- df[rowSums(is.na(df)) != ncol(df), ]
    # Remove samples (columns) with missing values (NA) across all proteins
    df <- df[, colSums(is.na(df)) != nrow(df)]
    # calculate number of rows and columns removed
    rem_row <- orig_rows - nrow(df)
    rem_col <- orig_col - ncol(df)
    message(paste0(rem_row, " empty row(s) removed."))
    message(paste0(rem_col, " empty column(s) removed."))
  } else {
    warning("Data frame may contain empty rows and/or columns.")
  }
  # Filter out some proteins based on specific columns. First check if the columns
  # are present in the data frame before removing rows based on the presence of
  # "+" signs.
  # Get number of rows in the df
  orig_rows_1 <- nrow(df)

  if (input_type == "maxquant") {

    #check if the data type is found in the prot_groups file
    stopifnot("data_type not found in prot_groups." = any(grepl(data_type, colnames(df))))

    if (filter_prot == TRUE) {
      if ("Only.identified.by.site" %in% colnames(df)) {
        df <- subset(
          df,
          df$Only.identified.by.site != "+"
        )
        message(paste0(
          orig_rows_1 - nrow(df),
          " protein(s) (rows) only identified by site removed."
        ))
      }

      # get number of rows
      orig_rows_2 <- nrow(df)

      if ("Reverse" %in% colnames(df)) {
        df <- subset(
          df,
          df$Reverse != "+"
        )
        message(paste0(
          orig_rows_2 - nrow(df),
          " reverse protein(s) (rows) removed."
        ))
      }

      # get number of rows
      orig_rows_3 <- nrow(df)

      if ("Potential.contaminant" %in% colnames(df)) {
        df <- subset(
          df,
          df$Potential.contaminant != "+"
        )
        message(paste0(
          orig_rows_3 - nrow(df),
          " protein potential contaminant(s) (rows) removed."
        ))
      }

      if ("Contaminant" %in% colnames(df)) {
        df <- subset(
          df,
          df$Contaminant != "+"
        )
        message(paste0(
          orig_rows_3 - nrow(df),
          " protein contaminant(s) (rows) removed."
        ))
      }

      # get number of rows
      orig_rows_4 <- nrow(df)

      if ("Unique.peptides" %in% colnames(df)) {
        df <- subset(
          df,
          df$Unique.peptides > uniq_pep
        )
        message(paste0(
          orig_rows_4 - nrow(df), " protein(s) identified by ",
          uniq_pep, " or fewer unique peptides removed."
        ))
      }
    } else {
      warning("Proteins have not been filtered")
    }

    # Extract majority protein group names
    maj_proteins <- df$Majority.protein.IDs

    if (data_type == "LFQ") {
      pattern <- paste0(data_type, ".", "intensity", ".", collapse = "")
    } else {
      pattern <- paste0(data_type, ".", collapse = "")
    }

    # Subset the data frame to only include the data_type columns
    samples <- df[, grepl(pattern, colnames(df))]

    # order dataframe columns by column name. Important for next steps involving
    # mapply.
    samples <- samples[, order(colnames(samples))]
    df <- as.matrix(samples)

    # remove data_type part from the column name
    raw_col <- gsub(pattern, "", colnames(df))
  }

  # If a standard table input is used
  if (input_type == "standard") {
    # Extract maj.protein names from the first column
    maj_proteins <- df[, 1]

    # Extract sample names from the column names
    raw_col <- colnames(df)[-1]

    # Create a matrix of intensities, remove the first column
    new_df <- df[, -1]
    df <- as.matrix(new_df)
  }



  # sort the design table by mq_label so that the order matches to that of
  # raw_col
  design <- design[order(as.character(design$mq_label)), ]

  # Compare the mq_label column in the design file with raw_col and replace
  # raw_col with the appropriate new_label.
  raw_col_edited <- mapply(gsub,
    design$mq_label,
    design$new_label,
    raw_col,
    USE.NAMES = FALSE
  )

  # add the newly edited column names to the matrix
  colnames(df) <- raw_col_edited

  # add majority protein names to the matrix
  rownames(df) <- c(maj_proteins)

  # Convert zeros to NA
  if (zero_na == TRUE) {
    # Convert matrix to dataframe and convert zeros to NAs
    df <- as.data.frame(df)
    df[df == 0] <- NA
    message("Zeros have been replaced with NAs.")
  } else {
    warning("Zeros have not been converted to NAs in the data frame")
  }

  # log2 transform the data
  if (log_tr == TRUE) {
    df <- log(df, base)
    message("Data have been log-transformed.")
  } else {
    warning("Intensities have not been log transformed")
  }
  return(df)
}
