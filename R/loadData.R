#Load data -----------------------------------------------------------
#' Load input data
#' @description This function loads the proteinGroups.txt file produced by
#' MaxQuant.
#' @author Chathurani Ranathunge
#' @export
load.data <- function(file.path = "") {
   input = read.csv(file.path,
              sep = "\t",
              stringsAsFactors = FALSE)
   return(input)
 }
