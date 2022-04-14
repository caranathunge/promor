
#Author: Chathurani Ranathunge
#email: ranathca@evms.edu

#' Load the data -----------------------------------------------------------
#' The protein intensity data is in the proteinGroups.txt file produced by MaxQuant.
#' @export

load.data <- function(file.path = "") {
   input = read.csv(file.path,
              sep = "\t",
              stringsAsFactors = FALSE)
   return(input)
 }
