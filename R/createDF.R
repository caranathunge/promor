
#Author: Chathurani Ranathunge
#email: ranathca@evms.edu

# Create data frame with LFQ intensities ----------------------------------
#This function creates a data frame with LFQ intensities from the filtered
#proteins from filter.prot function. Then converts zeros to NA.

create.df <- function(x, 
                      zeroNA = TRUE,
                      logtr = TRUE,
                      base = 2){
  #Extract majority protein group names
  maj_proteins <- x$Majority.protein.IDs 
  samples <- x[ , grepl( "LFQ.intensity" , names(x))]
  df <- as.matrix(samples)
  #remove LFQ intensity part from the column name
  colnames(df) <- gsub("LFQ.intensity.", "",colnames(df))
  #add majority protein names to the matrix
  rownames(df) <- c(maj_proteins)
  
  #Convert zeros to NA
  if (zeroNA == TRUE){
    #Convert matrix to dataframe and convert zeros to NAs
    df <- as.data.frame(df)
    df[df == 0] <- NA
  }else{
    warning("Zeros have not been converted to NAs in the data frame")
  }
  
  #log2 transform the data
  if(logtr == TRUE){
    df<- log(df, base)
    return(df)
  }else{
    warning("Intensities have not been log transformed")
    return(df)
  }
}  


