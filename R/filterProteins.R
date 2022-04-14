# Filter proteins ---------------------------------------------------------
#' Filter proteins
#' @author Chathurani Ranathunge
#' @description
#' This function filters out proteins with a "+" sign for "Reverse",
#' "Only identified by site", and "Potential contaminants" columns in the
#' proteinGroups.txt file.
#' It also removes proteins identified with less than two unique peptides.
#' @export
filter.prot <- function(x){
  raw_1 = subset(x,
                x$Only.identified.by.site !="+" &
                x$Reverse !="+" & x$Potential.contaminant !="+" &
                x$Unique.peptides > 2)
  return(raw_1)
}

# Filter missing data -----------------------------------------------------
#' Filter missing data
#' @author Chathurani Ranathunge
#' @description This function removes proteins (rows) that are absent (NA)
#' across all samples and samples (columns) with no protein data
#' @export

filter.NA <- function(x){
    #Remove proteins (rows) with missing values (NA) across all samples
    raw_2 <- x[rowSums(is.na(x)) != ncol(x), ]
    #Remove samples (columns) with missing values (NA) across all proteins
    raw_2<- raw_2[, colSums(is.na(raw_2)) != nrow(raw_2)]
    return(raw_2)
}


# Filter proteins by group level missing data -----------------------------
#' Filter proteins by group level missing data
#' @author Chathurani Ranathunge
#' @description This function filters out proteins with more than 33% missing
#' data in each group.
#' @export

filterby.groupNA <- function(x,
                             set.na = 0.33){

  #Extract group information from sample names in the dataframe x
  group <- factor(c(sapply(strsplit(colnames(x), "_"), getElement,1)))

  #Transpose the data frame. Columns are now proteins and rows are samples.
  x_trans <- as.data.frame(t(x))

  #Add a new column with the group information.
  #Group column is the rightmost column in the data frame.
  x_trans$Group <- group

  #Get the number of proteins
  n_proteins <- ncol(x_trans)-1

  #First calculate the mean  number of missing values for each group.
  #This outputs a data frame of mean NA: groups as rows and proteins as columns.
  df_na<-aggregate(is.na(x_trans[, 1:n_proteins]), list(x_trans$Group), mean)

  #Remove the first column that contains the group name in the data frame
  df_na[,1]<- NULL

  #Make a list of proteins with >33% NA
  rem_prot<-as.list(colnames(Filter(function(y) any(as.numeric(y) > set.na),
                                    df_na)))

  #Subset the dataframe by removing proteins in the list with >33% NA
  raw_3<-x[!rownames(x) %in% rem_prot, ]
  return(raw_3)
}


# Output proteins only present in one group -------------------------------

#' Proteins only found in one group
#' @author Chathurani Ranathunge
#' @description This function outputs a list of proteins absent in one group
#' and present in more than 66% of the samples in the other group.
#' @details Comparisons are limited to two groups.
#' @export
onegroup.only<- function(x,
                         abs.group = "",
                         pres.group = "",
                         set.na= 0.33){

  #Extract group information from sample names in the data frame x
  group <- factor(c(sapply(strsplit(colnames(x), "_"), getElement,1)))

  #Transpose the data frame: columns are now proteins and rows are samples.
  x_trans <- as.data.frame(t(x))

  #Add a new column with the group information.
  #Group column is the rightmost column in the data frame.
  x_trans$Group <- group

  #Number of proteins
  n_proteins <- ncol(x_trans)-1

  #Check if there are proteins with 100% missing data in one group and
  #>66% valid data in the other group
  group_only <- as.data.frame(
    ifelse(
      colSums(is.na(x_trans[x_trans$Group==abs.group, 1:n_proteins]))
      /nrow(x_trans[x_trans$Group==abs.group,])==1 &
      colSums(is.na(x_trans[x_trans$Group==pres.group, 1:n_proteins]))
      /nrow(x_trans[x_trans$Group==pres.group,])<=set.na,
      print("TRUE"), print("FALSE")))

  group_only$proteins <- rownames(group_only)

  # Print out a list of proteins only present in user specified group
  cat(group_only[group_only[1] == "TRUE",2],
      file=paste0("Group_",pres.group,"_only.txt"), sep="\n")

}
