
# Scatter plots : correlation between technical replicates ----------------
#' Correlation between technical replicates
#' @author Chathurani Ranathunge
#' @description This function makes scatter plots to visualize the
#' correlation between a given pair of technical replicates (1 vs 2)
#' for each sample.
#' @export
corr.plot<- function(x,
                     rep1,
                     rep2,
                     save = TRUE,
                     filetype= "pdf",
                     col = "red",
                     dot.size= 0.5,
                     text.size = 5,
                     nrow=4,
                     ncol=4,
                     label.size = 1,
                     dpi = 80){

#Separate the sample id from the column names, remove duplicates, paste _ and
#create a list. This is now a list of unique sample names used in the
#dataframe. Eg: _12_ instead of 12 because '12' can match 120, 121, 212 etc.
sample_name <- unlist(paste0("_", unique(sapply(strsplit(colnames(x), "_"),
                              '[',2)),
                              "_", sep=""))

#split the dataframe by each sample and output a list of dataframes for
#plotting. 'as.data.frame' forces even those with just one column or replicate
#into a dataframe rather than a row
sub_df<-lapply(sample_name,
                 function(y) as.data.frame(x[, grepl(y, names(x))]))

#Only keep the samples with at least 2 replicates for plotting TR1 v TR2
if(rep1 == 3  || rep2 == 3){
  plot_data<-Filter(function(c) ncol(c) == 3, sub_df)
}else{
  plot_data<-Filter(function(c) ncol(c) >= 2, sub_df)
}

#Create a list of scatter plots and print/save
plot_list<-lapply(1:length(plot_data), function(t)
  ggplot2::ggplot(plot_data[[t]],
  ggplot2::aes(x=plot_data[[t]][,rep1],y=plot_data[[t]][ ,rep2]))+
  ggplot2::geom_point(col= col, size = dot.size)+
  ggplot2::geom_text(label = sapply(strsplit(rownames(plot_data[[t]]), ";"),
                                    getElement,1),
                                    hjust="inward",
                                    vjust="inward",
                                    size= label.size,
                                    check_overlap = TRUE)+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size = text.size),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank())+
  ggplot2::ggtitle(gsub("\\_0\\d\\_","\\_", colnames(plot_data[[t]][1])))
)
if(save == TRUE){
ggplot2::ggsave(paste0("TR",rep1,"vs","TR",rep2, ".",filetype),
                marrangeGrob(grobs = plot_list, nrow=nrow, ncol=ncol, top=""),
                dpi = dpi)
}else{
  return(plot_list)
}
}


# Remove user-specified samples -------------------------------------------
#' Remove user-specified samples
#' @author Chathurani Ranathunge
#' @description This function removes user-specified samples (columns)
#' from the database
#' @export
rem.sample <- function(x, rem =""){
  x_rem <- x[, -grep(rem, colnames(x))]
  return(x_rem)
}


# Compute average intensity across tech.replicates for each sample --------
#' Compute average intensity
#' @author Chathurani Ranathunge
#' @description This function computes average intensities across
#' technical replicates for each sample.
#' @export
aver.techreps <- function(x){
  #Convert dataframe back to a matrix as dataframes don't allow multiple
  #columns with the same name. We need that feature to average over TRs.
  x_mat<- as.matrix(x)

  #substitute technical replicates with the sample name in the vector.
  #Now each sample, all technical replicates are labelled the same way.
  colnames(x_mat)<-gsub("\\_0\\d\\_","\\_", colnames(x_mat))

  #Average across technical replicates
  x_ave = limma::avearrays(x_mat, ID=colnames(x_mat), weights = NULL)
  return(x_ave)
}
