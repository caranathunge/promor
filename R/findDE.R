# Identify differentially expressed proteins with limma -------------------
#' Identify differentially expressed proteins between groups
#' @author Chathurani Ranathunge
#' @description This function performs differential expression analysis
#' on protein intensity data with limma.
#'
#' @import limma
#' @importFrom stats model.matrix
#' @importFrom utils write.table
#'
#' @param df A \code{norm.df} object.
#' @param save.output Logical. If \code{TRUE} (default) saves results from the
#' differential expression analysis in a text file labeled "limma_output.txt"
#' in the working directory.
#' @param save.tophits Logical. If \code{TRUE} (default) saves \code{n.top}
#' number of top hits from the differential expression analysis in a text file
#' labeled "TopHits.txt" in the working directory.
#' @param n.top The number of top differentially expressed proteins to save in
#' the "TopHits.txt" file. Default is \code{20}.
#'
#' @details \itemize{\item
#' It is important that the data is first log-transformed, ideally,
#' imputed, and normalized before performing differential expression analysis.
#' \item \code{save.output} saves the complete results table from the
#' differential expression analysis.
#' \item \code{save.tophits} first subsets the results to those with absolute
#' log fold change of more than 1, performs multiple correction with
#' "Benjamini Hochberg" method and outputs the top results based on lowest
#' p-value and adjusted p-value.}
#'
#' @return A \code{fit.df} object, which is similar to a \code{limma}
#' \code{fit} object.
#'
#' @seealso \itemize{\item\code{normalize.data}
#' \item\code{\link[limma:lmFit]{lmFit}},
#' \code{\link[limma:eBayes]{eBayes}},
#' \code{\link[limma:topTable]{topTable}}, and
#' \code{\link[limma:write.fit]{write.fit}} functions from the
#' \code{\link[limma]{limma}} package.}
#'
#' @examples
#' \dontrun{
#'
#' #Normalize an already imputed data set.
#' raw_nm <- normalize.data(raw_imp)
#'
#' #Perform differential expression analysis.
#' fit <- find.DEP(raw_nm)
#'
#' }
#'
#' @export
find.DEP <- function(df,
                     save.output = TRUE,
                     save.tophits = TRUE,
                     n.top = 20){

group <- factor(c(sapply(strsplit(colnames(df), "_"),
                         getElement,1)))

#create a design based on groups
design <- model.matrix(~group)
#Fit the model to the protein intensity data based on the experimental design
fit <- limma::lmFit(df, design)
fit <- limma::eBayes(fit,
                     robust = T,
                     trend = T)
if (save.output == TRUE){
  #Write the results of the DE analysis to a text file (tab-separated)
  limma::write.fit(fit,
            file = "limma_outout.txt")
  }

if(save.tophits == TRUE){
  results_DE<- limma::topTable(fit,
                               coef=2,
                               adjust.method = "BH",
                               n = length(fit$df.total))
  results_DE<- results_DE[abs(results_DE$logFC) > 1,]
  results_DE <- results_DE[order(results_DE$P.Value, results_DE$adj.P.Val),]
  write.table(results_DE[1: n.top,],
              file = "TopHits.txt",
              sep = "\t")
  print(results_DE[1: n.top,])
}
return(fit)
}

# Visualize differentially expressed proteins with volcano plots -------------------
#' Volcano plot
#' @author Chathurani Ranathunge
#' @description This function generates volcano plots to visualize
#' differentially expressed proteins between groups.
#' @import limma
#' @import ggplot2
#'
#' @param fit.df A \code{fit.df} object from performing \code{find.DEP}.
#' @param adj.method Method used for adjusting the p-values for multiple
#' testing. Default is \code{"BH."}.
#' @param cutoff Cutoff value for p-values and adjusted p-values. Default is
#' 0.05.
#' @param FC Minimum absolute log-fold change to use as threshold for
#' differential expression.
#' @param line.FC Logical. If \code{TRUE}(default), a dotted line will be shown
#' to indicate the \code{FC} threshold in the plot.
#' @param line.P Logical. If \code{TRUE}(default), a dotted line will be shown
#' to indicate the p-value \code{cutoff.}
#' @param file.name File name to save the plot. Default is "Volcano_plot."
#' @param file.type File type to save the plot. Default is \code{"pdf".}
#' @param plot.height Height of the plot. Default is 7.
#' @param plot.width Width of the plot. Default is 7.
#' @param sig Criteria to denote significance. Choices are \code{"adjP"}
#' (default) for adjusted p-value or \code{"P"} for p-value.
#' @param nsig.col Color of the dots representing significant hits. Default is
#' "Red."
#' @param sig.col Color of dots representing non-significant hits. Default is
#' "Black."
#' @param ln.col Line color for lines when \code{line.FC = TRUE} and/or
#' \code{line.P = TRUE.} Default is "Black."
#' @param label.top Logical. If \code{TRUE} (default), labels are added to the
#' dots to indicate protein names.
#' @param n.top The number of top hits to label with protein name when
#' \code{label.top = TRUE.} Default is \code{10}.
#' @param dpi Plot resolution. Default is \code{80.}
#' @param save Logical. If \code{TRUE} (default), saves a copy of the plot in the
#' working directory.
#'
#' @details \itemize{\item Volcano plots show log-2-fold change on the x-axis and
#' -log10(p-value) on the y-axis.\item \code{volcano.plot} requires a
#' \code{fit.df} object from performing differential expression analysis
#' with \code{fit.DEP.}
#' \item User has the option to choose the criterion to denote significance.}
#'
#' @return
#' A \code{ggplot2} plot object.
#'
#' @seealso
#' \itemize{
#' \item \code{find.DEP}
#' \item \code{\link[limma:topTable]{topTable}} and
#' \code{\link[limma:lmFit]{lmFit}} functions from the
#' \code{\link[limma]{limma}} package.
#' }
#' @examples
#' \dontrun{
#'
#' #Perform differential expression analysis.
#' fit <- find.DEP(raw_nm)
#'
#' #Create a volcano plot with default settings.
#' volcano.plot(fit)
#'
#' #Create a volcano plot with log fold change of 1 and p-value cutoff of 0.05.
#' volcanoplot(fit, cutoff = 0.05, sig = "P", FC = 1)
#'
#' }
#'
#' @export
volcano.plot <- function(fit.df,
                         adj.method = "BH",
                         cutoff = 0.05,
                         FC = 1,
                         line.FC = TRUE,
                         line.P = TRUE,
                         file.name = "Volcano_plot",
                         file.type = "pdf",
                         plot.height = 7,
                         plot.width = 7,
                         sig = "adjP",
                         nsig.col = "black",
                         sig.col = "red",
                         ln.col = "red",
                         label.top = FALSE,
                         n.top = 10,
                         dpi = 80,
                         save = TRUE){

#Set global variables to NULL
logFC <- P.Value <- DEP <- DEAP <- NULL


#Extract the required data from the fit object to make our own volcano plot
res_DE <- limma::topTable(fit.df,
                          adjust.method =  adj.method,
                          coef = colnames(fit.df)[2],
                          n = length(fit.df$df.total))

#Add a new column based on adj.Pval significance status
res_DE$DEAP <- res_DE$adj.P.Val < cutoff & abs(res_DE$logFC) > FC
#Add a new column based on P value significance status alone
res_DE$DEP <- res_DE$P.Value < cutoff & abs(res_DE$logFC) > FC

  if(sig == "P"){
    res_DE <- res_DE[order(-res_DE$DEP,res_DE$adj.P.Val),]
    DE_volcanoplot <- ggplot2::ggplot(res_DE,
                                      ggplot2::aes(x = logFC,
                                                   y = -log10(P.Value),
                                                   colour = DEP)) +
      ggplot2::geom_point(ggplot2::aes(colour = DEP),
                          alpha = 0.5,
                          size = 0.1)
    }else{
      res_DE <- res_DE[order(-res_DE$DEAP,res_DE$adj.P.Val),]
      DE_volcanoplot <- ggplot2::ggplot(res_DE,
                                        ggplot2::aes(x = logFC,
                                                     y = -log10(P.Value),
                                                     colour = DEAP)) +
        ggplot2::geom_point(ggplot2::aes(colour = DEAP),
                            alpha = 0.5,
                            size = 0.1)
    }



DE_volcanoplot <- DE_volcanoplot +
    ggplot2::xlab(expression("log" [2]* " fold change"))+
    ggplot2::ylab( expression("-log" [10]* "(P-value)"))+
    ggplot2::scale_colour_manual(values = c(nsig.col, sig.col)) +
  ggplot2::theme_bw()+
  ggplot2::theme(legend.position = "",
          text = element_text(size = 3),
          axis.line.x = element_line(size = 0.1),
          axis.line.y = element_line(size = 0.1),
          axis.ticks.x = element_line(size = 0.1),
          axis.ticks.y = element_line(size = 0.1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line (size = 0.1),
          panel.border = element_rect (size = 0.2))

if(line.FC == TRUE){
  DE_volcanoplot <- DE_volcanoplot +
    ggplot2::geom_vline(xintercept = c(-FC, FC),
                        colour = ln.col,
                        linetype = 2,
                        size = 0.2,
                        alpha = 0.5)
  }


if(line.P == TRUE){
  DE_volcanoplot <- DE_volcanoplot +
    ggplot2::geom_hline(yintercept = -log10(cutoff),
                        colour = ln.col,
                        linetype = 2,
                        size = 0.2,
                        alpha = 0.5)
  }

if(label.top == TRUE){
  DE_volcanoplot <- DE_volcanoplot +
    ggplot2::geom_text(data = res_DE[1:n.top,],
                       label = sapply(strsplit(
                         rownames(res_DE[1:n.top,]), ";"),getElement,1 ),
                       hjust="inward",
                       vjust="outward",
                       check_overlap = FALSE,
                       size = 0.8)
  }


if (save == TRUE){
    ggplot2::ggsave(paste0(file.name,".", file.type),
                    DE_volcanoplot,
                    dpi = dpi,
                    height = plot.height,
                    width = plot.width)
  return(DE_volcanoplot)
  }else{
    return(DE_volcanoplot)
  }
}

# Visualize DE proteins with a heat map -----------------------------------
#' Heatmap of differentially expressed proteins
#' @description This function generates a heatmap to visualize differentially
#' expressed proteins between groups
#'
#' @author Chathurani Ranathunge
#'
#' @import limma
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices hcl.colors
#' @importFrom stats reorder
#'
#' @param fit.df A \code{fit.df} object from performing \code{find.DEP}.
#' @param norm.df The \code{norm.df} object from which the \code{fit.df} object
#' was obtained.
#' @param cutoff Cutoff value for p-values and adjusted p-values. Default is
#' 0.05.
#' @param FC Minimum absolute log-fold change to use as threshold for
#' differential expression.
#' @param sig Criteria to denote significance. Choices are \code{"adjP"}
#' (default) for adjusted p-value or \code{"P"} for p-value.
#' @param palette Color palette for plots. Default is \code{"YlGnBu."}
#' @param save Logical. If \code{TRUE} (default), saves a copy of the plot in
#' the working directory.
#' @param file.name File name to save the plot. Default is "HeatmapDE."
#' @param file.type File type to save the plot. Default is \code{"pdf".}
#' @param dpi Plot resolution. Default is \code{80.}
#' @param plot.height Height of the plot. Default is 7.
#' @param plot.width Width of the plot. Default is 7.
#'
#' @details
#' By default the tiles in the heatmap are reordered by intensity values
#' along both axes (x axis = samples, y axis = proteins).
#'
#' @return A \code{ggplot2} plot object.
#'
#' @seealso
#' \itemize{
#' \item \code{find.DEP}
#' \item \code{\link[limma:topTable]{topTable}} and
#' \code{\link[limma:lmFit]{lmFit}} functions from the
#' \code{\link[limma]{limma}} package.}
#'
#' @examples
#' \dontrun{
#' #Perform differential expression analysis.
#' fit <- find.DEP(raw_nm)
#'
#' #Create a heatmap with default settings.
#' heatmap.DE(fit, raw_nm)
#'
#' #Create a heatmap with log fold change of 1 and p-value cutoff of 0.05.
#' heatmap.DE(fit, raw_nm, cutoff = 0.05, sig = "P", FC = 1)
#'
#' }
#'
#'
#' @export
heatmap.DE <- function(fit.df,
                       norm.df ,
                       cutoff = 0.05,
                       FC = 1,
                       sig = "adjP",
                       palette = "YlGnBu",
                       save = TRUE,
                       file.name = "HeatmapDE",
                       file.type = "pdf",
                       dpi = 80,
                       plot.height  = 7,
                       plot.width = 7){

#Binding the gloabl variables to a local function
logFC <- P.Value <- adj.P.Val <- Intensity <- protein <- NULL

#Extract the required data from the fit object
exp_DE <- limma::topTable(fit.df,
                       coef = colnames(fit.df)[2],
                       n = length(fit.df$df.total),
                       adjust.method = "BH")


#Pick the sig. proteins based on lowest p-value and highest logFC

if(sig == "P"){
  top_proteins <- rownames(subset(exp_DE,
                                  abs(logFC) > FC & P.Value < cutoff))
}else{
  top_proteins <- rownames(subset(exp_DE,
                                abs(logFC) > FC & adj.P.Val < cutoff))
}

#Check if there are sig. proteins before moving on to plotting
if (identical(top_proteins, character(0))){
  stop(message
       (paste0("No significant proteins found at ", sig, " < ", cutoff, "." )))
}else{

  #Extract intensity values for top proteins based on logFC and p-value cutoff
top_intensity <- norm.df[top_proteins,]

#Convert to long format for plotting
top_intMelted <- reshape2::melt(top_intensity)

#add column names
colnames(top_intMelted)<- c("protein", "sample", "Intensity")

#add new column based on cancer stage
top_intMelted$stage <- sapply(strsplit(
  as.character(top_intMelted[,"sample"]), '_'),
  getElement,1)

top_heatmap <- ggplot2::ggplot(top_intMelted,
                               ggplot2::aes(x = reorder(sample, -Intensity),
                                            y = reorder(protein, Intensity),
                                            fill = Intensity))+
  ggplot2::geom_tile(colour="white", size=0.1)+
  ggplot2::labs(x="", y="")+
  ggplot2::scale_y_discrete(expand=c(0, 0),
                            labels  = sapply(
                              strsplit(
                                as.character(top_intMelted[,"protein"]),';'),
                              getElement,1))+
  ggplot2::scale_fill_gradientn(colors = rev(hcl.colors(50, palette))) +
  ggplot2::theme_grey(base_size = 8)+
  ggplot2::theme(aspect.ratio = 1,
                 legend.position="right",
                 legend.direction="vertical",
                 legend.key.height = grid::unit(0.8, "cm"),
                 legend.key.width = grid::unit(0.2, "cm"),
                 legend.text = element_text(face = "bold"),
                 axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(),
                 strip.background = element_rect(fill = "grey97",
                                                 colour = "grey90",
                                                 size = 0.3),
                 plot.background = element_blank(),
                 panel.border = element_blank())+
    ggplot2::facet_wrap(.~stage,
                        scales = "free_x")

if (save == TRUE){
  ggplot2::ggsave(paste0(file.name, ".", file.type),
                  top_heatmap,
                  dpi = dpi,
                  height = plot.height,
                  width = plot.width)
  return(top_heatmap)
  }else{
    return(top_heatmap)
  }
}
}


