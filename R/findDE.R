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
#' @param data A \code{norm.df} object.
#' @param save.output Logical. If \code{TRUE} saves results from the
#' differential expression analysis in a text file labeled "limma_output.txt"
#' in the working directory.
#' @param save.tophits Logical. If \code{TRUE} saves \code{n.top}
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
#' "Benjamini Hochberg" method and outputs the top \code{n.top} results based
#' on lowest p-value and adjusted p-value.
#' \item If the number of hits with absolute log fold change of more than 1 is
#' less than \code{n.top}, \code{find_dep} prints only those with
#' log-fold change > 1 to "TopHits.txt"}
#'
#' @return A \code{fit.df} object, which is similar to a \code{limma}
#' \code{fit} object.
#'
#' @seealso \itemize{\item\code{normalize_data}
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
#' raw_nm <- normalize_data(raw_imp)
#'
#' #Perform differential expression analysis.
#' fit <- find_dep(raw_nm)
#'
#' }
#'
#' @export
find_dep <- function(data,
                     save.output = FALSE,
                     save.tophits = FALSE,
                     n.top = 20){

group <- factor(c(sapply(strsplit(colnames(data), "_"),
                         getElement,1)))

#create a design based on groups
design <- model.matrix(~group)
#Fit the model to the protein intensity data based on the experimental design
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit,
                     robust = T,
                     trend = T)
if (save.output == TRUE){
  #Write the results of the DE analysis to a text file (tab-separated)
  limma::write.fit(fit,
            file = "limma_outout.txt")
  }


  results_DE<- limma::topTable(fit,
                               coef=2,
                               adjust.method = "BH",
                               n = length(fit$df.total))
  results_DE<- results_DE[abs(results_DE$logFC) > 1,]
  results_DE <- results_DE[order(results_DE$P.Value, results_DE$adj.P.Val),]

if(save.tophits == TRUE){
  if(nrow(results_DE) < n.top){
  write.table(results_DE[1: nrow(results_DE),],
              file = "TopHits.txt",
              sep = "\t")
  print(results_DE[1: nrow(results_DE),])
  }else{
  write.table(results_DE[1: n.top,],
              file = "TopHits.txt",
              sep = "\t")
  print(results_DE[1: n.top,])
  }
}else{
  print(results_DE[1:n.top,])
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
#' @import ggrepel
#'
#' @param fit.df A \code{fit.df} object from performing \code{find_dep}.
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
#' @param dot.size Size of the dots/points. Default is \code{0.5}.
#' @param file.name File name to save the plot. Default is "Volcano_plot."
#' @param file.type File type to save the plot. Default is \code{"pdf".}
#' @param plot.height Height of the plot. Default is 7.
#' @param plot.width Width of the plot. Default is 7.
#' @param sig Criteria to denote significance. Choices are \code{"adjP"}
#' (default) for adjusted p-value or \code{"P"} for p-value.
#' @param nsig.col Color of the dots representing significant hits. Default is
#' "grey60."
#' @param sig.col Color of dots representing non-significant hits. Default is
#' "#ffb000."
#' @param ln.col Line color for lines when \code{line.FC = TRUE} and/or
#' \code{line.P = TRUE.} Default is "#419fb7."
#' @param text.size Text size for axis text, labels etc.
#' @param label.top Logical. If \code{TRUE} (default), labels are added to the
#' dots to indicate protein names.
#' @param n.top The number of top hits to label with protein name when
#' \code{label.top = TRUE.} Default is \code{10}.
#' @param dpi Plot resolution. Default is \code{80.}
#' @param save Logical. If \code{TRUE}, saves a copy of the plot in the
#' working directory.
#'
#' @details \itemize{\item Volcano plots show log-2-fold change on the x-axis and
#' -log10(p-value) on the y-axis.\item \code{volcano_plot} requires a
#' \code{fit.df} object from performing differential expression analysis
#' with \code{fit.DEP.}
#' \item User has the option to choose the criterion to denote significance.}
#'
#' @return
#' A \code{ggplot2} plot object.
#'
#' @seealso
#' \itemize{
#' \item \code{find_dep}
#' \item \code{\link[limma:topTable]{topTable}} and
#' \code{\link[limma:lmFit]{lmFit}} functions from the
#' \code{\link[limma]{limma}} package.
#' }
#' @examples
#' \dontrun{
#'
#' #Perform differential expression analysis.
#' fit <- find_dep(raw_nm)
#'
#' #Create a volcano plot with default settings.
#' volcano_plot(fit)
#'
#' #Create a volcano plot with log fold change of 1 and p-value cutoff of 0.05.
#' volcano_plot(fit, cutoff = 0.05, sig = "P", FC = 1)
#'
#' }
#'
#' @export
volcano_plot <- function(fit.df,
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
                         dot.size = 0.5,
                         nsig.col = "grey60",
                         sig.col = "#ffb000",
                         ln.col = "#419fb7",
                         text.size = 10,
                         label.top = FALSE,
                         n.top = 10,
                         dpi = 80,
                         save = FALSE){

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
                          size = dot.size)
    }else{
      res_DE <- res_DE[order(-res_DE$DEAP,res_DE$adj.P.Val),]
      DE_volcanoplot <- ggplot2::ggplot(res_DE,
                                        ggplot2::aes(x = logFC,
                                                     y = -log10(P.Value),
                                                     colour = DEAP)) +
        ggplot2::geom_point(ggplot2::aes(colour = DEAP),
                            alpha = 0.5,
                            size = dot.size)
    }



DE_volcanoplot <- DE_volcanoplot +
    ggplot2::xlab(expression("log" [2]* " fold change"))+
    ggplot2::ylab( expression("-log" [10]* "(P-value)"))+
    ggplot2::scale_colour_manual(values = c(nsig.col, sig.col)) +
  ggplot2::theme_bw()+
  ggplot2::theme(legend.position = "",
          text = element_text(size = text.size),
          axis.line.x = element_line(size = 0.1),
          axis.line.y = element_line(size = 0.1),
          axis.ticks.x = element_line(size = 0.1),
          axis.ticks.y = element_line(size = 0.1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line (size = 0.1),
          panel.border = element_rect(size = 0.2))

if(line.FC == TRUE){
  DE_volcanoplot <- DE_volcanoplot +
    ggplot2::geom_vline(xintercept = c(-FC, FC),
                        colour = ln.col,
                        linetype = 2,
                        size = 0.5,
                        alpha = 0.8)
  }


if(line.P == TRUE){
  DE_volcanoplot <- DE_volcanoplot +
    ggplot2::geom_hline(yintercept = -log10(cutoff),
                        colour = ln.col,
                        linetype = 2,
                        size = 0.5,
                        alpha = 0.8)
  }

if(label.top == TRUE){
  DE_volcanoplot <- DE_volcanoplot +
    ggrepel::geom_text_repel(data = res_DE[1:n.top,],
                       label = sapply(strsplit(
                         rownames(res_DE[1:n.top,]), ";"),
                         getElement,1 ),
                         size = text.size/4)
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
#' @param fit.df A \code{fit.df} object from performing \code{find_dep}.
#' @param norm.df The \code{norm.df} object from which the \code{fit.df} object
#' was obtained.
#' @param cutoff Cutoff value for p-values and adjusted p-values. Default is
#' 0.05.
#' @param FC Minimum absolute log-fold change to use as threshold for
#' differential expression. Default is 1.
#' @param sig Criteria to denote significance. Choices are \code{"adjP"}
#' (default) for adjusted p-value or \code{"P"} for p-value.
#' @param n.top Number of top hits to include in the heat map.
#' @param palette Color palette for plots. Default is \code{"YlGnBu."}
#' @param text.size Text size for axis text, labels etc.
#' @param save Logical. If \code{TRUE}, saves a copy of the plot in
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
#' \item \code{find_dep}
#' \item \code{\link[limma:topTable]{topTable}} and
#' \code{\link[limma:lmFit]{lmFit}} functions from the
#' \code{\link[limma]{limma}} package.}
#'
#' @examples
#' \dontrun{
#' #Perform differential expression analysis.
#' fit <- find_dep(raw_nm)
#'
#' #Create a heatmap with default settings.
#' heatmap_de(fit, raw_nm)
#'
#' #Create a heatmap with log fold change of 1 and p-value cutoff of 0.05.
#' heatmap_de(fit, raw_nm, cutoff = 0.05, sig = "P", FC = 1)
#'
#' }
#'
#'
#' @export
heatmap_de <- function(fit.df,
                       norm.df ,
                       cutoff = 0.05,
                       FC = 1,
                       sig = "adjP",
                       n.top = 20,
                       palette = "YlGnBu",
                       text.size = 10,
                       save = FALSE,
                       file.name = "HeatmapDE",
                       file.type = "pdf",
                       dpi = 80,
                       plot.height  = 7,
                       plot.width = 7){

#Binding the global variables to a local function
logFC <- P.Value <- adj.P.Val <- Intensity <- protein <- NULL

#Extract the required data from the fit object
exp_DE <- limma::topTable(fit.df,
                       coef = colnames(fit.df)[2],
                       n = length(fit.df$df.total),
                       adjust.method = "BH")

#Pick the sig. proteins based on lowest p-value and highest logFC

if(sig == "P"){
  top_proteins <- rownames(subset(exp_DE,
                                  abs(logFC) > FC & P.Value < cutoff,
                                 drop = FALSE))

#Or default: based on adj.P value
}else{
  top_proteins <- rownames(subset(exp_DE,
                                abs(logFC) > FC & adj.P.Val < cutoff,
                                drop = FALSE))
}

#Extract the top n.top hits from the top hit list
top_proteins <- top_proteins[1:n.top]

#Check if there are sig. proteins before moving on to plotting
if (identical(top_proteins, character(0))){
  stop(message
       (paste0("No significant proteins found at ", sig, " < ", cutoff, "." )))
}else{

#Extract intensity values for top proteins based on logFC and p-value cutoff
top_intensity <- subset(norm.df,
                      rownames(norm.df) %in% top_proteins,
                      drop = FALSE)
#Convert to long format for plotting
top_intMelted <- reshape2::melt(top_intensity)

#add column names
colnames(top_intMelted)<- c("protein", "sample", "Intensity")

#add new column based on cancer stage
top_intMelted$stage <- sapply(strsplit(
  as.character(top_intMelted[,"sample"]), '_'),
  getElement,1)

top_heatmap <- ggplot2::ggplot(top_intMelted,
                               #ggplot2::aes(x = reorder(sample, -Intensity),
                                            #y = reorder(protein, Intensity),
                                            #fill = Intensity))+
                               ggplot2::aes(x = sample,
                                            y = protein,
                                            fill = Intensity))+
  ggplot2::geom_tile(colour="white", size=0.1)+
  ggplot2::labs(x="", y="")+
  ggplot2::scale_y_discrete(expand = c(0, 0),
                            labels  = sapply(
                              strsplit(
                                as.character(top_intMelted[,"protein"]),';'),
                              getElement,1))+
  #ggplot2::scale_fill_gradientn(colors = rev(hcl.colors(10, palette))) +
  ggplot2::scale_fill_distiller(palette = palette, direction = 1)+
  ggplot2::theme_grey(base_size = 8)+
  ggplot2::theme(aspect.ratio = 1,
                 legend.position="right",
                 legend.direction="vertical",
                 legend.key.height = grid::unit(0.8, "cm"),
                 legend.key.width = grid::unit(0.2, "cm"),
                 legend.title = element_blank(),
                 legend.text = element_text(size = text.size*0.7,
                                            face = "bold"),
                 axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size = text.size * 0.7),
                 strip.background = element_blank(),
                 strip.text = element_text(size = text.size,
                                           hjust = 0.01,
                                           face = "bold",
                                           vjust = 0 ),
                 plot.background = element_blank(),
                 panel.border = element_blank(),
                 panel.spacing = unit(text.size * 0.001, "cm"),
                 panel.background = element_blank())+
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


