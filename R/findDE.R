
#Author: Chathurani Ranathunge
#email: ranathca@evms.edu

# Identify differentially expressed proteins with LIMMA -------------------

find.DEP <- function(x, 
                     save.output = TRUE,
                     filename = "LIMMA_output",
                     save.tophits = TRUE, 
                     ntop = 20){
  
group <- factor(c(sapply(strsplit(colnames(x), "_"), 
                         getElement,1)))
  
#create a design based on the cancer stage
design <- model.matrix(~group)
#Fit the model to the protein intensity data based on the experimental design
fit <- limma::lmFit(x, design)
fit <- limma::eBayes(fit, 
                     robust = T, 
                     trend = T)
if (save.output == TRUE){
  #Write the results of the DE analysis to a text file (tab-separated)
  limma::write.fit(fit,  
            file = paste0(filename,".txt"))          
  }

if(save.tophits == TRUE){
  results_DE<- limma::topTable(fit, 
                               coef=2, 
                               adjust.method = "BH",
                               n = length(fit$df.total))
  results_DE<- results_DE[abs(results_DE$logFC) > 1,]
  results_DE <- results_DE[order(results_DE$P.Value, results_DE$adj.P.Val),]
  write.table(results_DE[1: ntop,],
              file = "TopHits.txt",
              sep = "\t")
  print(results_DE[1: ntop,])
}
return(fit)  
}

# Visualize results from DE analysis with volcano plots -------------------

volcano.plot <- function(x, 
                         adjust.method = "BH",
                         cutoff = 0.05, 
                         FC = 1, 
                         line.FC = TRUE,
                         line.P = TRUE, 
                         filename = "Volcano_plot",
                         filetype = "pdf",
                         height = 2, 
                         width = 2, 
                         col.by = "adjP",
                         nsig.col = "black",
                         sig.col = "red",
                         ln.col = "red",
                         label.top = FALSE, 
                         label.n = 10, 
                         dpi = 80,
                         save= FALSE){

#Extract the required data from the fit object to make our own volcano plot
res_DE <- limma::topTable(x, 
                   coef = colnames(x)[2],
                   n = length(x$df.total))

#Add a new column based on adj.Pval significance status
res_DE$DEAP <- res_DE$adj.P.Val < cutoff & abs(res_DE$logFC) > FC
#Add a new column based on P value significance status alone
res_DE$DEP <- res_DE$P.Value < cutoff & abs(res_DE$logFC) > FC
  
  if(col.by == "P"){
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
    ggplot2::geom_vline(xintercept = c(-1, 1),
                        colour = ln.col, 
                        linetype = 2, 
                        size = 0.2, 
                        alpha = 0.5)  
  }


if(line.P == TRUE){
  DE_volcanoplot <- DE_volcanoplot + 
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        colour = ln.col, 
                        linetype = 2, 
                        size = 0.2, 
                        alpha = 0.5)
  }

if(label.top == TRUE){
  DE_volcanoplot <- DE_volcanoplot + 
    ggplot2::geom_text(data = res_DE[1:label.n,],
                       label = sapply(strsplit(
                         rownames(res_DE[1:label.n,]), ";"),getElement,1 ), 
                       hjust="inward",
                       vjust="outward",
                       check_overlap = FALSE, 
                       size = 0.8) 
  }


if (save == TRUE){
    ggplot2::ggsave(paste0(filename,".", filetype),
                    DE_volcanoplot, 
                    dpi = dpi, 
                    height = height, 
                    width = width)
  return(DE_volcanoplot)
  }else{
    return(DE_volcanoplot)
  }
}

# Visualize DE proteins with a heat map -----------------------------------

heatmap.DE <- function(x, 
                       y,
                       cutoff = 0.05, 
                       FC = 1, 
                       adj.P = 0.05, 
                       sig = "adj.P",
                       palette = "YlGnBu",
                       save = TRUE, 
                       filename = "HeatmapDE",
                       filetype = "pdf",
                       dpi = 80, 
                       height = 7,
                       width = 7){
  
  
#Extract the required data from the fit object
exp_DE <- limma::topTable(x, 
                       coef = colnames(x)[2],
                       n = length(x$df.total),
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
top_intensity <- y[top_proteins,]

#Convert to long format for plotting
top_intMelted <- reshape2::melt(top_intensity)

#add column names
colnames(top_intMelted)<- c("protein", "sample", "Intensity")

#add new column based on cancer stage
top_intMelted$stage <- sapply(strsplit(
  as.character(top_intMelted[,"sample"]), '_'),
  getElement,1)
  
top_heatmap <- ggplot2::ggplot(top_intMelted,
                               ggplot2::aes(x=reorder(sample, -Intensity),
                                            y=reorder(protein, Intensity),
                                            fill=Intensity))+
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
  ggplot2::ggsave(paste0(filename, ".", filetype),
                  top_heatmap, 
                  dpi = dpi, 
                  height = height, 
                  width = width)
  return(top_heatmap)
  }else{
    return(top_heatmap)
  }
}
}


