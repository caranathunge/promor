# Identify differentially expressed proteins with limma -------------------
#' Identify differentially expressed proteins between groups
#' @author Chathurani Ranathunge
#' @description This function performs differential expression analysis
#' on protein intensity data with limma.
#'
#' @import limma
#' @import statmod
#' @importFrom stats model.matrix
#' @importFrom utils write.table
#'
#' @param df A \code{norm_df} object or an \code{imp_df} object.
#' @param save_output Logical. If \code{TRUE} saves results from the
#' differential expression analysis in a text file labeled "limma_output.txt"
#' in the directory specified by \code{file_path}.
#' @param save_tophits Logical. If \code{TRUE} saves \code{n_top}
#' number of top hits from the differential expression analysis in a text file
#' labeled "TopHits.txt" in the directory specified by \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param lfc Minimum absolute log2-fold change to use as threshold for
#' differential expression.
#' @param adj_method Method used for adjusting the p-values for multiple
#' testing. Default is \code{"BH"} for "Benjamini-Hochberg" method.
#' @param cutoff Cutoff value for p-values and adjusted p-values. Default is
#' 0.05.
#' @param n_top The number of top differentially expressed proteins to save in
#' the "TopHits.txt" file. Default is \code{20}.
#'
#' @details \itemize{\item
#' It is important that the data is first log-transformed, ideally,
#' imputed, and normalized before performing differential expression analysis.
#' \item \code{save_output} saves the complete results table from the
#' differential expression analysis.
#' \item \code{save_tophits} first subsets the results to those with absolute
#' log fold change of more than 1, performs multiple correction with
#' the method specified in \code{adj_method} and outputs the top \code{n_top}
#' results based on lowest p-value and adjusted p-value.
#' \item If the number of hits with absolute log fold change of more than 1 is
#' less than \code{n_top}, \code{find_dep} prints only those with
#' log-fold change > 1 to "TopHits.txt".
#' \item If the \code{file_path} is not specified, text files will be saved in
#' a temporary directory.}
#'
#' @return A \code{fit_df} object, which is similar to a \code{limma}
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
#'
#' ## Perform differential expression analysis using default settings
#' fit_df1 <- find_dep(ecoli_norm_df)
#'
#' ## Change p-value and adjusted p-value cutoff
#' fit_df2 <- find_dep(ecoli_norm_df, cutoff = 0.1)
#'
#' @references Ritchie, Matthew E., et al. "limma powers differential expression
#' analyses for RNA-sequencing and microarray studies." Nucleic acids research
#' 43.7 (2015): e47-e47.
#'
#' @export
find_dep <- function(df,
                     save_output = FALSE,
                     save_tophits = FALSE,
                     file_path = NULL,
                     adj_method = "BH",
                     cutoff = 0.05,
                     lfc = 1,
                     n_top = 20) {
  # Extract group information from colnames
  group <- factor(c(sapply(
    strsplit(colnames(df), "_"),
    getElement, 1
  )))

  # create a design based on groups
  design <- model.matrix(~group)

  # Fit the model to the protein intensity data based on the experimental design
  fit <- limma::lmFit(df, design)
  fit <- limma::eBayes(fit,
    robust = T,
    trend = T
  )

  # Make a a list of DE results based on provided criteria
  dec_test <- limma::decideTests(fit,
    lfc = lfc,
    adjust.method = adj_method
  )

  # Set temporary file_path if not specified
  if (is.null(file_path)) {
    file_path <- tempdir()
  }

  # Write the results of the DE analysis to a text file (tab-separated)
  if (save_output == TRUE) {
    limma::write.fit(fit,
      file = paste0(file_path, "/limma_outout.txt"),
      adjust = adj_method,
      results = dec_test
    )
  }


  results_de <- limma::topTable(fit,
    coef = 2,
    adjust.method = adj_method,
    n = Inf
  )

  # add majority protein ids column
  results_de$majority_protein_id <- rownames(results_de)
  rownames(results_de) <- NULL

  # rearrange order of columns
  results_de <- results_de[, c(
    "majority_protein_id",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )]

  # extract proteins with absolute logfc > lfc
  results_de <- results_de[abs(results_de$logFC) > lfc, ]

  # extract sig. de. proteins and order from smallest to largest p. values
  results_de <- results_de[results_de$adj.P.Val < cutoff, ]
  results_de <- results_de[order(results_de$P.Value, results_de$adj.P.Val), ]

  if (nrow(results_de) == 0) {
    stop(
      message(
        paste0(
          "No differentially expressed proteins found at adj.P.value cutoff = ",
          cutoff
        )
      )
    )
  } else {
    message(paste0(
      nrow(results_de),
      " siginificantly differentially expressed proteins found."
    ))
  }

  if (save_tophits == TRUE) {
    if (nrow(results_de) < n_top) {
      write.table(results_de[seq_len(nrow(results_de)), ],
        file = paste0(file_path, "/TopHits.txt"),
        sep = "\t",
        quote = FALSE
      )
    } else {
      write.table(results_de[1:n_top, ],
        file = paste0(file_path, "/TopHits.txt"),
        sep = "\t",
        quote = FALSE
      )
    }
  }
  return(fit)
}

# Visualize differentially expressed proteins with volcano plots --------------
#' Volcano plot
#' @author Chathurani Ranathunge
#' @description This function generates volcano plots to visualize
#' differentially expressed proteins between groups.
#' @import limma
#' @import ggplot2
#' @import ggrepel
#' @import viridis
#'
#' @param fit_df A \code{fit_df} object from performing \code{find_dep}.
#' @param adj_method Method used for adjusting the p-values for multiple
#' testing. Default is \code{"BH"}.
#' @param cutoff Cutoff value for p-values and adjusted p-values. Default is
#' 0.05.
#' @param lfc Minimum absolute log2-fold change to use as threshold for
#' differential expression.
#' @param line_fc Logical. If \code{TRUE}(default), a dotted line will be shown
#' to indicate the \code{lfc} threshold in the plot.
#' @param line_p Logical. If \code{TRUE}(default), a dotted line will be shown
#' to indicate the p-value or adjusted p-value \code{cutoff.}
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis:viridis]{viridis}}
#' for available options.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' directory provided in \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param file_name File name to save the plot. Default is "Volcano_plot."
#' @param file_type File type to save the plot. Default is \code{"pdf".}
#' @param plot_height Height of the plot. Default is 7.
#' @param plot_width Width of the plot. Default is 7.
#' @param sig Criteria to denote significance. Choices are \code{"adjP"}
#' (default) for adjusted p-value or \code{"P"} for p-value.
#' @param text_size Text size for axis text, labels etc.
#' @param label_top Logical. If \code{TRUE} (default), labels are added to the
#' dots to indicate protein names.
#' @param n_top The number of top hits to label with protein name when
#' \code{label_top = TRUE.} Default is \code{10}.
#' @param dpi Plot resolution. Default is \code{80.}
#'
#' @details \itemize{\item Volcano plots show log-2-fold change on the x-axis,
#'  and based on the significance criteria chosen, either -log10(p-value) or
#'  -log10(adjusted p-value) on the y-axis.
#' \item \code{volcano_plot} requires a
#' \code{fit_df} object from performing differential expression analysis
#' with \code{find_dep.}
#' \item User has the option to choose criteria that denote significance.}
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
#'
#' ## Create a volcano plot with default settings.
#' volcano_plot(ecoli_fit_df)
#'
#' ## Change significance criteria and cutoff
#' volcano_plot(ecoli_fit_df, cutoff = 0.1, sig = "P")
#'
#' ## Label top 30 differentially expressed proteins and
#' ## change the color palette of the plot
#' volcano_plot(ecoli_fit_df, label_top = TRUE, n_top = 30, palette = "mako")
#'
#' @export
volcano_plot <- function(fit_df,
                         adj_method = "BH",
                         sig = "adjP",
                         cutoff = 0.05,
                         lfc = 1,
                         line_fc = TRUE,
                         line_p = TRUE,
                         palette = "viridis",
                         text_size = 10,
                         label_top = FALSE,
                         n_top = 10,
                         save = FALSE,
                         file_path = NULL,
                         file_name = "Volcano_plot",
                         file_type = "pdf",
                         plot_height = 7,
                         plot_width = 7,
                         dpi = 80) {
  # Set global variables to NULL
  logFC <- P.Value <- dep <- de_ap <- NULL


  # Extract the required data from the fit object to make our own volcano plot
  res_de <- limma::topTable(fit_df,
    adjust.method = adj_method,
    coef = colnames(fit_df)[2],
    n = length(fit_df$df.total)
  )

  # Add a new column based on adj.Pval significance status
  res_de$de_ap <- res_de$adj.P.Val < cutoff & abs(res_de$logFC) > lfc
  # Add a new column based on P value significance status alone
  res_de$dep <- res_de$P.Value < cutoff & abs(res_de$logFC) > lfc

  if (sig == "P") {
    res_de <- res_de[order(-res_de$dep, res_de$adj.P.Val), ]
    de_volcanoplot <- ggplot2::ggplot(
      res_de,
      aes(
        x = logFC,
        y = -log10(P.Value),
        color = dep
      )
    ) +
      ggplot2::geom_point(aes(color = dep),
        alpha = 0.7,
        size = text_size * 0.3
      )+
      ggplot2::ylab(expression("-log"[10] * "(P-value)"))

  } else {
    res_de <- res_de[order(-res_de$de_ap, res_de$adj.P.Val), ]
    de_volcanoplot <- ggplot2::ggplot(
      res_de,
      aes(
        x = logFC,
        y = -log10(adj.P.Val),
        color = de_ap
      )
    ) +
      ggplot2::geom_point(aes(color = de_ap),
        alpha = 0.7,
        size = text_size * 0.3
      )+
      ggplot2::ylab(expression("-log"[10] * "(adj.P-value)"))
  }



  de_volcanoplot <- de_volcanoplot +
    ggplot2::xlab(expression("log"[2] * " fold change")) +
    viridis::scale_color_viridis(
      discrete = TRUE,
      direction = 1,
      option = palette,
      begin = 0.2,
      end = 0.8
    ) +
    promor_theme() +
    ggplot2::theme(
      legend.position = "",
      panel.grid.major = element_line(
        size = 0.1,
        color = "grey80"
      )
    )

  if (line_fc == TRUE) {
    de_volcanoplot <- de_volcanoplot +
      ggplot2::geom_vline(
        xintercept = c(-lfc, lfc),
        color = "grey60",
        linetype = 2,
        size = 0.5,
        alpha = 0.8
      )
  }


  if (line_p == TRUE) {
    de_volcanoplot <- de_volcanoplot +
      ggplot2::geom_hline(
        yintercept = -log10(cutoff),
        color = "grey60",
        linetype = 2,
        size = 0.5,
        alpha = 0.8
      )
  }

  if (label_top == TRUE) {
    de_volcanoplot <- de_volcanoplot +
      ggrepel::geom_text_repel(
        data = res_de[1:n_top, ],
        label = sapply(
          strsplit(
            rownames(res_de[1:n_top, ]), ";"
          ),
          getElement, 1
        ),
        size = text_size / 4
      )
  }

  # Set temporary file_path if not specified
  if (is.null(file_path)) {
    file_path <- tempdir()
  }

  if (save == TRUE) {
    ggplot2::ggsave(paste0(file_path, "/", file_name, ".", file_type),
      de_volcanoplot,
      dpi = dpi,
      height = plot_height,
      width = plot_width
    )
    return(de_volcanoplot)
  } else {
    return(de_volcanoplot)
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
#' @importFrom stats reorder
#' @import viridis
#'
#' @param fit_df A \code{fit_df} object from performing \code{find_dep}.
#' @param df The \code{norm_df} object or the \code{imp_df} object from which
#' the \code{fit_df} object was obtained.
#' @param adj_method Method used for adjusting the p-values for multiple
#' testing. Default is \code{"BH"}.
#' @param cutoff Cutoff value for p-values and adjusted p-values. Default is
#' 0.05.
#' @param lfc Minimum absolute log2-fold change to use as threshold for
#' differential expression. Default is 1.
#' @param sig Criteria to denote significance. Choices are \code{"adjP"}
#' (default) for adjusted p-value or \code{"P"} for p-value.
#' @param n_top Number of top hits to include in the heat map.
#' @param palette Viridis color palette option for plots. Default is
#' \code{"viridis"}. See
#' \code{\link[viridis:viridis]{viridis}}
#' for available options.
#' @param text_size Text size for axis text, labels etc.
#' @param save Logical. If \code{TRUE} saves a copy of the plot in the
#' directory provided in \code{file_path}.
#' @param file_path A string containing the directory path to save the file.
#' @param file_name File name to save the plot. Default is "HeatmapDE."
#' @param file_type File type to save the plot. Default is \code{"pdf".}
#' @param dpi Plot resolution. Default is \code{80.}
#' @param plot_height Height of the plot. Default is 7.
#' @param plot_width Width of the plot. Default is 7.
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
#'
#' ## Build a heatmap of differentially expressed proteins using the provided
#' ## example fit_df and norm_df data objects
#' heatmap_de(covid_fit_df, covid_norm_df)
#'
#' ## Create a heatmap with P-value of 0.05 and log fold change of 1 as
#' ## significance criteria.
#' heatmap_de(covid_fit_df, covid_norm_df, cutoff = 0.05, sig = "P")
#'
#' ## Visualize the top 30 differentially expressed proteins in the heatmap and
#' ## change the color palette
#' heatmap_de(covid_fit_df, covid_norm_df,
#'   cutoff = 0.05, sig = "P", n_top = 30,
#'   palette = "magma"
#' )
#'
#' @export
heatmap_de <- function(fit_df,
                       df,
                       adj_method = "BH",
                       cutoff = 0.05,
                       lfc = 1,
                       sig = "adjP",
                       n_top = 20,
                       palette = "viridis",
                       text_size = 10,
                       save = FALSE,
                       file_path = NULL,
                       file_name = "HeatmapDE",
                       file_type = "pdf",
                       dpi = 80,
                       plot_height = 7,
                       plot_width = 7) {
  # Binding the global variables to a local function
  logFC <- P.Value <- adj.P.Val <- intensity <- protein <- NULL

  # convert norm_df or imp_df object into a matrix
  norm_df <- as.matrix(df)

  # Extract the required data from the fit object
  exp_de <- limma::topTable(fit_df,
    coef = colnames(fit_df)[2],
    n = length(fit_df$df.total),
    adjust.method = adj_method
  )

  # Pick the sig. proteins based on lowest p-value and highest logFC

  if (sig == "P") {
    top_proteins <- rownames(subset(exp_de,
      abs(logFC) > lfc & P.Value < cutoff,
      drop = FALSE
    ))

    # Or default: based on adj.P value
  } else {
    top_proteins <- rownames(subset(exp_de,
      abs(logFC) > lfc & adj.P.Val < cutoff,
      drop = FALSE
    ))
  }

  # Extract the top n_top hits from the top hit list
  top_proteins <- top_proteins[1:n_top]

  # Check if there are sig. proteins before moving on to plotting
  if (identical(top_proteins, character(0))) {
    stop(message
    (paste0("No significant proteins found at ", sig, " < ", cutoff, ".")))
  } else {
    # Extract intensity values for top proteins based on logFC and p-val cutoff
    top_intensity <- subset(norm_df,
      rownames(norm_df) %in% top_proteins,
      drop = FALSE
    )
    # Convert to long format for plotting
    top_int_melted <- reshape2::melt(top_intensity)

    # add column names
    colnames(top_int_melted) <- c("protein", "sample", "intensity")

    # add new column based on cancer stage
    top_int_melted$stage <- sapply(
      strsplit(
        as.character(top_int_melted[, "sample"]), "_"
      ),
      getElement, 1
    )

    top_heatmap <- ggplot2::ggplot(
      top_int_melted,
      ggplot2::aes(
        x = sample,
        y = reorder(protein, intensity),
        fill = intensity
      )
    ) +
      ggplot2::geom_tile(
        colour = "white", size = 0.2,
        stat = "identity"
      ) +
      ggplot2::labs(
        x = "",
        y = ""
      ) +
      ggplot2::scale_y_discrete(
        expand = c(0, 0),
        labels = sapply(
          strsplit(
            as.character(top_int_melted[, "protein"]), ";"
          ),
          getElement, 1
        )
      ) +
      viridis::scale_fill_viridis(
        discrete = FALSE,
        direction = -1,
        option = palette,
        begin = 0,
        end = 1
      ) +
      promor_facet_theme() +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = grid::unit(text_size * 0.08, "cm"),
        legend.key.height = grid::unit(text_size * 0.02, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(
          size = text_size * 0.7,
          face = "bold"
        ),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = text_size * 0.7),
        plot.background = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(text_size * 0.05,
          units = "cm"
        ),
        panel.background = element_blank()
      ) +
      ggplot2::facet_grid(. ~ stage, scales = "free")

    # Set temporary file_path if not specified
    if (is.null(file_path)) {
      file_path <- tempdir()
    }

    if (save == TRUE) {
      ggplot2::ggsave(paste0(file_path, "/", file_name, ".", file_type),
        top_heatmap,
        dpi = dpi,
        height = plot_height,
        width = plot_width
      )
      return(top_heatmap)
    } else {
      return(top_heatmap)
    }
  }
}
