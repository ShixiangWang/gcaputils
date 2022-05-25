#' Draw Oncoprint for fCNA Profile Visualization
#'
#' See [oncoprint](https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html)
#' and [heatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#size-of-the-heatmap)
#' to modify the result plot in depth.
#'
#' @param data a `data.frame` to draw oncoprint, rows represent
#' genes/cytobands and columns represent samples.
#' @param genes a list of genes. Must be subsets of `fCNA$data$gene_id`.
#' @param top_n top N genes to show.
#' @param samples a vector to subset samples.
#' @param show_column_names see `?ComplexHeatmap::oncoPrint`.
#' @param remove_empty_columns see `?ComplexHeatmap::oncoPrint`.
#' @param remove_empty_rows see `?ComplexHeatmap::oncoPrint`.
#' @param heatmap_legend_param see `?ComplexHeatmap::oncoPrint`.
#' @param col see `?ComplexHeatmap::oncoPrint`.
#' @param ... Other parameters passing to `ComplexHeatmap::oncoPrint`.
#'
#' @return a `Heatmap` object, you can go further modify it.
#' @export
#' @seealso [gcap::fCNA] for building object.
#'
#' @examples
#' \donttest{
#' library(gcap)
#' if (require("ComplexHeatmap") && require("IDConverter")) {
#'   data("ascn")
#'   data <- ascn
#'
#'   # Create fake data
#'   set.seed(1234)
#'   data$sample <- sample(LETTERS[1:10], nrow(data), replace = TRUE)
#'   rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'
#'   data2 <- rv$getGeneSummary(return_mat = TRUE)
#'   gcap.plotProfile(data2)
#'
#'   rv$convertGeneID()
#'   data2 <- rv$getGeneSummary(return_mat = TRUE)
#'   ht <- gcap.plotProfile(data2,
#'     samples = c("B", "A", "C", "D", "F"),
#'     genes = rownames(data2)[1:10],
#'     top_annotation = ComplexHeatmap::HeatmapAnnotation(
#'       cbar = ComplexHeatmap::anno_oncoprint_barplot(),
#'       foo1 = c("g1", "g1", "g2", "g2", "g3"),
#'       bar1 = ComplexHeatmap::anno_points(1:5)
#'     )
#'   )
#'   ht
#'   ComplexHeatmap::draw(ht, merge_legends = TRUE)
#'
#'   data2 <- rv$getCytobandSummary(return_mat = TRUE)
#'   gcap.plotProfile(data2)
#' }
#' }
#' @testexamples
#' expect_type(ht, "S4")
gcap.plotProfile <- function(data,
                             genes = NULL,
                             samples = NULL,
                             top_n = NULL,
                             show_column_names = TRUE,
                             remove_empty_columns = FALSE,
                             remove_empty_rows = TRUE,
                             heatmap_legend_param = list(title = "fCNA"),
                             col = c(
                               "noncircular" = "#0066CC",
                               "circular" = "#CC0033"
                             ),
                             ...) {
  stopifnot(inherits(data, "data.frame") | inherits(data, "matrix"))
  .check_install("ComplexHeatmap")

  if (!is.null(genes)) {
    data <- data[genes, , drop = FALSE]
    if (nrow(data) == 0) {
      warning("No gene left to plot", immediate. = TRUE)
      return(NULL)
    }
  }

  if (!is.null(samples)) {
    data <- data[, samples, drop = FALSE]
  }

  if (!is.null(top_n)) {
    top_n <- min(top_n, nrow(data))
    data <- data[seq_len(top_n), ]
  }

  alter_fun <- list(
    background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
    circular = ComplexHeatmap::alter_graphic("rect", fill = col["circular"]),
    noncircular = ComplexHeatmap::alter_graphic("rect", fill = col["noncircular"])
  )

  data2 <- as.data.frame(lapply(data, factor_to_chrs))
  data2[is.na(data2)] <- ""
  rownames(data2) <- rownames(data)
  data <- data2

  ht <- ComplexHeatmap::oncoPrint(data,
    alter_fun = alter_fun, col = col,
    heatmap_legend_param = heatmap_legend_param,
    remove_empty_columns = remove_empty_columns,
    remove_empty_rows = remove_empty_rows,
    show_column_names = show_column_names,
    ...
  )
  return(ht)
}
