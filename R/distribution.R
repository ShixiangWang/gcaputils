#' Draw Gene Amplicon Distribution By User Specified Class
#'
#' @param fCNA a `fCNA` object or a `data.frame` with at least
#' 3 columns "sample", "class" and "by".
#' @param x a column name in `fCNA$sample_summary`.
#' @param x_size font size of x axis text.
#' @param fill if `TRUE`, show the percentage instead of count.
#' @param set_order set order for `by`.
#' @param set_label set labels to add sample count for `by`.
#' @param palette color palette.
#' @param plot if `FALSE`, return data instead of plot.
#' @param by_gene sow distribution for genes.
#' @param genelist only listed gene will be shown.
#' @param bar_width set the width of bar, when it is `1`, the border can be removed.
#' @param ... other parameters passing to `ggplot2::geom_bar`.
#'
#' @return a ggplot object.
#' @export
#' @examples
#' set.seed(1234)
#' data <- data.frame(
#'   sample = sample(LETTERS[1:10], 100, replace = TRUE),
#'   class = sample(c("nofocal", "noncircular", "circular"), 100, replace = TRUE),
#'   by = sample(1:4, 100, replace = TRUE)
#' )
#' p <- gcap.plotDistribution(data)
#' p
#' @testexamples
#' expect_is(p, "ggplot")
gcap.plotDistribution <- function(fCNA,
                                  x = NULL,
                                  x_size = 8,
                                  set_order = TRUE,
                                  set_label = TRUE,
                                  fill = TRUE,
                                  palette = c("#CCCCCC", "#0066CC", "#FFCCCC", "#CC0033"),
                                  plot = TRUE,
                                  by_gene = FALSE,
                                  genelist = NULL,
                                  bar_width = 0.9,
                                  ...) {
  stopifnot(inherits(fCNA, "fCNA") | inherits(fCNA, "data.frame"))
  .check_install("ggplot2")

  if (!is.data.frame(fCNA)) {
    if (!by_gene) {
      data <- fCNA$sample_summary[, c("sample", "class", x), with = FALSE]
      colnames(data)[3] <- "by"
    } else {
      data <- fCNA$data[, c("sample", "gene_class", "gene_id"), with = FALSE]
      if (!is.null(genelist)) {
        data = data[gene_id %in% genelist]
      }
      colnames(data) = c("sample", "class", "by")
    }
  } else {
    data <- data.table::as.data.table(fCNA)[, c("sample", "class", "by"), with = FALSE]
  }

  if (by_gene) {
    class_lvls <- c("noncircular", "circular")
    data[, class := factor(class, class_lvls)]
  } else {
    data[, class := set_default_factor(class)]
    class_lvls <- levels(data$class)
  }
  dt_n <- data[, .N, by = list(by, class)]

  if (set_label) {
    dt_s <- dt_n[, list(N = sum(N, na.rm = TRUE)), by = list(by)]
    dt_s$label <- paste0(dt_s$by, " (N=", dt_s$N, ")")
    labels <- dt_s$label
    names(labels) <- dt_s$by
  } else {
    labels <- NULL
  }

  if (fill) {
    dt_n <- dt_n[, list(class, N = N / sum(N)), by = list(by)]
  }

  if (set_order) {
    by_order <- dt_n[, list(N = ifelse("circular" %in% class, N[class == "circular"], if (fill) 0 else 0L)),
                     by = list(by)]
    by_order <- by_order[order(-N)]$by
    dt_n[, by := factor(by, levels = by_order)]
  }


  if (length(class_lvls) == 3 && identical(palette, c("#CCCCCC", "#0066CC", "#FFCCCC", "#CC0033"))) {
    palette <- palette[-3]
  }

  if (plot) {
    p <- ggplot2::ggplot(dt_n, ggplot2::aes(by, N, fill = class)) +
      ggplot2::geom_bar(stat = "identity", width = bar_width, ...) +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion()) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::scale_x_discrete(position = "top", labels = labels) + # guide = ggplot2::guide_axis(n.dodge = 2)
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          size = x_size,
          vjust = 0,
          hjust = 0,
          color = "black",
          face = "bold"
        )
      ) +
      ggplot2::labs(x = NULL, y = NULL)

    p
  } else {
    dt_n
  }
}
