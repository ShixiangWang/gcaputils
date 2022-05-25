#' Draw Gene Amplicon Distribution By User Specified Class
#'
#' @param fCNA a `fCNA` object or a `data.frame` with at least
#' 3 columns "sample", "class" and "by".
#' @param x a column name in `fCNA$sample_summary`.
#' @param x_size font size of x axis text.
#' @param fill if `TRUE`, show the percentage instead of count.
#' @param palette color palette.
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
                                  merge_circular = TRUE,
                                  fill = TRUE,
                                  palette = c("#CCCCCC", "#0066CC", "#FFCCCC", "#CC0033"),
                                  ...) {
  stopifnot(inherits(fCNA, "fCNA") | inherits(fCNA, "data.frame"))
  .check_install("ggplot2")

  if (!is.data.frame(fCNA)) {
    data <- fCNA$sample_summary[, c("sample", "class", x), with = FALSE]
    colnames(data)[3] <- "by"
  } else {
    data <- data.table::as.data.table(fCNA)[, c("sample", "class", "by"), with = FALSE]
  }

  if (merge_circular) {
    data[, class := set_default_factor(class)]
    class_lvls <- levels(data$class)
  } else {
    class_lvls <- c("nofocal", "noncircular", "possibly_circular", "circular")
    data[, class := factor(class, class_lvls)]
  }

  dt_n <- data[, .N, by = list(by, class)]

  if (fill) {
    dt_n <- dt_n[, list(class, N = N / sum(N)), by = list(by)]
  }

  if (length(class_lvls) == 3 && identical(palette, c("#CCCCCC", "#0066CC", "#FFCCCC", "#CC0033"))) {
    palette <- palette[-3]
  }

  p <- ggplot2::ggplot(dt_n, ggplot2::aes(by, N, fill = class)) +
    ggplot2::geom_bar(stat = "identity", ...) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion()) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_discrete(position = "top") + # guide = ggplot2::guide_axis(n.dodge = 2)
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
}
