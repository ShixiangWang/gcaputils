#' Plot CN Mean vs. Frequency Dotplot
#'
#' @param fCNA a [fCNA](gcap::fCNA) object.
#' @param by the level of focal amplicons.
#' @param filter a filter expression based on `cn` (CN mean) and `N` (frequency)
#' to add text labels.
#' @param include amplicon type to include for plotting.
#' @param unique if `TRUE`, count samples instead of genes.
#' @param prob_cutoff,gap_cn thresholds for classify 'circular' and 'noncircular'.
#' @param ... other parameters passing to `ggrepel::geom_label_repel()`.
#'
#' @return a ggplot.
#' @export
#' @examples
#' \donttest{
#' library(gcap)
#' if (require("ggrepel") && require("cowplot")) {
#'   data("ascn")
#'   data <- ascn
#'
#'   # Create fake data
#'   set.seed(1234)
#'   data$sample <- sample(LETTERS[1:10], nrow(data), replace = TRUE)
#'   rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'
#'   p <- gcap.dotplot(rv, filter = cn > 60)
#'   p
#'   p2 <- gcap.dotplot(rv, filter = cn > 60 | N > 15, by = "band")
#'   p2
#'   p3 <- gcap.dotplot(rv, filter = cn > 60 | N > 50, by = "chr")
#'   p3
#'   p4 <- gcap.dotplot(rv, filter = cn > 60 | N > 3, by = "band", unique = TRUE)
#'   p4
#' }
#' }
#' @testexamples
#' expect_is(p, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "ggplot")
#' expect_is(p4, "ggplot")
gcap.dotplot <- function(fCNA,
                         filter,
                         by = c("gene_id", "band", "chr"),
                         unique = FALSE,
                         include = c("circular", "possibly_circular"),
                         prob_cutoff = 0.5,
                         gap_cn = 4L,
                         ...) {
  .check_install("ggrepel")
  .check_install("cowplot")
  stopifnot(inherits(fCNA, "fCNA") | is.data.frame(fCNA))
  if (is.data.frame(fCNA)) {
    data <- fCNA
    stopifnot(all(c("amplicon_type", "gene_id" %in% colnames(data))))
  } else {
    data <- fCNA$getGeneSummary(prob_cutoff, gap_cn, return_record = TRUE)
  }
  by <- match.arg(by)
  if (by == "chr") {
    data$chr <- gsub("(.*):(.*)", "\\1", data$band)
  }
  if (!unique) {
    genes_summary <- data[
      , .(
        cn = mean(total_cn[amplicon_type %in% include], na.rm = TRUE),
        N = sum(amplicon_type %in% include, na.rm = TRUE)
      ),
      by = by
    ][N > 0]
  } else {
    genes_summary <- data[
      , .(
        cn = mean(total_cn[amplicon_type %in% include], na.rm = TRUE),
        N = length(unique(sample[amplicon_type %in% include]))
      ),
      by = by
    ][N > 0]
  }
  genes_summary <- genes_summary[!is.na(genes_summary[[by]])]

  e <- substitute(filter)
  p <- ggplot2::ggplot(data = genes_summary, ggplot2::aes(x = N, y = cn)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2, col = "black") +
    ggrepel::geom_label_repel(
      ggplot2::aes_string(label = by),
      data = subset(genes_summary, eval(e)),
      ...
    ) +
    ggplot2::labs(
      x = if (unique) "Frequency" else "Gene frequency",
      y = "Copy number mean"
    ) +
    cowplot::theme_cowplot()
  p
}
