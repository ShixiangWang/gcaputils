#' Plot Genome-level Heatmap
#'
#' If `highlight` is set, make sure use `pdf()` or other device for
#' rendering plot.
#'
#' @inheritParams gcap.dotplot
#' @param highlights gene list or region (in `data.frame`) to highlight.
#' @param total_n if not `NULL`, this is used for calculating percentage.
#' @param genome_build genome version.
#' @param chrs chromosome names.
#' @param type data type to plot.
#' @param group a group column for comparison, default is "sample".
#' @param draw if `TRUE`, draw plot.
#' @param bycol if `FALSE`, draw in row manner.
#' @param w window size.
#' @param title title for heatmap legend.
#' @param top_col color setting.
#' @param max_val maximum value in heatmap.
#' @references https://jokergoo.github.io/ComplexHeatmap-reference/book/genome-level-heatmap.html
#'
#' @return Nothing.
#' @export
#' @examples
#' \donttest{
#' library(gcap)
#' if (require("circlize") && require("scales")) {
#'   data("ascn")
#'   data <- ascn
#'
#'   # Create fake data
#'   set.seed(1234)
#'   data$sample <- sample(LETTERS[1:10], nrow(data), replace = TRUE)
#'   rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'
#'   gcap.plotGenomeHeatmap(rv, highlights = c("TP53", "EGFR", "ERBB2"))
#' }
#' }
gcap.plotGenomeHeatmap <- function(fCNA,
                                   type = c(
                                     "prob", "circ_freq",
                                     "noncirc_freq", "circ_cn",
                                     "noncirc_cn"
                                   ),
                                   group = NULL,
                                   highlights = NULL,
                                   total_n = NULL,
                                   draw = TRUE,
                                   bycol = TRUE,
                                   sort = FALSE,
                                   w = 1e6,
                                   title = "auto",
                                   top_col = "red",
                                   max_val = NULL,
                                   genome_build = c("hg38", "hg19"),
                                   chrs = paste0("chr", 1:22)) {
  .check_install("GenomicRanges", bioc = TRUE)
  .check_install("circlize", bioc = TRUE)
  .check_install("ComplexHeatmap", bioc = TRUE)
  .check_install("EnrichedHeatmap", bioc = TRUE)

  type <- match.arg(type)
  if (title == "auto") {
    title <- type
  }

  library(circlize)
  library(GenomicRanges)
  library(EnrichedHeatmap)
  library(ComplexHeatmap)

  genome_build <- match.arg(genome_build)

  chr_df <- read.chromInfo(species = genome_build)$df
  chr_df <- chr_df[chr_df$chr %in% chrs, ]
  chr_gr <- GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
  chr_gr

  chr_window <- makeWindows(chr_gr, w = w)
  chr_window

  dt_s <- fCNA$sample_summary
  dt <- fCNA$data
  if (!is.null(group)) stopifnot(group %in% colnames(fCNA$sample_summary))

  ref_file <- system.file(
    "extdata", paste0(genome_build, "_target_genes.rds"),
    package = "gcap", mustWork = TRUE
  )
  ref <- readRDS(ref_file)
  colnames(ref)[1:3] <- c("chr", "start", "end")
  if (!startsWith(dt$gene_id[1], "ENSG")) {
    message("detected you have transformed ENSEMBL ID, also transforming gene annotation data")
    opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
    options(IDConverter.datapath = opts)
    ref$gene_id <- IDConverter::convert_hm_genes(ref$gene_id, genome_build = genome_build)
    ref <- ref[!is.na(ref$gene_id), ]
  }

  if (is.null(group)) group <- "sample"
  if (!is.null(group) && group == "whole") {
    dt_s$whole <- "whole"
  }
  dt_s <- dt_s[, unique(c("sample", group)), with = FALSE]
  dt_n <- dt_s[, .N, by = list(group = dt_s[[group]])]

  dt2 <- merge(dt, dt_s, by = "sample")
  data.table::setnames(dt2, group, "group")

  # Note, only amplicon data used for a group
  # i.e., non-amplicon data are not included in summary
  if (type == "prob") {
    dt_gene <- dt2[, list(value = mean(prob, na.rm = TRUE)), by = list(group, gene_id)]
  } else if (type == "circ_cn") {
    dt_gene <- dt2[gene_class == "circular",
      list(value = mean(total_cn, na.rm = TRUE)),
      by = list(group, gene_id)
    ]
  } else if (type == "noncirc_cn") {
    dt_gene <- dt2[gene_class == "noncircular",
      list(value = mean(total_cn, na.rm = TRUE)),
      by = list(group, gene_id)
    ]
  } else {
    if (type == "circ_freq") {
      dt_gene <- dt2[gene_class == "circular",
        list(value = .N),
        by = list(group, gene_id)
      ]
    } else {
      dt_gene <- dt2[gene_class == "noncircular",
        list(value = .N),
        by = list(group, gene_id)
      ]
    }
    if (is.null(total_n)) {
      dt_gene <- merge(dt_gene, dt_n, by = "group")
      dt_gene[, value := value / N]
      dt_gene$N <- NULL
    } else {
      dt_gene[, value := value / total_n]
    }
  }
  subgroup <- unique(dt_gene$group)
  colnames(dt_gene) <- c("group", "gene_id", "value")

  dt_gene <- merge(ref, dt_gene, by = "gene_id")
  if (sort) {
    dt_sort <- dt2[, list(N = length(unique(band[gene_class == "circular"]))), by = list(group)][order(N)]
    dt_gene[, group := factor(group, dt_sort$group)]
    subgroup = dt_sort$group
  }
  
  bed_list <- lapply(split(dt_gene, dt_gene$group), function(dt) {
    dt[, list(chr, start, end, value)]
  })

  num_mat <- NULL
  for (i in seq_along(bed_list)) {
    bed <- as.data.frame(bed_list[[i]])
    gr_cnv <- GRanges(seqnames = bed[, 1], ranges = IRanges(bed[, 2], bed[, 3]))
    num_mat <- cbind(num_mat, average_in_window(chr_window, gr_cnv, bed[, 4]))
  }

  # visualize is a list of gene symbols that we want to mark in the plot.
  # gr3 contains genomic positions for the genes as well as their symbols.
  if (!is.null(highlights)) {
    if (is.data.frame(highlights)) {
      # A data.frame as bed with 4 columns
      stopifnot(ncol(highlights) == 4)
      bed3 <- as.data.frame(highlights)
    } else {
      if (is.character(highlights)) {
        if (startsWith(highlights[1], "ENSG") || !startsWith(ref$gene_id[1], "ENSG")) {
          bed3 <- as.data.frame(ref[gene_id %in% highlights])
        } else {
          message("detected you are using gene symbol, transforming gene annotation data")
          opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
          options(IDConverter.datapath = opts)
          ref$gene_id <- IDConverter::convert_hm_genes(ref$gene_id, genome_build = genome_build)
          bed3 <- as.data.frame(ref[gene_id %in% highlights])
        }
      } else {
        s
        stop("unsupported input for 'highlights'")
      }
    }
    gr3 <- GRanges(seqnames = bed3[, 1], ranges = IRanges(bed3[, 2], bed3[, 2]))
    gr3$gene <- bed3[[4]]

    mtch <- as.matrix(findOverlaps(chr_window, gr3))
    at <- mtch[, 1]
    labels <- mcols(gr3)[mtch[, 2], 1]
  }

  chr <- as.vector(seqnames(chr_window))
  chr <- factor(chr, levels = chrs)

  if (type == "prob") {
    if (is.null(max_val)) {
      colors <- colorRamp2(
        c(0, 0.5, round(max(dt_gene$value, na.rm = TRUE) + 0.005, 2)),
        c("blue", "white", top_col)
      )
    } else {
      colors <- colorRamp2(c(0, 0.5, max_val), c("blue", "white", top_col))
    }
  } else {
    if (is.null(max_val)) {
      colors <- colorRamp2(
        c(0, round(max(dt_gene$value, na.rm = TRUE) + 0.005, 2)),
        c("white", top_col)
      )
    } else {
      colors <- colorRamp2(c(0, max_val), c("white", top_col))
    }
  }

  # ht_opt$TITLE_PADDING <- unit(c(4, 4), "points")
  if (bycol) {
    ht_list <- Heatmap(num_mat,
      name = title, col = colors,
      # heatmap_legend_param = list(color_bar = "discrete", at = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
      row_split = chr, cluster_rows = FALSE,
      column_split = subgroup,
      cluster_column_slices = FALSE, cluster_columns = FALSE,
      top_annotation = HeatmapAnnotation(
        group = subgroup,
        show_legend = FALSE,
        show_annotation_name = FALSE
      ),
      column_title_rot = 90, column_title_gp = gpar(fontsize = 8),
      row_title_rot = 0, row_title_gp = gpar(fontsize = 8), border = TRUE,
      row_gap = unit(0, "points")
    )

    if (!is.null(highlights)) {
      ht_list <- ht_list + rowAnnotation(
        label = anno_mark(
          at = at, labels = labels,
          labels_gp = gpar(fontsize = 8)
        )
      )
    }
  } else {
    ht_list <- Heatmap(
      t(num_mat),
      name = title,
      col = colors,
      column_split = chr, cluster_columns = FALSE, show_row_dend = FALSE,
      row_split = subgroup, cluster_row_slices = FALSE,
      left_annotation = rowAnnotation(
        group = subgroup, show_annotation_name = FALSE,
        annotation_legend_param = list(
          group = list(direction = "horizontal", title_position = "lefttop", nrow = 1)
        )
      ),
      column_title_gp = gpar(fontsize = 10), border = TRUE,
      column_gap = unit(0, "points"),
      column_title = ifelse(seq_along(chrs) %% 2 == 0,
        paste0("\n", chrs), paste0(chrs, "\n")
      ),
      heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop")
    )

    if (!is.null(highlights)) {
      ht_list <-
        HeatmapAnnotation(label = anno_mark(at = at, labels = labels, side = "top")) %v%
        ht_list
    }
  }

  if (draw) {
    if (bycol) {
      draw(ht_list, merge_legend = TRUE)
    } else {
      draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE)
    }
  } else {
    ht_list
  }
}

# Copy from https://jokergoo.github.io/ComplexHeatmap-reference/book/genome-level-heatmap.html
average_in_window <- function(window, gr, v, method = "weighted", empty_v = NA) {
  if (missing(v)) v <- rep(1, length(gr))
  if (is.null(v)) v <- rep(1, length(gr))
  if (is.atomic(v) && is.vector(v)) v <- cbind(v)

  v <- as.matrix(v)
  if (is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }

  if (length(empty_v) == 1) {
    empty_v <- rep(empty_v, ncol(v))
  }

  u <- matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))

  mtch <- as.matrix(findOverlaps(window, gr))
  intersect <- pintersect(window[mtch[, 1]], gr[mtch[, 2]])
  w <- width(intersect)
  v <- v[mtch[, 2], , drop = FALSE]
  n <- nrow(v)

  ind_list <- split(seq_len(n), mtch[, 1])
  window_index <- as.numeric(names(ind_list))
  window_w <- width(window)

  if (is.character(v)) {
    for (i in seq_along(ind_list)) {
      ind <- ind_list[[i]]
      if (is.function(method)) {
        u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
      } else {
        tb <- tapply(w[ind], v[ind], sum)
        u[window_index[i], ] <- names(tb[which.max(tb)])
      }
    }
  } else {
    if (method == "w0") {
      gr2 <- reduce(gr, min.gapwidth = 0)
      mtch2 <- as.matrix(findOverlaps(window, gr2))
      intersect2 <- pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])

      width_intersect <- tapply(width(intersect2), mtch2[, 1], sum)
      ind <- unique(mtch2[, 1])
      width_setdiff <- width(window[ind]) - width_intersect

      w2 <- width(window[ind])

      for (i in seq_along(ind_list)) {
        ind <- ind_list[[i]]
        x <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
        u[window_index[i], ] <- (x * width_intersect[i] + empty_v * width_setdiff[i]) / w2[i]
      }
    } else if (method == "absolute") {
      for (i in seq_along(ind_list)) {
        u[window_index[i], ] <- colMeans(v[ind_list[[i]], , drop = FALSE])
      }
    } else if (method == "weighted") {
      for (i in seq_along(ind_list)) {
        ind <- ind_list[[i]]
        u[window_index[i], ] <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
      }
    } else {
      if (is.function(method)) {
        for (i in seq_along(ind_list)) {
          ind <- ind_list[[i]]
          u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }

  return(u)
}