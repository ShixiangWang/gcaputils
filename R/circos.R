#' Plot Circos
#'
#' @inheritParams gcap.dotplot
#' @param highlights gene/cytoband list to highlight.
#' @param total_n if not `NULL`, this is used for calculating percentage.
#' @param only if not `NA`, only shows a type of amplicon.
#' @param clust_distance a distance as cutoff for different clusters.
#' Default is 1e7, i.e. 10Mb. Note 100 Mb is set to genes on different
#' chromosomes, so please don't set value larger than that.
#' @param col length-2 colors for circular and noncircular.
#' @param genome_build genome version.
#' @param chrs chromosome names.
#' @param ideogram_height ideogram height at default.
#'
#' @return Nothing.
#' @export
#' @examples
#' \donttest{
#' library(gcap)
#' if (require("circlize") && require("scales")) {
#'   data("ascn")
#'   data = ascn
#'
#'   # Create fake data
#'   set.seed(1234)
#'   data$sample = sample(LETTERS[1:10], nrow(data), replace = TRUE)
#'   rv = gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'
#'   gcap.plotCircos(rv)
#'
#'   # Select genes to highlight in plot
#'   r1 = rv$getGeneSummary()
#'   gcap.plotCircos(rv, r1[circular == 1]$gene_id[1:10])
#'
#'   gcap.plotCircos(rv, rv$data[, .(label = .N), by = .(band)][1:10])
#'   gcap.plotCircos(
#'     rv,
#'     rv$data[, .(label = .N, cluster = gsub("(.*):.*", "\\1", band)),
#'       by = .(band)
#'     ][1:10]
#'   )
#' }
#' }
#' @testexamples
#' expect_equal(5, 5)
gcap.plotCircos = function(fCNA,
                            highlights = NULL,
                            total_n = NULL,
                            clust_distance = 1e7,
                            col = c("#FF000080", "#0000FF80"),
                            genome_build = c("hg38", "hg19"),
                            chrs = paste0("chr", 1:22),
                            ideogram_height = 1,
                            only = c(NA_character_, "circular", "noncircular"),
                            ...) {
  .check_install("circlize")
  .check_install("scales")
  genome_build = match.arg(genome_build)
  only = match.arg(only)
  target_genes = readRDS(file.path(
    system.file("extdata", package = "gcap"),
    paste0(genome_build, "_target_genes.rds")
  ))
  data = fCNA$getGeneSummary()
  if (!startsWith(data$gene_id[1], "ENSG")) {
    message("detected you have transformed ENSEMBL ID, also transforming gene annotation data")
    opts = getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
    options(IDConverter.datapath = opts)
    target_genes$gene_id = IDConverter::convert_hm_genes(target_genes$gene_id, genome_build = genome_build)
    target_genes = target_genes[!is.na(target_genes$gene_id), ]
  }
  data_bed = merge(data, target_genes, by = "gene_id", all.x = TRUE, sort = FALSE)
  data_bed[, amplicon_type := ifelse(circular > 0, "circular", "noncircular")]
  data_bed = data_bed[!is.na(data_bed$gene_id) & data_bed$chrom %in% chrs,
    c(
      "chrom", "start", "end",
      "gene_id", "amplicon_type"
    ),
    with = FALSE
  ]
  colnames(data_bed)[1] = "chr"
  if (nrow(data_bed) < 1) {
    message("no data to plot, please check!")
    return(invisible(NULL))
  }

  bed_list = split(data_bed, data_bed$amplicon_type)
  bed_list = lapply(bed_list, function(x) {
    if (!is.null(total_n)) {
      message("calculate percentage")
      x[, .(freq = .N / total_n), by = .(chr, start, end, gene_id)]
    } else {
      x[, .(freq = .N), by = .(chr, start, end, gene_id)]
    }
  })

  gap_after = c(rep(1, length(chrs) - 1), 12)
  circlize::circos.par(
    "start.degree" = 90,
    "track.height" = if (!is.na(only)) 0.5 else 0.25, # % of the circle radius
    "gap.after" = gap_after
  )

  track_height = 0.15
  circlize::circos.initializeWithIdeogram(
    species = genome_build,
    chromosome.index = chrs,
    track.height = circlize::convert_height(
      if (!is.null(highlights)) 0 else track_height,
      "mm"
    ), ideogram.height = circlize::convert_height(
      if (!is.null(highlights)) 0 else ideogram_height,
      "mm"
    ),
    plotType = if (!is.null(highlights)) NULL else c("ideogram", "axis", "labels")
  )
  on.exit(circlize::circos.clear())

  if (!is.null(highlights)) {
    # circlize::circos.genomicInitialize(
    #   circlize::read.cytoband(species = genome_build, chromosome.index = chrs)$df,
    #   plotType = "labels"
    # )
    if (is.data.frame(highlights)) {
      message("found input a data.frame for highlighting")
      if ("gene_id" %in% colnames(highlights)) {
        data = data.table::as.data.table(highlights)
        highlights = data$gene_id
      } else if ("band" %in% colnames(highlights)) {
        data = data.table::as.data.table(highlights)
        highlights = data$band
      } else {
        stop("when input a data.frame, 'gene_id' or 'band' column must exist")
      }
    } else {
      data = data.table::data.table()
    }

    if (!"band" %in% colnames(data)) {
      bed = unique(data_bed[
        data_bed$gene_id %in% highlights,
        c("chr", "start", "end", "gene_id")
      ])
      if (ncol(data) > 0) {
        bed = merge(bed, data, by = "gene_id", all.x = TRUE, sort = FALSE)
        data.table::setcolorder(bed, c("chr", "start", "end", "gene_id"))
      }

      if (nrow(bed) < 1) {
        message("no data for your selected genes, please check")
        return(invisible(NULL))
      }

      if (is.null(clust_distance) & ncol(data) == 0) {
        bed_col = as.numeric(factor(bed$gene_id))
      } else if ("cluster" %in% colnames(bed)) {
        message("found 'cluster' column, use it for color mapping")
        bed_col = as.numeric(factor(bed$cluster))
      } else {
        message("clustering highlight genes by 'hclust' average method with distance data from gene centers")
        message("distance cutoff: ", clust_distance)
        bed_col = clusterGPosition(bed, clust_distance)
        message("done")
      }

      if ("label" %in% colnames(bed)) {
        message("label detected, add it to gene id")
        bed$gene_id = paste0(bed$gene_id, " (", bed$label, ")")
      }

      data.table::setnames(bed, "gene_id", "id")
    } else {
      .check_install("sigminer")
      cytobands = sigminer::get_genome_annotation("cytobands", genome_build = genome_build)
      data.table::setDT(cytobands)
      colnames(cytobands)[1] = "chr"
      cytobands$stain = NULL
      cytobands[, band := paste(chr, band, sep = ":")]

      bed = merge(data, cytobands, by = "band", sort = FALSE)
      data.table::setcolorder(bed, c("chr", "start", "end", "band"))

      if ("cluster" %in% colnames(bed)) {
        message("found 'cluster' column, use it for color mapping")
        bed_col = as.numeric(factor(bed$cluster))
      } else {
        bed_col = as.numeric(factor(bed$band))
      }

      if ("label" %in% colnames(bed)) {
        message("label detected, add it to band id")
        bed$band = paste0(bed$band, " (", bed$label, ")")
      }

      data.table::setnames(bed, "band", "id")
    }

    ssize = max(nchar(highlights)) / 8
    circlize::circos.genomicLabels(bed,
      labels = bed$id, side = "outside",
      cex = 0.5 / (ssize),
      connection_height = circlize::mm_h(5 / ssize),
      col = bed_col,
      line_col = bed_col
    )
    circlize::circos.genomicIdeogram(
      species = genome_build,
      track.height = 0.02 * ideogram_height
    )
    circlize::circos.track(
      track.index = circlize::get.current.track.index(),
      panel.fun = function(x, y) {
        circlize::circos.genomicAxis(h = "top", direction = "outside")
      }
    )
  }

  if (is.na(only)) {
    draw_freq_track(bed_list$circular, col[1])
    draw_track_y_axis(chrs)
    draw_freq_track(bed_list$noncircular, col[2])
    draw_track_y_axis(chrs)
  } else if (only == "circular") {
    draw_freq_track(bed_list$circular, col[1])
    draw_track_y_axis(chrs)
  } else {
    draw_freq_track(bed_list$noncircular, col[1])
    draw_track_y_axis(chrs)
  }

  # draw_bi_track(bed_cn, col)
  # draw_track_y_axis(chrs, at = scales::breaks_pretty(5)(c(0, cnrange[2])))
}

draw_track_y_axis = function(chrs, at = NULL) {
  circlize::circos.yaxis(
    at = at,
    labels.cex = 0.3,
    tick.length = circlize::convert_x(
      0.3,
      "mm", circlize::get.cell.meta.data("sector.index"),
      circlize::get.cell.meta.data("track.index")
    ),
    side = "right", sector.index = chrs[length(chrs)]
  )
}

draw_freq_track = function(data, col) {
  data$gene_id = NULL # Only keep freq as value
  circlize::circos.genomicTrack(
    data,
    ylim = c(0, max(data$freq, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      # circlize::circos.genomicRect(region, value,
      #   ytop.column = 1,
      #   ybottom = 0,
      #   col = col,
      #   border = NA
      # )
      circlize::circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = col, border = NA)
    }
  )
}

draw_bi_track = function(data, col) {
  circlize::circos.genomicTrack(
    data,
    # track.height = if (!is.null(highlights)) 0.2 else 0.3, # ylim = c(0, 0.1),
    panel.fun = function(region,
                         value, ...) {
      circlize::circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = col, border = NA)

      # circlize::circos.genomicRect(
      #   region, value,
      #   ytop = value$circular,
      #   ybottom = value$noncircular,
      #   col = ifelse(value$circular > 0, col[1], col[2]), border = NA
      # )
      circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#000040")
    }
  )
}
