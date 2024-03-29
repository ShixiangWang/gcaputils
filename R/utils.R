.check_install = function(pkg, bioc = FALSE, ...) {
  install_func = if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) {
    message("installing required package {pkg}")
    install_func(pkg, ...)
  }
}

#' Cluster genomic regions by distance from region center
#' @param dt a genomic region data with first 3 columns for
#' chromosome, start and end.
#' @param distance a distance value for clustering,
#' regions with average distance lower than this value would be clustered
#' into one. Should not greater than 1e8.
#' @param no_annotation if `TRUE`, use builtin data for annotation. This is
#' useful when input a list of genes.
#' @param genome_build genome build version.
#' @param simplify if `FALSE`, return a `data.table` instead of a vector.
#' @return cluster list.
#' @export
#' @examples
#' library(data.table)
#' region = data.table(
#'   chr = c("chr1", "chr2", "chr1"),
#'   start = c(1032, 992, 10000), end = c(4242, 1e6, 400023)
#' )
#' clusterGPosition(region, distance = 1000)
#' clusterGPosition(region)
#' clusterGPosition(data.table(gene_id = c("KRAS", "EGFR", "MYC", "ERBB2", "GRB7")),
#'   no_annotation = TRUE, simplify = FALSE
#' )
#' clusterGPosition(data.table(gene_id = c("KRAS", "EGFR", "MYC", "ERBB2", "GRB7")),
#'   no_annotation = TRUE, simplify = TRUE
#' )
clusterGPosition = function(dt, distance = 1e7,
                             no_annotation = FALSE,
                             genome_build = c("hg38", "hg19"),
                             simplify = TRUE) {
  dt = if (data.table::is.data.table(dt)) data.table::copy(dt) else data.table::as.data.table(dt)
  genome_build = match.arg(genome_build)
  if (no_annotation) {
    stopifnot("gene_id" %in% colnames(dt))
    target_genes = readRDS(file.path(
      system.file("extdata", package = "gcap"),
      paste0(genome_build, "_target_genes.rds")
    ))
    if (!startsWith(dt$gene_id[1], "ENSG")) {
      message("detected you have transformed ENSEMBL ID, also transforming gene annotation data")
      opts = getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
      options(IDConverter.datapath = opts)
      target_genes$gene_id = IDConverter::convert_hm_genes(target_genes$gene_id, genome_build = genome_build)
      target_genes = target_genes[!is.na(target_genes$gene_id), ]
    }
    dt = merge(dt, target_genes, by = "gene_id", all.x = TRUE, sort = FALSE)
    data.table::setcolorder(dt, c("chrom", "start", "end"))
  }
  colnames(dt)[1:3] = c("chr", "start", "end")
  dt[, `:=`(x = as.integer(factor(chr)), y = (start + end) / 2)]
  dst = stats::as.dist(calc_dist(as.matrix(dt[, .(x, y)])))
  cls = stats::cutree(stats::hclust(dst, method = "average"), h = distance)
  if (simplify) {
    cls
  } else {
    dt$cluster = cls
    dt
  }
}

# refine_class = function(fcna, min_n = 2) {
#   v1 = unique(fcna$data[, list(sample, band, gene_class)])[
#     , list(class = ifelse(sum(gene_class == "circular") > min_n,
#                           "circular", "noncircular")), by = list(sample)]
#   v2 = fcna$sample_summary[, -"class"]
#   rv = merge(v2, v1, by = "sample", all.x = TRUE)
#   rv[, class := ifelse(is.na(class), "nofocal", class)]
#   return(rv)
# }


#' Set default factor level for fCNA class
#'
#' @param class a vector of fCNA class.
#' @param ref_level the reference level of factor.
#' @param default a default value for input `NA`s.
#'
#' @return a vector.
#' @export
#'
#' @examples
#' set_default_factor(c("nofocal", "noncircular", "circular", "nofocal", "possibly_circular"))
#' set_default_factor(c("nofocal", "noncircular", "circular", "nofocal", "possibly_circular", NA))
#' set_default_factor(c("nofocal", "noncircular", "circular", "nofocal", "possibly_circular", NA), 
#'                    default = "nofocal")
set_default_factor = function(class, ref_level = "nofocal", default = NA_character_) {
  class = data.table::fcase(class %in% c("circular", "possibly_circular"), "circular",
    class == "noncircular", "noncircular",
    class == "nofocal", "nofocal",
    default = default
  )
  factor(class, levels = c(ref_level, setdiff(c("nofocal", "noncircular", "circular"), ref_level)))
}

factor_to_chrs = function(x) {
  if (inherits(x, "factor")) as.character(x) else x
}

utils::globalVariables(
  c(
    "chr", "start", "end", ".", "x", "y",
    "cn", "N", "total_cn", "amplicon_type",
    "gene_id", "circular", "noncircular",
    "time", "status", "band"
  )
)
