gcap.plotDistHeatmap <- function(fCNA,
                                 group = NULL,
                                 highlights = NULL,
                                 total_n = NULL,
                                 draw = TRUE,
                                 w = 1e4,
                                 title = "numeric matrix",
                                 col = c("#FF000080", "#0000FF80"),
                                 genome_build = c("hg38", "hg19"),
                                 chrs = paste0("chr", 1:22),
                                 ...) {
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

    # bed1 = generateRandomBed(nr = 1000, nc = 10) # generateRandomBed() is from circlize package
    # # convert to a GRanes object
    # gr1 = GRanges(seqnames = bed1[, 1], ranges = IRanges(bed1[, 2], bed1[, 3]))
    # num_mat = average_in_window(chr_window, gr1, bed1[, -(1:3)])

    dt <- fCNA$data
    
    ref_file <- system.file(
        "extdata", paste0(genome_build, "_target_genes.rds"),
        package = "gcap", mustWork = TRUE
    )
    ref <- readRDS(ref_file)
    colnames(ref)[1:3] <- c("chr", "start", "end")

    if (is.null(group)) {
        group = "sample"
        dt_gene = dt[, .N, by = list(sample, gene_id, gene_class)]
        if (!is.null(total_n)) {
            dt_gene[, value := N / total_n]
        } else {
            dt_gene[, value := N / length(unique(dt_gene$sample))]
        }
        dt_gene = dt_gene[gene_class == "circular", list(sample, gene_id, value)]
        subgroup = unique(dt_gene$sample)
    } else {
        subgroup <- rep(c("A", "B"), each = 5)  
    }
    colnames(dt_gene) = c("group", "gene_id", "value")

    dt_gene = merge(ref, dt_gene, by = "gene_id")
    
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
    bed3 <- generateRandomBed(nr = 40, nc = 0)
    gr3 <- GRanges(seqnames = bed3[, 1], ranges = IRanges(bed3[, 2], bed3[, 2]))
    gr3$gene <- paste0("gene_", 1:length(gr3))

    mtch <- as.matrix(findOverlaps(chr_window, gr3))
    at <- mtch[, 1]
    labels <- mcols(gr3)[mtch[, 2], 1]

    chr <- as.vector(seqnames(chr_window))
    chr_level <- chrs
    chr <- factor(chr, levels = chrs)

    ht_opt$TITLE_PADDING <- unit(c(4, 4), "points")
    ht_list <- Heatmap(num_mat,
        name = title, col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
        row_split = chr, cluster_rows = FALSE, show_column_dend = FALSE,
        column_split = subgroup, cluster_column_slices = FALSE,
        top_annotation = HeatmapAnnotation(subgroup = subgroup, annotation_name_side = "left"),
        row_title_rot = 0, row_title_gp = gpar(fontsize = 10), border = TRUE,
        row_gap = unit(0, "points")
    )

    if (!is.null(highlights)) ht_list <- ht_list + rowAnnotation(label = anno_mark(at = at, labels = labels))

    if (draw) {
        draw(ht_list, merge_legend = TRUE)
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