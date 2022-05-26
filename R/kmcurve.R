#' Draw K-M Curve for fCNA Survival Comparison by Sample or Cytoband or Gene
#'
#' @param fCNA a [fCNA](gcap::fCNA) object.
#' @param surv_data survival data, eithor a 3-column `data.frame` to store
#' sample, time and status, or a length-2 string to specify the colnames
#' representing time and status in `fCNA$sample_summary`.
#' - sample must be identical to sample ID in `fCNA`.
#' - time must be numeric.
#' - status must be 0 or 1.
#' @param mat a gene/cytoband-by-sample matrix like `data.frame`.
#' @param ID a list of gene or cytoband IDs.
#' @param merge_circular if `TRUE`, merge 'circular' and 'possibly_circular'
#' as one class.
#' @param focus focal amplication type you focus on.
#' Can be 'fCNA' or 'circular'. If 'fCNA' selected,
#' noncircular and circular genes/cytobands are included to classify samples.
#' @param palette plot color palette.
#' @param class_col column name in `sample_summary` field for classification.
#' If you set to other column (you want to run survival analysis with custom column),
#' parameters like `merge_circular`, `ID`, `focus`
#' etc. will be omitted.
#' @param ending_time survival analysis ending time. If a numeric ending
#' is typed, all survival data longer than the ending time will be rewritten.
#' @param ... other parameters passing to `survminer::ggsurvplot`.
#'
#' @return a plot.
#' @export
#' @seealso [gcap.plotProfile] for plot landscape of fCNA, [gcap::fCNA] for building object.
#' @examples
#' \donttest{
#' library(gcap)
#' if (require("survminer") && require("IDConverter")) {
#'   data("ascn")
#'   data <- ascn
#'
#'   # Create fake data
#'   set.seed(1234)
#'   data$sample <- sample(LETTERS[1:10], nrow(data), replace = TRUE)
#'   rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'   rv$convertGeneID()
#'
#'   surv_data <- data.frame(
#'     sample = rv$sample_summary$sample,
#'     time = 3000 * abs(rnorm(nrow(rv$sample_summary))),
#'     status = sample(c(0, 1), nrow(rv$sample_summary), replace = TRUE)
#'   )
#'   p <- gcap.plotKMcurve(rv, surv_data)
#'   p
#'
#'   p2 <- gcap.plotKMcurve(rv, surv_data,
#'     ID = "MYC",
#'     mat = rv$getGeneSummary(return_mat = TRUE)
#'   )
#'
#'   p2
#' }
#' }
#' @testexamples
#' expect_s3_class(p, "ggsurvplot")
#' expect_s3_class(p2, "ggsurvplot")
gcap.plotKMcurve <- function(fCNA,
                             surv_data,
                             merge_circular = TRUE,
                             mat = NULL,
                             ID = NULL,
                             focus = c("fCNA", "circular"),
                             palette = c("grey", "#0066CC", "#CC0033"),
                             class_col = "class",
                             ending_time = NULL,
                             ...) {
  stopifnot(inherits(fCNA, "fCNA") | is.data.frame(fCNA))
  .check_install("survminer")
  focus <- match.arg(focus)
  if (is.data.frame(fCNA)) {
    if (ncol(fCNA) != 4) {
      stop("when input fCNA is a data.frame, columns should be sample, class, time and status")
    }
  }

  if (!is.data.frame(fCNA)) {
    if (is.character(surv_data)) {
    surv_data <- fCNA$sample_summary[, c("sample", surv_data), with = FALSE]
  }
  colnames(surv_data)[2:3] <- c("time", "status")
  }

  if (is.null(ID)) {
    if (is.data.frame(fCNA)) {
      data = data.table::as.data.table(fCNA)
      colnames(data) = c("sample", "class", "time", "status")
    } else {
      data <- fCNA$sample_summary[, c("sample", class_col), with = FALSE]
      colnames(data)[2] <- "class"
    }
    if (merge_circular & class_col == "class") {
      data[, class := set_default_factor(class)]
    } else if (class_col == "class") {
      data[, class := factor(class, c("nofocal", "noncircular", "possibly_circular", "circular"))]
    }
  } else {
    if (is.null(mat)) {
      stop("When you want to specify the genes/cytobands, please input gene/cytoband-by-sample matrix to 'mat'")
    }
    # Extract class based on gene/cytoband

    # labels for AMP samples
    if (focus == "fCNA") {
      types <- c("noncircular", "possibly_circular", "circular")
      labels <- c("fCNA-", "fCNA+")
    } else if (merge_circular) {
      types <- c("possibly_circular", "circular")
      labels <- c("circular-", "circular+")
    } else {
      types <- "circular"
      labels <- c("circular-", "circular+")
    }
    data <- mat[ID, , drop = FALSE]
    data <- plyr::ldply(data,
      function(x) if (any(x %in% types)) labels[2] else labels[1],
      .id = "sample"
    )
    colnames(data)[2] <- "class"
    if (length(table(data$class)) <= 1) {
      stop("cannot genrate two groups for comparison")
    }

    data$class <- factor(data$class, levels = labels)
  }

  if (!is.data.frame(fCNA)) {
    data <- merge(data, data.table::as.data.table(surv_data),
    by = "sample",
    all.x = TRUE
  )
  }

  if (!is.null(ending_time)) {
    data[, status := ifelse(time >= ending_time, 0, status)]
    data[, time := ifelse(time >= ending_time, ending_time, time)]
  }

  fit <- survminer::surv_fit(survival::Surv(time, status) ~ class, data = data)
  print(survminer::surv_pvalue(fit = fit))

  cls_lvls <- sub("class=", "", names(fit$strata))
  if (identical(palette, c("grey", "#0066CC", "#CC0033"))) {
    if (length(cls_lvls) == 2) palette <- palette[2:3]
  }

  p <- survminer::ggsurvplot(fit,
    pval = TRUE, data = data,
    palette = palette,
    risk.table = TRUE,
    legend.labs = cls_lvls,
    ...
  )
  p
}
