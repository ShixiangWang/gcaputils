# Generated by roxytest: Do not edit by hand!

# File R/circos.R: @testexamples

test_that("Function gcap.plotCircos() @ L44", {
  
  
  library(gcap)
  if (require("circlize") && require("scales")) {
    data("ascn")
    data <- ascn
  
    # Create fake data
    set.seed(1234)
    data$sample <- sample(LETTERS[1:10], nrow(data), replace = TRUE)
    rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
  
    gcap.plotCircos(rv)
  
    # Select genes to highlight in plot
    r1 <- rv$getGeneSummary(return_record = TRUE)
    gcap.plotCircos(rv, r1[amplicon_type == "circular"]$gene_id[1:10])
  
    gcap.plotCircos(rv, r1[, .(label = .N), by = .(band)][1:10])
    gcap.plotCircos(
      rv,
      r1[, .(label = .N, cluster = gsub("(.*):.*", "\\1", band)),
        by = .(band)
      ][1:10]
    )
  }
  
  expect_equal(5, 5)
})

