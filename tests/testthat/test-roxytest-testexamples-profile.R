# Generated by roxytest: Do not edit by hand!

context("File R/profile.R: @testexamples")

test_that("Function gcap.plotProfile() @ L58", {
  
  
  library(gcap)
  if (require("ComplexHeatmap") && require("IDConverter")) {
    data("ascn")
    data = ascn
  
    # Create fake data
    set.seed(1234)
    data$sample = sample(LETTERS[1:10], nrow(data), replace = TRUE)
    rv = gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
  
    data2 = rv$getGeneSummary(return_mat = TRUE)
    gcap.plotProfile(data2)
  
    rv$convertGeneID()
    data2 = rv$getGeneSummary(return_mat = TRUE)
    ht = gcap.plotProfile(data2,
      samples = c("B", "A", "C", "D", "F"),
      genes = rownames(data2)[1:10],
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        cbar = ComplexHeatmap::anno_oncoprint_barplot(),
        foo1 = c("g1", "g1", "g2", "g2", "g3"),
        bar1 = ComplexHeatmap::anno_points(1:5)
      )
    )
    ht
    ComplexHeatmap::draw(ht, merge_legends = TRUE)
  
    data2 = rv$getCytobandSummary(return_mat = TRUE)
    gcap.plotProfile(data2)
  }
  
  expect_type(ht, "S4")
})

