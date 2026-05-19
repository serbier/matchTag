test_that("MTA initializes correctly with SNP parameters", {
  mta <- MTA$new(
    alleleID = "120515446-16-G/T",
    alleleSequence = "TGCAGAGGAAGAATCATGGGTTAGCTTGTTCAAATCTTCGGGTAACTAAACAAATAAAACGACTTGTGC",
    favorableAllele = "T"
  )

  expect_s3_class(mta, "R6")
  expect_equal(mta$alleleID, "120515446-16-G/T")
  expect_equal(mta$alleleSequence, "TGCAGAGGAAGAATCATGGGTTAGCTTGTTCAAATCTTCGGGTAACTAAACAAATAAAACGACTTGTGC")
  expect_equal(mta$favorableAllele, "T")
  expect_equal(mta$mtaType, "SNP")
  expect_equal(mta$alleles, c("G", "T"))
  expect_equal(mta$loc_tag_pos, 16)
})
